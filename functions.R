require(PaleoSpec)
require(zoo)
require(RColorBrewer)
require(ncdf4)
require(maps, include.only = 'map') 
require(polyclip, include.only = 'polysimplify') 
require(stats, include.only = 'filer')
require(utils, include.only='data')

#--------------------------------------------------------------#
#' @title Spectral Gain 
#' @description Computes the spectral gain by dividing the output by the input spectrum, published under 
#' https://github.com/paleovar/TimescaleDependency/blob/main/processing/functions_processing.R and based on functions from the PaleoSpec package
#' @param specList list that contains at least two objects of class "spec"
#' @param input index of the input spectrum 
#' @param output index of the output spectrum
#' @param iRemoveLowest number of lowest frequencies to remove (e.g. to remove detrending bias)
#' @return object of class "spec"
#' @export
transferSpec <- function (specList, input=input.spec, output=output.spec, iRemoveLowest = 1)
  #, weights = rep(1, length(specList)
{
  remove.lowestFreq <- function (spec, iRemove){ 
    {
      if (iRemove == 0) 
        index = seq(spec$spec)
      else index <- (-(1:iRemove))
      spec$spec <- spec$spec[index]
      spec$freq <- spec$freq[index]
      spec$dof <- spec$dof[index]
      return(spec)
    }
  }
  get.fend.existing <- function (x){
    return(max(x$freq[!is.na(x$spec)]))
  }
  get.fstart.existing <- function (x) {
    return(min(x$freq[!is.na(x$spec)]))
  }
  get.df <- function (x){
    return(mean(diff(x$freq)))
  }
  AddConfInterval_Fdist <- function(transferspec, var.dof1, var.dof2, pval = 0.2){ 
    {
      if (!(length(transferspec$spec) == length(var.dof1)) && (!length(transferspec$spec) == 
                                                               length(var.dof2))) {
        stop("same lengths must be provided")
      }
      if (!(is.numeric(transferspec$spec)) || !(is.numeric(var.dof1)) || 
          !is.numeric(var.dof2)) {
        stop("non-numeric arguments")
      }
      res <- matrix(NA, nrow = length(transferspec$spec), ncol = 2)
      for (i in 1:length(transferspec$spec)) {
        QF <- qf(p = c(pval/2, (1 - pval/2)), df1 = var.dof1[i], 
                 df2 = var.dof2[i])
        tmp <- QF * transferspec$spec[i]
        res[i, ] <- tmp
      }
      transferspec$lim.1 <- res[,2]
      transferspec$lim.2 <- res[,1]
      class(transferspec) <- "spec"
      return(transferspec)
    }
  }
  specList <- lapply(specList, remove.lowestFreq, iRemove = iRemoveLowest)
  freqRef <- seq(from = min(unlist(lapply(specList, get.fstart.existing))), 
                 to = max(unlist(lapply(specList, get.fend.existing))), 
                 by = min(unlist(lapply(specList, get.df))))
  specList.interpolated <- list()
  for (i in 1:length(specList)) specList.interpolated[[i]] <- SpecInterpolate(freqRef, 
                                                                              specList[[i]])
  NSpectra <- length(specList.interpolated)
  result <- list(freq = specList.interpolated[[1]]$freq, spec = rep(0, 
                                                                    length(specList.interpolated[[1]]$spec)))
  specMatrix <- matrix(NA, NSpectra, length(specList.interpolated[[1]]$spec))
  dofMatrix <- matrix(NA, NSpectra, length(specList.interpolated[[1]]$spec))
  for (i in 1:length(specList.interpolated)) {
    if (sum((result$freq - specList.interpolated[[i]]$freq)^2) > 
        0.1) 
      stop("Different spectra length or resolutions")
    specMatrix[i, ] <- specList.interpolated[[i]]$spec
    dofMatrix[i, ] <- specList.interpolated[[i]]$dof
  }
  
  var.dof1 <- dofMatrix[output, ]
  var.dof2 <- dofMatrix[input,]
  result$spec <- mapply('/', specMatrix[output, ], specMatrix[input, ]) #na.rm=TRUE
  result <- AddConfInterval_Fdist(result, dofMatrix[output, ], dofMatrix[input,])
  class(result) <- "spec"
  return(list(spec = result))
}

#----------------------------#
#' @title Assign climate zones to latitudes

#' @param d dimensions of the matrix to be assigned to zones
#' @param latitude values of latitude bands
#' @return indices that assign matrix to zones
#' @export
get_zones <- function(d, latitude){
  zones <- list()
  zones$tropics <- c(-23.5, 23.5)
  zones$mid_lats_n <- c(23.5, 66.5)
  zones$mid_lats_s <- rev(-1*zones$mid_lats_n)
  zones$hight_lats_n <- c(66.5, 90)
  zones$hight_lats_s <- rev(-1*zones$hight_lats_n)

  zones_idx <- lapply(zones, function(x) which(latitude >= x[1] & latitude <= x[2]))
  lons <- 1:d[1]

  get_ind_lats <- function(lats, lons=1:d[1]){
    spectra.ind <- c()
    for (lat in lats){
      for (lon in lons){
        spectra.ind_new <- (lat-1)*length(lons)+lon
        spectra.ind <- c(spectra.ind, spectra.ind_new)
      }
    }
    return(spectra.ind)
  }

  zones_ind <- lapply(zones_idx, function(x) get_ind_lats(x))
  return(zones_ind)
}

#----------------------------#
#' @title Timeseries at location 
#' @description Extracts time series at certain location and filters for NA values
#' @param custom_TS object with $data parameter that contains the field of observations
#' @param location vector that contains the location indices
#' @return timeseries
#' @export
to.ts <- function(custom_TS, location)
{
  stopifnot(length(location)==2)
  
  # Approx single na values
  x <- na.approx(custom_TS$data[location[1],location[2],], na.rm=F, maxgap=1)
  x <- try(na.contiguous(x), silent=T)
  
  if(inherits(x, "try-error")) return(rep(0,9))
  
  if(length(x) < 0.5*length(custom_TS$data[location[1],location[2],])) warning("Only 50% of timeseries used!")
  return(ts(x, start=custom_TS$tstart, deltat=custom_TS$deltat/360))
}

#----------------------------#
#' @title Correlation length estimate 
#' @description Get estimate of correlation length defined by 1/e decay using bootstrapping and the autocorrelation function acf
#' @param x vector
#' @return  correlation length
#' @export
get_corr_length <- function(x)
{
  x <- na.contiguous(x)
  l <- length(x)
  boot.l <- 100
  stopifnot(3*boot.l < l)
  
  # Bootstrap 100 correlation lenghts
  i <- sample.int(l-boot.l,1000,replace=T)
  return(sapply(i, function(i.boot) min(which(acf(x[i.boot:(i.boot+boot.l)], plot=F)$acf < 0.3678794))))
}

#----------------------------#
#' @title Color scale for image() function
#' @description This function creates a color scale for use with the image() function. 
#' Input parameters should be consistent with those used in the corresponding image plot. 
#' The "axis.pos" argument defines the side of the axis. The "add.axis" argument defines
#' whether the axis is added (default: TRUE) or not (FALSE).
#' @param z input values for color scale
#' @param zlim limit of color scale
#' @param col Color palette
#' @param breaks breaks of color scale axis
#' @param axis.pos plotting parameter, default 1 
#' @param add.axis plotting parameter, defult TRUE
#' @return 
#' @export
image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, axis.pos=1, add.axis=TRUE, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  if(axis.pos %in% c(1,3)){ylim<-c(0,1); xlim<-range(breaks)}
  if(axis.pos %in% c(2,4)){ylim<-range(breaks); xlim<-c(0,1)}
  plot(1,1,t="n",ylim=ylim, xlim=xlim, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(axis.pos %in% c(1,3)){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(axis.pos %in% c(2,4)){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
  box()
  if(add.axis) {axis(axis.pos)}
}

#----------------------------#
#' @title Polygon from Matrix
#' @description #create polygon from matrix from www.menugget.blogspot.com
#' @param x seq(0,1,,dim(z)[1])
#' @param y seq(0,1,,dim(z)[2])
#' @param z matrix 
#' @param n defaul NULL
#' @return 
#' @export
matrix.poly <- function(x, y, z=mat, n=NULL){
  if(missing(z)) stop("Must define matrix 'z'")
  if(missing(n)) stop("Must define at least 1 grid location 'n'")
  if(missing(x)) x <- seq(0,1,,dim(z)[1])
  if(missing(y)) y <- seq(0,1,,dim(z)[2])
  poly <- vector(mode="list", length(n))
  for(i in seq(length(n))) {
    ROW <- ((n[i]-1) %% dim(z)[1]) +1
    COL <- ((n[i]-1) %/% dim(z)[1]) +1
    
    dist.left <- (x[ROW]-x[ROW-1])/2
    dist.right <- (x[ROW+1]-x[ROW])/2
    if(ROW==1) dist.left <- dist.right
    if(ROW==dim(z)[1]) dist.right <- dist.left
    
    dist.down <- (y[COL]-y[COL-1])/2
    dist.up <- (y[COL+1]-y[COL])/2
    if(COL==1) dist.down <- dist.up
    if(COL==dim(z)[2]) dist.up <- dist.down
    
    xs <- c(x[ROW]-dist.left, x[ROW]-dist.left, x[ROW]+dist.right, x[ROW]+dist.right)
    ys <- c(y[COL]-dist.down, y[COL]+dist.up, y[COL]+dist.up, y[COL]-dist.down)
    poly[[i]] <- data.frame(x=xs, y=ys)
  }
  return(polysimplify(poly))
}

#----------------------------#
#' @title Mean Spectrum
#' @description Wrapper for the PaleoSpec::MeanSpectrum function. We introduced a correction to the number of records, 
#' in case spectra where only partially overlapping.
#' Part of PTBox https://palmodapp.cloud.dkrz.de/index.php/ptbox/
#' @param speclist list of spectra
#' @param iRemoveLowest number of lowest frequencies to remove (e.g. to remove detrending bias)
#' @param weights vector of weights (same length as elements in speclist)
#' @return list(spec,nRecords) spec=average spectrum, nRecords = number of records contributing to each spectral estimate
#' @export
MeanSpec <- function(specList, iRemoveLowest = 1, weights = rep(1, length(specList))){
  meanspec <- PaleoSpec::MeanSpectrum(specList, iRemoveLowest, weights)
  meanspec$spec$spec <- meanspec$spec$spec * length(specList)/meanspec$nRecord
  meanspec$spec <- AddConfInterval(meanspec$spec)
  return(meanspec)
}

#----------------------------#
#' @title Linear detrending
#' 
#' @param y input vector
#' @return detrended outout vector
#' @export
detrend <- function(y){
  fit <- lm(y ~ index(y))
  return(y - index(y)*fit$coefficients[2] - fit$coefficients[1])
}

#----------------------------#
#' @title Transform lons and lats to hadcm3 indices
#' 
#' @description Grids have to be in the right interval, see ancil$grids.*
#'
#' @param lons Longitudes of map for dots to be added in. Have to be in ascending order and must be part of land or sea grid
#' @param lats Latitudes of map for dots to be added in. Have to be in ascending order and must be part of land or sea grid
#' 
#' @return A list of HadCM3 indices
#'
#' @export
transform_grid <- function(lons,lats)
{
  # load("ancil.Rdata") if necessary
  
  stopifnot(sum(lons > 180) == 0)
  stopifnot(!is.na(lons), !is.na(lats))
  
  if(length(lons) > 1)
  {
    if (sort(lons)[2]-sort(lons)[1] == 3.75)
    {
      land <- T
    } else if (sort(lons)[2]-sort(lons)[1] == 1.25) {
      land <- F
    } else {
      print("Unknown Grid.")
      stop()
    }
  } else {
    return( interpolate_grid_to_index(lons,lats) )
  }
  
  if (land) # Land
  {
    l <- length(ancil$grid.land$lon)
    shifted_lons <- (ancil$grid.land$lon-180)[c((l/2+1):l,1:(l/2))]
    shifted_lats <- ancil$grid.land$lat
  } else {
    l <- length(ancil$grid.sea$lon)
    shifted_lons <- (ancil$grid.sea$lon-180)[c((l/2+1):l,1:(l/2))]
    shifted_lats <- ancil$grid.sea$lat
  }
  
  lons <- match(shifted_lons,lons)
  lats <- match(shifted_lats,lats)
  
  return( list(lons=lons[!is.na(lons)], lats=lats[!is.na(lats)]) )
}

#----------------------------#
#' @title Plot windfields
#' 
#' @description Give u/v wind fields to plot vectorplot with strength through color matrix.
#'
#' 
#' @param fld_u u wind field
#' @param fld_v v wind field
#' @param zlim Limits for wind speed
#' @param col color palette for wind speeds
#' @param stretch stretch factor for vectors
#' @param strength.cutoff cutoff value below which vectors are replaced by dots
#' @param arrowhead.size size multiplier for arrowheads
#' @param by only print by-th vector instead of every vector
#' @export
#'
#' @examples
plot_windfield <- function(fld_u, fld_v, zlim, col, stretch, strength.cutoff, arrowhead.size=1, by=2)
{
  # load("ancil.Rdata") if necessary
  strength_raw <- sqrt(fld_u*fld_u+fld_v*fld_v)
  
  d <- dim(fld_u)
  
  lons <- ancil$grid.land$lon-180
  lats <- rev(ancil$grid.land$lat)
  grid <- transform_grid(lons, lats)
  grid$lats <- grid$lats[-1]
  lats <- lats[-73]+1.25
  lons <- lons
  
  fld_u <- array(fld_u[grid$lons,grid$lats], d[1]*d[2])
  fld_v <- array(fld_v[grid$lons,grid$lats], d[1]*d[2])
  strength <- fld_u*fld_u+fld_v*fld_v
  
  stopifnot(by > 1)
  mask.number <- rep(F,by*(by+by%%2)*96)
  for (i in 1:(by*(by+by%%2)))
  {
    i.ceil <- ceiling(i/by)
    if (i %% by == 0) 
    {
      offset <- rep(F,(i.ceil%%2)*round(0.5*by))
      mask.number[((i-1)*96+1):(i*96)] <- c(offset,array(c(T,rep(F,by-1)),96))[1:96]
    }
  }
  mask.number <- array(mask.number, d[1]*d[2])
  mask.na <- (!is.na(fld_u) & !is.na(fld_v) & 
              fld_v != Inf & fld_u != Inf &
              fld_v != -Inf & fld_u != -Inf)
  mask.strength <- strength > strength.cutoff & !is.na(strength) & strength != Inf & strength != -Inf
  mask <- mask.na & mask.strength & mask.number
  
  strength_raw <- array(strength_raw, c(96,73))
  strength_raw[,73] <- strength_raw[,72]
  plot_matrix(array(strength_raw, c(96,73)),
              zlim=zlim,
              col=col)
  u <- rep(lons,d[2])[mask.number & !mask.strength & mask.na]
  v <- rep(lats,each=d[1])[mask.number & !mask.strength & mask.na]
  points(u,v,cex=0.5,pch=20)
  u <- rep(lons,d[2])[mask]
  v <- rep(lats,each=d[1])[mask]
  uu <- u+stretch*fld_u[mask]
  vv <- v+stretch*fld_v[mask]
  arrows(x0=u, x1=uu,
         y0=v, y1=vv
         ,length=0.03*arrowhead.size,lwd=2)
}

#----------------------------#
#' @title Add dots to plot_matrix
#' 
#' @description Lons and Lats have to be in the right interval, check ancil$grid.land/ancil$grid.sea
#'
#' @param cond 2D boolean array where dots should be added
#' @param lons Longitudes of map for dots to be added in. Have to be in ascending order and must be part of land or sea grid
#' @param lats Latitudes of map for dots to be added in. Have to be in ascending order and must be part of land or sea grid
#' @param col Color of the dots
#'
#' @export
add_dots<-function(cond, lons, lats, col=rgb(0,0,0,alpha=7,maxColorValue = 10))
{
  # load("ancil.Rdata") if necessary
  d <- dim(cond)
  stopifnot(length(dim(cond)) == 2)
  
  if (missing(lons) | missing(lats))
  {
    if(d[1] == 96)
    {
      lons <- ancil$grid.land$lon-180
      lats <- rev(ancil$grid.land$lat)
    } else if (d[1] == 288) {
      lons <- ancil$grid.sea$lon-180
      lats <- ancil$grid.sea$lat
    } else {
      lons <- NA
      lats <- NA
    }
  }
  
  ind <- which(cond, arr.ind=T)
  
  grid <- transform_grid(lons, lats)
  
  lines(lons[grid$lons[ind[,1]]], lats[grid$lats[ind[,2]]], type="p", pch=19, cex=0.1, c=col)
}

#----------------------------#
#' @title Plot line with shaded area used for conf intervals
#'
#' @param x x-axis values
#' @param y mean values
#' @param upper upper bound for shading
#' @param lower lower bound for shading
#' @param col color
#' @param lty line type
#' 
#' @export
plot_mean_conf <- function(x, y, upper, lower, col, lty)
{
  polygon(c(x, rev(x))
          , c(lower, rev(upper))
          , col=rgb(t(col2rgb(col))
                    , alpha=100
                    , maxColorValue=255)
          , border = NA)
  lines(x
        , y
        , type='l'
        , lw=2
        , lty=lty
        , col=rgb(t(col2rgb(col))
                  , alpha=255
                  , maxColorValue=255))
}

#----------------------------#
#' @title  Plot a matrix on a grid
#' 
#' @description Lons and Lats have to be in the right interval, check ancil$grid.land/ancil$grid.sea
#' 
#' @param fld 2D matrix to plot
#' @param zlim Array of length 2 for range for colorpalette
#' @param lons Longitudes of map for dots to be added in. Have to be in ascending order and must be part of land or sea grid
#' @param lats Latitudes of map for dots to be added in. Have to be in ascending order and must be part of land or sea grid
#' @param col Colorpalette with breaks
#' @param NA_col Color of NA values, leave out for transparent NA values
#' @param add Add to plot or write to new plot
#' @param ... Additional arguments to pass to image()
#'
#' @import RColorBrewer
#'
#' @export
plot_matrix <- function (fld, zlim, lons, lats, col = colorRampPalette(rev(brewer.pal(nbreakpoints - 1, "YlGnBu"))), 
                         NA_col = NULL, add, ...)
{
  # load("ancil.Rdata") if necessary
  d <- dim(fld)
  if (missing(lons) | missing(lats)) {
    if (d[1] == 96) {
      lons <- ancil$grid.land$lon - 180
      lats <- rev(ancil$grid.land$lat)
    }
    else if (d[1] == 288) {
      lons <- ancil$grid.sea$lon - 180
      lats <- ancil$grid.sea$lat
    }
    else {
      lons <- NA
      lats <- NA
    }
  }
  grid <- transform_grid(lons, lats)
  if(missing(add)) add <- !is.null(NA_col)
  if (!is.null(NA_col)) 
    image(lons, lats, matrix(1, d[1], d[2]), col = NA_col, 
          xlab = "", ylab = "", axes = FALSE, add=add, ...)
  image(lons, lats, fld[grid$lons, grid$lats], zlim = zlim, 
        col = col, add = add, breaks = seq(zlim[1], 
                                           zlim[2], length.out = length(col) + 1), xlab = "", 
        ylab = "", axes = FALSE, ...)
}

#----------------------------#
#' @title  Plot lettering inside plot
#' 
#' @description Note: This function will reset a plot, only use at the very end
#' 
#' @param letter String to paste into top-left corner
#' @param cex cex param to set for letter
#'
#' @export
plot_letter <- function (letter, cex) 
{
  cex.old <- par("cex")
  if(missing(cex)) cex <- cex.old
  par(new = T, cex=cex)
  
  plot(1, axes = F, type = "n", xlab = "", ylab = "")
  box <- par("usr")
  fin <- par("fin")
  width <- box[2] - box[1]
  height <- box[4] - box[3]
  
  offset.y <- 0.15
  offset.x <- 0.04
  whitespace.multiplier <- c(1,1)[(strwidth(letter) < 0.10)+1]
  whitespace.left <- 0.2*whitespace.multiplier
  whitespace.right <- 1.3*whitespace.multiplier
  whitespace.topbottom <- 0.7*whitespace.multiplier
  
  rect(box[1] + offset.x * width/fin[1] - whitespace.left * strwidth(letter), 
       box[4] - offset.y * height/fin[2] + whitespace.topbottom * strheight(letter), 
       box[1] + offset.x * width/fin[1] + whitespace.right * strwidth(letter), 
       box[4] - offset.y * height/fin[2] - whitespace.topbottom * strheight(letter), 
       col = "white", lty = 0)
  text(x = box[1] + offset.x * width/fin[1], y = box[4] - offset.y * 
         height/fin[2], letter, adj = 0, font = 2)
  par(cex=cex.old)
}

#----------------------------#
#' @title Add colorbar as standalone plot
#'
#' @param zlim Array with limits of colorbar
#' @param col Colorpalette with breaks
#' @param add.axis Boolean whether or not axis should be added
#' @param axis.pos Position of axis to add
#' @param cex cex of axis
#' @param ... Additional parameters to pass to plot
#'
#' @export
plot_colorbar <- function(zlim, col = heat.colors(12), add.axis=TRUE, axis.pos=1, cex=1, ...)
{
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  if(axis.pos %in% c(1,3)){ylim<-c(0,1); xlim<-range(breaks)}
  if(axis.pos %in% c(2,4)){ylim<-range(breaks); xlim<-c(0,1)}
  plot(1,1,t="n",ylim=ylim, xlim=xlim, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(axis.pos %in% c(1,3)){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(axis.pos %in% c(2,4)){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
  box()
  if(add.axis) {
    axis(axis.pos,labels=F)
    mtext(axTicks(axis.pos),side=axis.pos, at=axTicks(axis.pos),font=2,line=2,cex=cex)
  }
}

#----------------------------#
#' @title Plot Axis
#' 
#' @description Add axis to plot using mtext, with option to set line and adjust
#' 
#' @param side Side of plot
#' @param at Position of ticks
#' @param ticks Label of ticks
#' @param label Label of axis
#' @param lineVal line for ticks
#' @param lineLab line for label of axis
#' @param adj adj of label of axis
#' @param font font of text
#' @param col color of axis and ticks
#' 
#' @export
plot_axis <- function (side, at = axTicks(side), ticks = axTicks(side), label = "", 
                       lineVal = 0.5, lineLab = 1.5, adj = 0.5, font = 1, col = "black") 
{
  if(par("pin")[(side+1)%%2+1] < sum(strwidth(ticks, "inches")))
  {
    if(length(ticks)%%2==1) i <- array(c(F,T),length(ticks))
    else i <- array(c(T,F),length(ticks))
    ticks <- ticks[i]
    at <- at[i]
  }
  axis(side, at = at, labels = F, col = col)
  mtext(side = side, at = at, ticks, line = lineVal, font = font, 
        col = col)
  mtext(side = side, label, line = lineLab, adj = adj, font = font)
}

#----------------------------#
#' @title Generate custom timeseries list from nc file
#' 
#' @param nc An `nc_open` nc file
#' @param varname The variable name as a string
#'
#' @export
nc_var_to_TS <- function(nc, varname)
{
  x <- list()
  x$data <- ncvar_get(nc, varname)
  x$time <- ncvar_get(nc, "t")
  x$tstart <- x$time[1]
  deltat <- x$time[2]-x$time[1]
  
  if(x$time[length(x$time)] == x$tstart + (length(x$time)-1)*deltat)
  {
    x$deltat <- deltat
  } else {
    warning(paste(nc$filename,"has irregular timesteps, not saving deltat"))
  }
  
  return(x)
}

#----------------------------#
#' @title Find array indices from hadcm land array corresponding to supplied lon/lat
#'
#' @param lon Longitude
#' @param lat Latitude
#' 
#' @return Array indices for HadCM3 land grid corresponding to lon/lat. Can 
#'
#' @examples
#' interpolate_grid_to_index(-42.32,75.1)
#'
#' @export
interpolate_grid_to_index <- function(lon,lat)
{
  if (lon > 180) lon <- lon-360
  lats <- seq(90,-90,-2.5)
  lons <- c(seq(0,180,3.75),seq(-176.25,-3.75,3.75))
  
  return( c(which.min(abs(lons - lon)),which.min(abs(lats - lat))) )
}

#----------------------------#
#' @title Variance computation by integration of the spectrum
#' @description Simple computation of variance on frequency intervals by integration of the spectrum, akin to PaleoSpec::GetVarFromSpectra(), 
#' with consideration of degrees of freedom
#' @param specs list of objects of class "spec"
#' @param tsc.in target timescales
#' @return  returns list with variance and degrees of freedom
#' @export
getvaranddof_from_spec <- function(specs, tsc.in){
  out <- lapply(specs, function(x) PaleoSpec::GetVarFromSpectra(x, c(1/tsc.in[[2]], 1/tsc.in[[1]])))
  
  var.out <- sapply(out,function(x){x$var})
  dof.out <- sapply(out,function(x){x$dof})
  
  return(list(var=as.numeric(var.out), dof=as.numeric(dof.out)))
}

#----------------------------#
#' @title Wrapper to prepare list of spectra for variance estimate
#' @description relies on PaleoSpec::MakeEquidistant() and PaleoSpec::SpecMTM()
#' @param tslist  list with timeseries / objects / vectors
#' @param DT.res list with target resolutions
#' @return list with objects of class "spec"
#' @export
get_spectra_for_var <- function(tslist,DT.res)
{
  tslist.int<-list()
  for (m in 1:length(tslist)){
    if (length(tslist[[m]])<=1) {tslist.int[[m]]<-NA} else {
      tslist.int[[m]] <- MakeEquidistant(index(tslist[[m]]),coredata(tslist[[m]]),dt=DT.res[m])
    }
  }
  names(tslist.int) <- names(tslist)
  specs <- lapply(tslist.int, function(x) SpecMTM(na.omit(x)))
  
  return(specs)
}

#----------------------------#
#' @title Weighted mean
#' @param x input values
#' @param weights assigned weights
#' @param na.rm default FALS, otherwise NAs to be ignored in mean
#' @return value
#' @export
weighted_mean <- function(x, weights, na.rm=FALSE)
{
  if (na.rm) {
    weights <- weights[i <- !is.na(x)]
    x <- x[i]
  }
  return(sum(x*weights)/sum(weights))
}

#----------------------------#
#' @title Get nearest index
#' 
#' @param arr An array where x should be searched
#' @param x Value whose nearest index in arr should be searched
#'
#' @export
getNearestIndex <- function(arr, x)
{ 
  return(which.min(abs(arr - x)))
}

#----------------------------#
#' @title Generate polygons to plot where matrix of global values fulfills condition
#' 
#' @description Creates a list of polygons to plot with polygon(). 
#' This process is very computing-intensive. 
#' It is advised to save the 
#' resulting polygons.
#'
#' @param cond A boolean array condition
#' @param lons Longitude of grid
#' @param lats Latitudes of grid
#' 
#' @return A simplified list of polygons to plot with polygon()
#'
#' @examples
#' arr <- array(c(1,0),c(6,5)) >= 0.5
#' generate_polygon_area(arr, lons=seq(-7.5,7.5, 3.75), lats=seq(-5,5,2.5))
#'
#' @importFrom polyclip polysimplify
#'
#' @export
generate_polygon_area <- function(cond, lons, lats)
{
  # load("ancil.Rdata") if necessary
  d <- dim(cond)
  stopifnot(length(d) == 2)
  
  if (missing(lons) | missing(lats))
  {
    if(d[1] == 96)
    {
      lons <- ancil$grid.land$lon-180
      lats <- rev(ancil$grid.land$lat)
    } else if (d[1] == 288) {
      lons <- ancil$grid.sea$lon-180
      lats <- ancil$grid.sea$lat
    } else {
      lons <- NA
      lats <- NA
    }
  }
  
  grid <- transform_grid(lons, lats)
  
  ind.cond <- which(cond[grid$lons, grid$lats])
  
  poly <- vector(mode="list", length(ind.cond))
  for(i in seq(length(ind.cond))) 
  {
    ROW <- ((ind.cond[i]-1) %% d[1]) + 1
    COL <- ((ind.cond[i]-1) %/% d[1]) + 1
    
    dist.left <- (lons[ROW]-lons[ROW-1])/2
    dist.right <- (lons[ROW+1]-lons[ROW])/2
    if(ROW==1) dist.left <- dist.right
    if(ROW==d[1]) dist.right <- dist.left
    
    dist.down <- (lats[COL]-lats[COL-1])/2
    dist.up <- (lats[COL+1]-lats[COL])/2
    if(COL==1) dist.down <- dist.up
    if(COL==d[2]) dist.up <- dist.down
    
    xs <- c(lons[ROW]-dist.left, lons[ROW]-dist.left, lons[ROW]+dist.right, lons[ROW]+dist.right)
    ys <- c(lats[COL]-dist.down, lats[COL]+dist.up, lats[COL]+dist.up, lats[COL]-dist.down)
    poly[[i]] <- data.frame(x=xs, y=ys)
  }
  return(polysimplify(poly))
}

#----------------------------#
#' @title Compute area-weighted global mean of sea or land arrays
#'
#' @description If a 3D array is passed, fldmean will assume the third dimension as time and return a timeseries of fldmeans.
#'
#' @param fld 2/3D array of land or sea values
#' @param na.rm A boolean to define whether NA's should be ignored
#' @param mask 2D index mask
#' 
#' @return A single area-weighted global mean or a timeseries of area-weighted global means.
#'
#' @importFrom utils data
#'
#' @export
fldmean <- function(fld, na.rm = F, mask = array(1, dim(fld)[1:2]))
{
  # load("ancil.Rdata") if necessary
  
  d <- dim(fld)
  
  # If we don't ignore NA's, stop if there are NA's
  if(!na.rm) stopifnot(sum(is.na(fld[mask])) == 0)
  
  if(length(d) == 3) #2D timeseries
  {
    if (d[1] == 288 & d[2] == 144) return(apply(fld, 3, .fldmean, sea=T, mask))
    if (d[1] == 96 & d[2] == 73) return(apply(fld, 3, .fldmean, sea=F, mask))
  }
  
  if(length(d) == 2) #2D field
  {           
    if (d[1] == 288 & d[2] == 144) return(.fldmean(fld,T,mask))
    if (d[1] == 96 & d[2] == 73) return(.fldmean(fld,F,mask))
  }
  
  stop()
}

.fldmean <- function(fld, sea=F, mask)
{
  # load("ancil.Rdata") if necessary
  mask <- which(!is.na(fld) & mask, arr.ind=T)
  
  if(sea)
  {
    return(sum(fld[mask]*ancil$weights.sea[mask])/sum(ancil$weights.sea[mask]))
  } else {
    return(sum(fld[mask]*ancil$weights.land[mask])/sum(ancil$weights.land[mask]))
  }
}

#----------------------------#
#' @title Fix missing months in HadCM3 simulations
#'
#' @description Attempts to fill skipped timesteps in monthly .nc files using the seasonality of the previous three years. Resulting .nc filed will be named \code{folders}*.nc.fixed, a log is written into \code{folders}fix.log.
#'
#' @param folders Path to the folder containing the .nc files
#' @param pattern [optional] Exact name of the .nc file to be fixed, if left empty, every .nc file not yet processed will be fixed.
#'
#' @examples
#' \donttest{
#' fix_monthly_data("~/data/2dst/","xmzke.nc")
#' }
#' 
#' @import ncdf4
#' @importFrom stats filter
#' 
#' @export
fix_monthly_data <- function(folders, pattern=".nc$")
{
  path = c()
  for (folder in folders)
  {
    path <- c(path,paste(folder,list.files(path=folder, pattern=pattern),sep=""))
  }
  
  time <- list()
  for (i in 1:length(path))
  {
    if(file.exists(paste(path[i],".fixed",sep=""))) next
    
    nc <- nc_open(path[i])
    
    # Find skips
    time <- ncvar_get(nc, "t")
    dt <- filter(time, c(1.,-1.))
    
    # Remove repeated timesteps
    rem <- dt==0 & !is.na(dt)
    if(sum(rem)>0)
    {
      ind.rem <- which(rem)
      
      system(paste("cdo delete,timestep=",paste(ind.rem,collapse=",")," ",path[i]," ",path[i],".new",sep=""))
      
      nc_close(nc)
      system(paste("mv ",path[i],".new ",path[i],sep=""))
      nc <- nc_open(path[i])
      
      # Find skips again
      time <- ncvar_get(nc, "t")
      dt <- filter(time, c(1.,-1.))
    }
    
    skips <- dt!=30 & !is.na(dt)
    ind.skips <- which(skips)
    len.skips <- dt[skips]/30-1
    
    if(file.exists(paste(path[i],".fixed",sep=""))) next
    
    if(sum(skips) == 0) {
      system(paste("cp ",path[i]," ",path[i],".fixed",sep=""))
      next
    }
      
    # Generte fill .nc to write data that will be mergetimed later
    system(paste("cdo -b f64 seltimestep,1/",sum(len.skips)," ",
                 path[i]," ",path[i],".fill",sep=""))
    
    # open fill file
    nc_fill <- nc_open(paste(path[i],".fill",sep=""), write=T)
    
    # Fill time
    fill <- array(NA, sum(len.skips))
    for (j in 1:sum(skips))
    {
      from <- sum(len.skips[1:j])-len.skips[j]+1
      fill[from:(from+len.skips[j]-1)] <- time[ind.skips[j]] + 30*(1:len.skips[j])
    }
    ncvar_put(nc_fill, "t", fill)
    
    for (var in c(names(nc_fill$var)))
    {
      print(paste("Filling var ",var,sep=""))
      ts <- ncvar_get(nc, var)
      n_dim <- length(dim(ts))
      n_timesteps <- dim(ts)[n_dim]
      # Set dimensions of fill array
      if (n_dim == 1) 
      {
        fill <- array(NA, sum(len.skips))
      } else {
        fill <- array(NA, c(dim(ts)[1:(n_dim-1)],sum(len.skips)))
      }
      
      for (j in 1:sum(skips))
      {
        # Assertions
        if (ind.skips[j] < 3*12)
        {
          print(paste("Too early skip found in var",path[i],"at index",ind.skips[j],"(filling NA)"))
          if (var == c(names(nc_fill$var))[1])
          {
            system(paste("echo \"",path[i]," NA skip @ ",ind.skips[j]," for ",len.skips[j],"\" >> ",folders,"fix.log",sep=""))
          }
          next
        }
        if (len.skips[j] > 11)
        {
          print(paste("Too big skip found in",path[i],"at index",ind.skips[j],"(filling NA)"))
          if (var == c(names(nc_fill$var))[1])
          {
            system(paste("echo \"",path[i]," NA skip @ ",ind.skips[j]," for ",len.skips[j],"\" >> ",folders,"fix.log",sep=""))
          }
          next
        }
        
        if (var == c(names(nc_fill$var))[1])
        {
          system(paste("echo \"",path[i]," filled skip @ ",ind.skips[j]," for ",len.skips[j],"\" >> ",folders,"fix.log",sep=""))
        }
        
        from <- sum(len.skips[1:j])-len.skips[j]+1
        # Extrapolate last value before skip using a seasonality "multiplier" averaged from three years prior
        if (n_dim == 1)
        {
          seasonality <- apply(array(ts[ind.skips[j]+(-36:-1)], c(12,3)), 1, mean)
          seasonality <- (seasonality/seasonality[1])[2:12]
          seasonality[is.nan(seasonality)] <- 0.0
          fill[from:(from+len.skips[j]-1)] <- seasonality[1:len.skips[j]]*ts[ind.skips[j]]
        } else {
          # Get seasonality for generalized dimensions via stacking of index matrices
          ts.dim <- dim(ts)[1:(n_dim-1)]
          ts.index <- array(c(rep(F,prod(ts.dim)*ind.skips[j]-37),
                              rep(T,prod(ts.dim)*36),
                              rep(F,prod(ts.dim)*(n_timesteps-ind.skips[j]+1))), 
                            c(ts.dim,n_timesteps))
          # Average over last dimension (the three prior years that were stacked)
          seasonality <- apply(array(ts[ts.index], c(ts.dim,12,3)), 1:n_dim, mean)
          # Normalize seasonality to values at first month
          seasonality.index <- array(c(rep(T,prod(ts.dim)),
                                       rep(F,prod(ts.dim)*11)), 
                                     c(ts.dim,12))
          seasonality <- (seasonality/array(seasonality[seasonality.index],c(ts.dim,12)))[!seasonality.index]
          seasonality[is.nan(seasonality)] <- 0.0
          
          # Copy over 
          ts.index <- array(c(rep(F,prod(ts.dim)*(ind.skips[j]-1)),
                              rep(T,prod(ts.dim)),
                              rep(F,prod(ts.dim)*(n_timesteps-ind.skips[j]))), 
                            c(ts.dim,n_timesteps))
          fill.index <- array(c(rep(F,prod(ts.dim)*(from-1)),
                                rep(T,prod(ts.dim)*len.skips[j]),
                                rep(F,prod(ts.dim)*(sum(len.skips)-from-len.skips[j]+1))), 
                              c(ts.dim,sum(len.skips)))
          fill[fill.index] <- array(ts[ts.index], c(ts.dim,len.skips[j])) * array(seasonality, c(ts.dim,len.skips[j]))
        }
      }
      ncvar_put(nc_fill, var, fill)
    }
    nc_close(nc_fill)
    nc_close(nc)
  
    if(file.exists(paste(path[i],".fixed",sep="")))
    {
      system(paste("rm ",path[i],".fixed",sep=""))
    }
    system(paste("cdo -b f64 -mergetime ",path[i],"* ",path[i],".fixed",sep=""))
    system(paste("rm ",path[i],".fill",sep=""))
  }
}

#----------------------------#
#' @title Convert HadCM3 AD to Days Since
#'
#' @param x AD to convert
#'
#' @return double representing DaysSince
#' 
#' @export
ADToDaysSince <- function(x) { return((x-200)*360) }


#----------------------------#
#' @title Convert HadCM3 Days Since to AD
#'
#' @param x DaysSince to convert
#'
#' @return double representing AD
#' 
#' @export
#' 
DaysSinceToAD <- function(x) { return(x/360+200) }

#----------------------------#
#' @title Draw hatched polygons generated by generate_polygon_area
#'  
#' @param polys Polygons generated by generate_polygon_area()
#' @param double_hatched Boolean whether or not polygons should by hatched once or twice
#' @param col Color of the polygon
#'
#' @export
add_polygons <- function(polys, double_hatched=F, col=rgb(0,0,0,alpha=3,maxColorValue = 10))
{
  for (i in seq(polys))
  {
    polygon(polys[[i]], density=20, angle=45, lwd=1, col=col)
    if (double_hatched) polygon(polys[[i]], density=20, angle=-45, lwd=1, col=col)
  }
}

#----------------------------#
#' @title Add map contours to plot
#'
#' @description Lons and Lats have to be in the right interval, check ancil$grid.land/ancil$grid.sea
#' 
#' @param lons Longitudes of map for dots to be added in. Have to be in ascending order and must be part of land or sea grid
#' @param lats Latitudes of map for dots to be added in. Have to be in ascending order and must be part of land or sea grid
#' @param lgm Boolean, if LGM map should be added (or PI)
#' @param land Boolean, if LGM land resolution land-sea mask (or ocean resolution) should be added
#' @param map Boolean, if Land Sea Mask or modern continent contours should be used
#' @param col Color of the map contours
#'
#' @importFrom maps map
#'
#' @export
add_map <- function(lons, lats, lgm=F, land=T, map="model", col="black")
{
  #data(ancil)
  
  if(map == "model") 
  {  
    if (missing(lons) | missing(lats))
    {
      lons <- ancil$grid.land$lon-180
      lats <- rev(ancil$grid.land$lat)
    }
    grid <- transform_grid(lons, lats)
    
    if(lgm) state <- "LGM_grid-"
    else state <- "PI_grid-"
    if(land) state <- paste(state,"land",sep="")
    else state <- paste(state,"ocean",sep="")
    
    contour(lons, lats, ancil$lsm[[state]][grid$lons,grid$lats], levels=c(0.5), 
            method="simple", c=col, drawlabels=F, lwd=1, add=T, xlab="", ylab="", axes=FALSE)
  }
  if(map == "world") map("world",interior=F,add=T)
}
