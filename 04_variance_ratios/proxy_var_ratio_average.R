#THIS SCRIPT WILL BE CALLED BY F4.R

tsc.in<-list("2-5" = c(2,5),
             "5-50" = c(5,50),
             "50-200" = c(50,200),
             "200-500" = c(200,500))


if(!exists("interpolate")){
  interpolate <- T
}
if(!exists("average")){
  average <- T
}

get_hadcm_bilint_ts <- function(run, list.locs)
{
  #' get_hadcm_ts
  #'
  #' @param run
  #' @param locs 
  #'
  #' @return list of time series for run at nearest gridbox for list of locs
  #' @export
  tser <- list()
  
  l.locs <- lapply(list.locs, function(x){
    if(x[[1]] < 0){x[[1]]<- x[[1]]+360}
    return(x)} 
  )
  loc.names <- names(l.locs)
  
  # Load stationary (linearly detrended + deseasonalized w.r.t. whole TS) TAS
  nc <- nc_open(paste(data.dir, "surface_temperature/stationary/",run,".nc",sep=""))
  temp.2d <- ncdf4::ncvar_get(nc, "temp_1")
  tstart <- nc_var_to_TS(nc, "temp_1")$tstart
  lon <- ncdf4::ncvar_get(nc, "longitude")
  lat <- ncdf4::ncvar_get(nc, "latitude")
  nc_close(nc)
  
  # For every location...
  for (i.loc in 1:length(loc.names))
  {
    print(i.loc)
    y <- ts(
        apply(temp.2d, 3, function(data_tmp){
          interp.surface(list(x=lon, y=lat, z=data_tmp), data.frame(x=l.locs[[i.loc]][1], y=l.locs[[i.loc]][2]))}),
      start=tstart/360, deltat=1/12)
    
    if (!is.null(y)) tser[[loc.names[i.loc]]] <- y 
    else tser[[loc.names[i.loc]]] <- NA
    
  }
  return(tser)
}

ind.ok <- lapply(QC.all, function(state) {unique(unlist(sapply(state, function(x) x$ind.ok)))})
locs <- lapply(ind.ok, function(state) {lapply(state, function(i) c(meta$Lon[i], meta$Lat[i]))})
for (state in c("LGM","PI")) names(locs[[state]]) <- meta$Name[ind.ok[[state]]]

generate <- F
if(generate){
  #1.2 GET HadCM3 @ PPROXY LOCATION USING "LOCS"
  hadcm.ts <- list()
  for (run_type in c("LGM","LGM*","PI","PI*"))
  {
    print(run_type)
    for (run in names(run_dict)[sapply(run_dict, function(x) {x == run_type})])
    {
      print(run)
      hadcm.ts[[run]] <- get_hadcm_bilint_ts(run, locs[[sub('\\*','',run_type)]])
      stopifnot(any(names(hadcm.ts[[run]]) == names(locs[[sub('\\*','',run_type)]])))
    }
  }
  save("hadcm.ts", file="04_variance_ratios/hadcm_bilint_ts.RData")
} else {
  load("04_variance_ratios/hadcm_bilint_ts.RData")
}

#stop()
#2. GET VARRATIONS
#2.1 GET VARRATIOS AT LOCS THAT ARE PRESENT IN PROXY AND MODEL
vars.proxy <- list()
vars.hadcm <- list()
var.ratio <- list()
var.ratio.CI <- list()
spec.hadcm <- list()
meanspec.hadcm<- list()
spec.proxy <- list()

for (state in c("PI","PI*","LGM","LGM*"))
  {
    print(state)
    vars.proxy[[state]] <- list()
    vars.hadcm[[state]] <- list()
    var.ratio[[state]] <- list()
    var.ratio.CI[[state]] <- list()
    spec.hadcm[[state]] <- list()
    meanspec.hadcm[[state]] <- list()
    spec.proxy[[state]] <- list()
  
    state.uf <- sub('\\*','',state)
    
    runs <- names(run_dict)[sapply(run_dict, function(x) {x == state})]
    
    for (tsc in names(tsc.in))
    {
      print(tsc)
      meanspec.hadcm[[state]][[tsc]] <- list()
      if(interpolate && tsc == tail(names(tsc.in),1) && state == "PI"){runs <- runs[which(runs!="xnagh")]}
      if(interpolate && tsc == tail(names(tsc.in),1) && state == "LGM"){runs <- runs[which(runs!="xmzki")]}
      ind.ok <- QC.all[[state.uf]][[tsc]]$ind.ok
      if(length(ind.ok) == 0) next
      
      proxy.ts <- lapply(QC.all[[state.uf]][[tsc]][["WBT.out"]][ind.ok],function(x){x$tsplit[[1]]})
      DT.proxy <- sapply(proxy.ts,function(x){mean(diff(time(x)))})
      locs <- QC.all[[state.uf]][[tsc]]$metafilt$Name
      
      if(average){
        var.ratio[[state]][[tsc]] <- array(NA,length(locs))
        var.ratio.CI[[state]][[tsc]] <- matrix(NA,nrow=length(locs),ncol=2)
            
        DT.hadcm <- sapply(hadcm.ts[[runs[[1]]]][locs],function(x){if(!is.null(x)) mean(diff(time(x))) else NA})
        stopifnot(length(DT.hadcm) == length(DT.proxy))
        DT <- apply(simplify2array(list(DT.hadcm,DT.proxy)),1,max,na.rm=T)
        spec.proxy[[state]][[tsc]] <- get_spectra_for_var(proxy.ts,DT)
        vars.proxy[[state]][[tsc]] <- getvaranddof_from_spec(spec.proxy[[state]][[tsc]], tsc.in[[tsc]])
          
        for(i.run in 1:length(runs)){
            if(!interpolate){DT <- rep(1, length(locs))}
            spec.hadcm[[state]][[tsc]][[i.run]] <- get_spectra_for_var(hadcm.ts[[runs[[i.run]]]][locs],DT)
        }
     
        for(i in 1:length(locs)){
          l <- list()
          for(i.run in 1:length(runs)){
            l <- c(l, spec.hadcm[[state]][[tsc]][[i.run]][i])
          }
        meanspec.hadcm[[state]][[tsc]][[unique(names(l))]] <- MeanSpec(l, iRemoveLowest = 0)$spec
        rm(l)
        }
          
        vars.hadcm[[state]][[tsc]] <- getvaranddof_from_spec(meanspec.hadcm[[state]][[tsc]],tsc.in[[tsc]])
        var.ratio[[state]][[tsc]][1:length(locs)] <- vars.hadcm[[state]][[tsc]]$var/vars.proxy[[state]][[tsc]]$var
          
        for (i in 1:length(locs)){
        var.ratio.CI[[state]][[tsc]][i,] <- PaleoSpec:::ConfRatio(as.numeric(var.ratio[[state]][[tsc]])[i],
                                                                      vars.hadcm[[state]][[tsc]]$dof[i],
                                                                      vars.proxy[[state]][[tsc]]$dof[i],
                                                                      pval=0.1)
        }
        rm(DT.hadcm,DT.proxy)
      
        }else{
        var.ratio[[state]][[tsc]] <- array(NA,length(locs)*length(runs))
        var.ratio.CI[[state]][[tsc]] <- matrix(NA,nrow=length(locs)*length(runs),ncol=2)
        proxy.dof <- array(NA,length(locs)*length(runs))
        hadcm.dof <- array(NA,length(locs)*length(runs))
        
        for (i.run in 1:length(runs)){
          DT.hadcm <- sapply(hadcm.ts[[runs[i.run]]][locs],function(x){if(!is.null(x)) mean(diff(time(x))) else NA})
          stopifnot(length(DT.hadcm) == length(DT.proxy))
          DT <- apply(simplify2array(list(DT.hadcm,DT.proxy)),1,max,na.rm=T)
          
          spec.proxy[[state]][[tsc]] <- get_spectra_for_var(proxy.ts,DT)
          vars.proxy[[state]][[tsc]] <- getvaranddof_from_spec(spec.proxy[[state]][[tsc]], tsc.in[[tsc]])
          if(!interpolate){DT <- rep(1, length(locs))}
          spec.hadcm[[state]][[tsc]][[i.run]] <- get_spectra_for_var(hadcm.ts[[runs[[i.run]]]][locs],DT)
          vars.hadcm[[state]][[tsc]][[i.run]] <- getvaranddof_from_spec(spec.hadcm[[state]][[tsc]][[i.run]],tsc.in[[tsc]])
          ind.ratio.run <- 1:length(locs)+(i.run-1)*length(locs)
          var.ratio[[state]][[tsc]][ind.ratio.run] <- vars.proxy[[state]][[tsc]]$var/vars.hadcm[[state]][[tsc]][[i.run]]$var
          proxy.dof[ind.ratio.run] <- vars.proxy[[state]][[tsc]]$dof
          hadcm.dof[ind.ratio.run] <- vars.hadcm[[state]][[tsc]][[i.run]]$dof
          for (i in 1:(length(locs)*length(runs))){
            i <- (i-1)%%length(locs)+1
            var.ratio.CI[[state]][[tsc]][i+(i.run-1)*length(locs),] <- PaleoSpec:::ConfRatio(as.numeric(var.ratio[[state]][[tsc]])[i],
                                                                                               vars.hadcm[[state]][[tsc]][[i.run]]$dof[i],
                                                                                               vars.proxy[[state]][[tsc]]$dof[i],
                                                                                               pval=0.1)
            }
          }
          rm(ind.ratio.run,DT.hadcm,DT.proxy)
        }
    }
    rm(proxy.ts, locs)
  }
  
rm(spec.hadcm)
gc()

#---STATISTICS---#
hadcm3nc <- nc_open("hadcm_weights.nc")
cell_weights <- ncvar_get(hadcm3nc, "cell_weights")
nc_close(hadcm3nc)

changes <- list()
lower.bounds <- list()
upper.bounds <- list()

for (state in c("PI","LGM"))
{
  for (tsc in names(tsc.in))
  {
    n.all <- length(var.ratio[[state]][[tsc]])
    n <- length(QC.all[[state]][[tsc]]$ind.ok)
    if(n == 0) next
    
    loc.weights <- sapply(1:n.all,
                          function(i)
                          {
                            ind <- interpolate_grid_to_index(
                              QC.all[[state]][[tsc]]$metafilt$Lon[(i-1)%%n+1],
                              QC.all[[state]][[tsc]]$metafilt$Lat[(i-1)%%n+1]
                            )
                            return(cell_weights[ind[1],ind[2]])
                          }
    )
      a <- var.ratio[[paste0(state,"*")]][[tsc]]
      b <- var.ratio[[state]][[tsc]]
      a[which(a<1)] <- 1/(a[which(a<1)])
      b[which(b<1)] <- 1/(b[which(b<1)])
      imp <- length(which(((log10(b) - log10(a)) >0) == T))
      change_log <- weighted_mean(log10(b) - log10(a), loc.weights, na.rm=T)
      change <- 10**change_log

       error.a <- var.ratio.CI[[paste0(state,"*")]][[tsc]]
       error.b <- var.ratio.CI[[state]][[tsc]]
       error.a[,1][which(error.a[,1] < 1)] <- 1/error.a[,1][which(error.a[,1] < 1)] 
       error.a[,2][which(error.a[,2] < 1)] <- 1/error.a[,2][which(error.a[,2] < 1)] 
       error.b[,1][which(error.a[,1] < 1)] <- 1/error.a[,1][which(error.a[,1] < 1)] 
       error.b[,2][which(error.a[,2] < 1)] <- 1/error.a[,2][which(error.a[,2] < 1)] 
       error.a <- apply(error.a, 1, max)
       error.b <- apply(error.b, 1, max)
       
       error.change <- weighted_mean(log10(error.b) + log10(error.a), loc.weights, na.rm=T)
       error.change <- weighted_mean(sqrt((log10(error.b)/(b*log(10)))**2 + (log10(error.a)/(a*log(10)))**2), loc.weights, na.rm=T)
       upper.bound <- 10**(change_log+error.change)
       lower.bound <- 10**(change_log-error.change)
      
    print(paste0("Going from ",state," to ",state,"* at ",tsc,"yrs decreased the",
                 " ratio bias by a factor of ",round(change,3)," /w confidence range(",round(lower.bound,3),",",round(upper.bound,3), ") \\%",
                 " with ",n.all/n," runs with ",n," samples each and overall ",imp, " improved samples."))
    
    changes[[state]][[tsc]] <- change
    upper.bounds[[state]][[tsc]] <- lower.bound
    lower.bounds[[state]][[tsc]] <- upper.bound
    }
}
