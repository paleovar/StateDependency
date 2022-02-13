#THIS SCRIPT WILL BE CALLED FROM F3.R

#---DATA---#
hadcm3nc <- nc_open("hadcm_weights.nc")
cell_weights <- ncvar_get(hadcm3nc, "cell_weights")
longitude <- ncvar_get(hadcm3nc, "longitude")
latitude <- ncvar_get(hadcm3nc, "latitude")
d <-  c(length(longitude),length(latitude)) 
nc_close(hadcm3nc)

zones_ind <- get_zones(d, latitude)

runs <- c("xmzke","xmzkg","xmzkh",
          "xmzka","xmzkd","xmzki",
          "xnagd","xnagf","xnagg",
          "xnagb","xnage","xnagh")
forced <- c(rep(T,3),rep(F,3),rep(T,3),rep(F,3))#,T,T)
lgm <- c(rep(T,6),rep(F,6))

savefile <- "03_spectral_analysis/regional_vs_global_sat.RData"

### Load 2d surface temperatures
if (!file.exists(savefile))
{
  global.spectra <- list()
  regional.spectra <- list()
  zonal.spectra <- list()
  # Get data
  SP.2d <- list()
  SP.global <- list()
  SP.2d_lats <- list()
  for (i.run in 1:length(runs))
  {
    print(i.run)
    run <- runs[i.run]
    #data is detrended + deseasonalized by subtracting average seasonality
    nc <- nc_open(paste(paste0(data.dir, "/surface_temperature/stationary/",run,".nc"),sep=""))
    temp.2d <- nc_var_to_TS(nc, "temp_1")
    nc_close(nc)
    
    d <- dim(temp.2d$data)
    
    #Compute spectra
      y <- ts(
        na.approx(
          fldmean(
            temp.2d$data, na.rm=T
          ),
          na.rm=F, maxgap=1), 
        start=temp.2d$tstart/360, deltat=1/12)

    SP.global[[run]] <- SpecMTM(na.contiguous(y))
    
    ## regional spectra
    # Loop over gridpoints and generate flat list of spectra for every run
    spectra <- list()
    spectra.mask <- array(T, d[1]*d[2])
      for (lat in 1:d[2])
      {
        for (lon in 1:d[1])
        {
          spectra.ind <- (lat-1)*d[1]+lon
          # If too many NA or land
          if(sum(is.na(temp.2d$data[lon,lat,])) > 0.1*d[3])
          {
            spectra.mask[spectra.ind] <- F
            spectra[[spectra.ind]] <- NA
          } else {
            y <- ts(na.approx(temp.2d$data[lon,lat,], na.rm=F, maxgap=1),
                    start=temp.2d$tstart/360, deltat=1/12)
            spectra[[spectra.ind]] <- SpecMTM(na.contiguous(y))
          }
        }
      }
    SP.2d[[run]] <- MeanSpec(specList = spectra[spectra.mask], weights=as.vector(cell_weights)[spectra.mask])
    
    zones_mask <- lapply(zones_ind, function(x){
      tmp.mask <- spectra.mask
      tmp.mask[which(!1:dim(tmp.mask) %in% x)] <- F
      return(tmp.mask)
    })
    
    SP.2d_lats[[run]] <- lapply(zones_mask, function(x){
      MeanSpec(specList = spectra[x], weights=as.vector(cell_weights)[x])
    }
    )

    rm(temp.2d, nc, spectra)
    gc()
  }
  
  # Assemble global spectra
  global.spectra[["lgm_uf"]] <- MeanSpec(specList = SP.global[ lgm & !forced ])$spec
  global.spectra[["lgm_f"]]  <- MeanSpec(specList = SP.global[ lgm &  forced ])$spec
  global.spectra[["pi_uf"]]  <- MeanSpec(specList = SP.global[!lgm & !forced ])$spec
  global.spectra[["pi_f"]]   <- MeanSpec(specList = SP.global[!lgm &  forced ])$spec
  
  # Assemble regional spectra
  SP.2d <- lapply(SP.2d, function(x) {x$spec$dof <- x$spec$dof/x$nRecord; return(x$spec)})
  regional.spectra[["lgm_uf"]] <- MeanSpec(specList = SP.2d[ lgm & !forced ])$spec
  regional.spectra[["lgm_f"]]  <- MeanSpec(specList = SP.2d[ lgm &  forced ])$spec
  regional.spectra[["pi_uf"]]  <- MeanSpec(specList = SP.2d[!lgm & !forced ])$spec
  regional.spectra[["pi_f"]]   <- MeanSpec(specList = SP.2d[!lgm &  forced ])$spec
  
  SP.2d_lats <- lapply(SP.2d_lats, function(x){lapply(x, function(x){x$spec$dof <- x$spec$dof/x$nRecord; return(x$spec)})}) 
  for(i in names(SP.2d_lats[[1]])){
    zonal.spectra[["lgm_uf"]][[i]] <- MeanSpec(lapply(SP.2d_lats[lgm & !forced], function(x) x[[i]]))$spec
    zonal.spectra[["lgm_f"]][[i]] <- MeanSpec(lapply(SP.2d_lats[lgm&forced], function(x) x[[i]]))$spec
    zonal.spectra[["pi_uf"]][[i]] <- MeanSpec(lapply(SP.2d_lats[!lgm & !forced], function(x) x[[i]]))$spec
    zonal.spectra[["pi_f"]][[i]] <- MeanSpec(lapply(SP.2d_lats[!lgm & forced], function(x) x[[i]]))$spec
  }

  save("global.spectra","regional.spectra", "SP.2d", "zonal.spectra", file=savefile)
} else {
  load(savefile)
}