data.dir = "data/HadCM3/"


#---FUNCTIONS---#
#if(!file.exists("03_spectral_analysis/ebm_regional_vs_global_sat.RData") | !file.exists("03_spectral_analysis/ebm_regional_vs_global_sat.RData")){
#    source("03_spectral_analysis/TransEBM_regional_vs_global.R")
#}

ebm_nc_to_TS <- function(nc, cut=50*12) 
{
  x <- list()
  x$data <- ncvar_get(nc, "temperature")
  x$time <- ncvar_get(nc, "time")
  x$tstart <- round(x$time[1]/100)
  deltat <- 1/12
  x$time <- seq(x$tstart, by=deltat, length.out=length(x$time))
  d <- dim(x$data)
  if(length(d)==3){
    x$data <- x$data[,,cut:length(x$time)]
  }
  if(length(d)==1){
    x$data <- x$data[cut:length(x$time)]
  }
  x$time <- x$time[cut:length(x$time)]
  return(x)
}

###### DATA #######
for(seaice in c(T,F)){
    if(seaice){
        savefile <- "03_spectral_analysis/ebm_regional_vs_global_sat.RData"
    } else{
        savefile <- "03_spectral_analysis/data/ebm_timmean_regional_vs_global_sat.RData"
    }

    ### Load 2d surface temperatures
    if(!file.exists(savefile)){
        if(seaice){
            EBM_data <- paste0(data.dir, "/EBM1.2/")
            runs <- c("xmzke","xmzkg","xmzkh",
                        "xmzka","xmzkd","xmzki",
                        "xnagd","xnagf","xnagg",
                        "xnagb","xnage","xnagh")
            forced <- c(rep(T,3),rep(F,3),rep(T,3),rep(F,3))#,T,T)
            lgm <- c(rep(T,6),rep(F,6))
        } else {
            EBM_data <- paste0(data.dir, "/EBM1.2/timmean/")
            runs <- c("xmzke","xmzkg","xmzkh",
                        "xnagd","xnagf","xnagg")
            forced <- c(rep(T,6))
            lgm <- c(rep(T,3),rep(F,3))
        }

        ebmnc <- nc_open(paste0(data.dir, "/EBM1.2/xmzkg/xmzkg_ebm_full_cdo.nc"))
        cell_weights <- ncvar_get(ebmnc, "temperature")
        longitude <- ncvar_get(ebmnc, "longitude")
        latitude <- ncvar_get(ebmnc, "latitude")
        d <- dim(ncvar_get(ebmnc, "temperature"))
        nc_close(ebmnc)

        zones_ind <- get_zones(d, latitude)
        load("data/TransEBM_cell_weights.RData")

    global.spectra <- list()
    regional.spectra <- list()
    zonal.spectra <- list()
    SP.2d <- list()
    SP.2d_lats <- list()
    SP.global <- list()
    for (i.run in 1:length(runs))
    {                     
        print(i.run)
        run <- runs[i.run]
        # ../stationary is detrended + deseasonalized by subtracting average seasonality
        nc <- nc_open(paste(paste0(EBM_data, "/EBM1.2/stationary/",run,".nc"),sep=""))
        temp.2d <- ebm_nc_to_TS(nc) 
        nc_close(nc)
        
        nc <- nc_open(paste(paste0(EBM_data, "/EBM1.2/stationary/",run,"_fldmean.nc"),sep=""))
        temp.2d_fldmean <- ebm_nc_to_TS(nc) 
        nc_close(nc)
        
        #Compute spectra
        y <- ts(
            na.approx(
            temp.2d_fldmean$data), 
            start=temp.2d$tstart, deltat=1/12) 

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
                y <- ts(na.approx(temp.2d$data[lon,lat,]),
                        start=temp.2d$tstart, deltat=1/12)
                spectra[[spectra.ind]] <- SpecMTM(na.contiguous(y))
            }
            }
        }
        
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
    if(seaice){global.spectra[["lgm_uf"]] <- MeanSpec(specList = SP.global[ lgm & !forced ])$spec}
    global.spectra[["lgm_f"]]  <- MeanSpec(specList = SP.global[ lgm &  forced ])$spec
    if(seaice){global.spectra[["pi_uf"]]  <- MeanSpec(specList = SP.global[!lgm & !forced ])$spec}
    global.spectra[["pi_f"]]   <- MeanSpec(specList = SP.global[!lgm &  forced ])$spec

    # Assemble regional spectra
    SP.2d <- lapply(SP.2d, function(x) {x$spec$dof <- x$spec$dof/x$nRecord; return(x$spec)}) 
    if(seaice){regional.spectra[["lgm_uf"]] <- MeanSpec(specList = SP.2d[ lgm & !forced ])$spec}
    regional.spectra[["lgm_f"]]  <- MeanSpec(specList = SP.2d[ lgm &  forced ])$spec
    if(seaice){regional.spectra[["pi_uf"]]  <- MeanSpec(specList = SP.2d[!lgm & !forced ])$spec}
    regional.spectra[["pi_f"]]   <- MeanSpec(specList = SP.2d[!lgm &  forced ])$spec
    
    # Assemble zonal spectra
    SP.2d_lats <- lapply(SP.2d_lats, function(x){lapply(x, function(x){x$spec$dof <- x$spec$dof/x$nRecord; return(x$spec)})}) 
    for(i in names(SP.2d_lats[[1]])){
        if(seaice){zonal.spectra[["lgm_uf"]][[i]] <- MeanSpec(lapply(SP.2d_lats[lgm & !forced], function(x) x[[i]]))$spec}
        zonal.spectra[["lgm_f"]][[i]] <- MeanSpec(lapply(SP.2d_lats[lgm&forced], function(x) x[[i]]))$spec
        if(seaice){zonal.spectra[["pi_uf"]][[i]] <- MeanSpec(lapply(SP.2d_lats[!lgm & !forced], function(x) x[[i]]))$spec}
        zonal.spectra[["pi_f"]][[i]] <- MeanSpec(lapply(SP.2d_lats[!lgm & forced], function(x) x[[i]]))$spec
    }
    
    save("global.spectra","regional.spectra", "SP.2d", "zonal.spectra", file=savefile)
    }
}
