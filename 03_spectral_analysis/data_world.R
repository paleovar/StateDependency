#THIS SCRIPT WILL BE CALLED FROM FS4.R
load("lsm_regions.RData")
ind.world <- which(lsm$world, arr.ind=T)
run_types <- c("LGM","LGM*","PI","PI*")

#exclude xmzki and xnagh on 200-500 yr timescale
run_dict$xmzki <- NULL
run_dict$xnagh <- NULL

path <- paste0(data.dir, "/surface_temperature/stationary/")
ncvar_name <- ncvarnames
  
if(tscale1 != F & tscale2 !=F){
    savefile <- paste("03_spectral_analysis/spectral_analysis_world_",ncvar_name,"_", tscale1, "_", tscale2, ".RData",sep="")
} else {
    print("define tscale1 and tscale2 to estimate scaling coefficients")
}

if(!file.exists(savefile))
  {
    scaling <- list()
    scaling_sdev <- list()
    
    # Get scaling, sd of scaling and mean spectrum for every gridbox in antarctica
    for (run_type in run_types)
    {
      gc()
      dat <- list()
      scaling[[run_type]] <- array(NA, dim(lsm[[1]]))
      scaling_sdev[[run_type]] <- array(NA, dim(lsm[[1]]))
      print(run_type)
      for (i.ind in 1:dim(ind.world)[1])
      {
        spec_tmp <- list()
        ind <- ind.world[i.ind,]
        for (run in names(run_dict)[sapply(run_dict, function(x) {x == run_type})])
        {
          if(!(run %in% names(dat))) {
            dat[[run]] <- nc_var_to_TS(nc_open(paste(path,run,".nc",sep="")), ncvar_name)
            if(ncvar_name == "precip") dat[[run]]$data <- dat[[run]]$data*31104000.0
          }
          spec_tmp[[run]] <- SpecMTM(to.ts(dat[[run]], ind))
          if(sum(!is.na(spec_tmp[[run]]$spec)) == 0) spec_tmp[[run]] <- NULL
        }
        if(length(spec_tmp) == 0) next
        spec <- MeanSpec(spec_tmp)$spec
        if(tscale1 !=F & tscale2 != F){
          fit <- SlopeFit(spec, 1/tscale1, 1/tscale2, bDebug=F)
        }
        scaling[[run_type]][ind[1],ind[2]] <- fit$slope
        scaling_sdev[[run_type]][ind[1],ind[2]] <- fit$slopesd 
    }
    }
    save("scaling","scaling_sdev", file=savefile)
  } else {
    load(savefile)
  }
