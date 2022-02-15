##THIS SCRIPT WILL BE CALLED FROM F2.R and FS3.R

#---FUNCTIONS---#
run_types <- c("LGM*","PI*")

getNearestIndex.vec <- Vectorize(getNearestIndex, vectorize.args = "x")
iTo2Dx <- function(i,dim2) return((i-1)%/%dim2+1)
iTo2Dy <- function(i,dim2) return((i-1)%%dim2+1)

getEpochData <- function(all_data, var, path, ncvarname, mult=1, boot=F)
{
  if(!(var %in% names(all_data)))
  {
    all_data[[var]] <- list()
    all_data[[var]][["mean.active"]] <- list()
    all_data[[var]][["mean"]] <- list()
    if(boot)
    {
      all_data[[var]][["lower-0.99"]] <- list()
      all_data[[var]][["upper-0.99"]] <- list()
    }
    all_data[[var]][["sd"]] <- list()
    for (run_type in run_types)
    {
      print(run_type)
      for (run in names(run_dict)[sapply(run_dict, function(x) {x == run_type})])
      {
        print(run)
        nc <- nc_open(paste(path,run,".nc",sep=""))
        raw <- nc_var_to_TS(nc, ncvarname)
        raw$data <- raw$data*mult
        nc_close(nc)
        
        print("filter")
        raw$data <- aperm(apply(raw$data, c(1,2), filter, filter=rep(1/12,12)),c(2,3,1)) #moving average 12 months
        
        d <- dim(raw$data)
        
        ind.volc.active <- unique(getNearestIndex.vec(raw$time, volc.active.DaysSince))
        ind.volc.active <- ind.volc.active[!(ind.volc.active == length(raw$time)) & !(ind.volc.active == 1)]
        
        print("mean.active")
        all_data[[var]][["mean.active"]][[run]] <- apply(raw$data[,,ind.volc.active], c(1,2), mean, na.rm=T)
        print("mean")
        all_data[[var]][["mean"]][[run]] <- apply(raw$data[,,], c(1,2), mean, na.rm=T)
        print("sd")
        all_data[[var]][["sd"]][[run]] <- apply(raw$data[,,], c(1,2), sd, na.rm=T)
        
        if(boot)
        {
          print("booting")
          perc <- mclapply(1:(d[1]*d[2]), function(i.grid) {
            b <- tsboot(raw$data[iTo2Dx(i.grid,d[2]),iTo2Dy(i.grid,d[2]),],
                        function(u,i.boot) mean(u[i.boot], na.rm=T), R = 400, sim="fixed", l=48)
            ci <- try(boot.ci(b, type = c("perc"), conf=0.99), silent=T)
            
            if(inherits(ci, "try-error"))
            {
              warning("Too many NA in time-series, returning NA")
              return(list(percent=rep(NA,5)))
            } else  {
              return(ci)
            }
          })
          all_data[[var]][["lower-0.99"]][[run]] <- t(array(do.call(c,
                                                             lapply(perc, function(x) x$percent[4])), d[2:1]))
          all_data[[var]][["upper-0.99"]][[run]] <- t(array(do.call(c,
                                                             lapply(perc, function(x) x$percent[5])), d[2:1]))
        }
      }
    }
    rm(raw)
    gc()
  }
  
  return(all_data)
}

#---DATA---#
if(!exists("eruptions"))
{
  load("data/categorized_eruptions.RData")
  volc_mean <- apply(D,1,mean)
}

if(!exists("booting")){
  booting <- F
  print("no booting chosen")
}
if(booting){booted = "boot"}else{booted = "wo_boot"}

volc.active.DaysSince <- ADToDaysSince(as.double(names(which(volc_mean > 0.13))))

if(!exists("all_data")) 
{
  if(!file.exists(paste0("02_surface_climate_response/048_stdanom_all_data_", booted,".RData"))){
    all_data <- list()
    print("temp")
    all_data <- getEpochData(all_data, "Surface Temperature", paste0(data.dir, "surface_temperature/surfa"), "temp_1", boot=booting) #T
    print("Precip")
    all_data <- getEpochData(all_data, "Precipitation",  paste0(data.dir, "precipitation/preci"), "precip",mult=31104000.0, boot=booting) #T  #31,104,000 seconds (=1.0 years)
    print("sea level")
    all_data <- getEpochData(all_data, "Sea Level Pressure", paste0(data.dir, "sea_level_pressure/sea_l"), "p", boot=booting)#T
    print("wind v")
    all_data <- getEpochData(all_data, "Southerly Winds", paste0(data.dir, "wind/wind_"), "v", boot=booting) #T
    print("wind u")
    all_data <- getEpochData(all_data, "Westerly Winds", paste0(data.dir, "wind/wind_"), "u", boot=booting) #T
    
    save("all_data", file=paste0("02_surface_climate/response/048_stdanom_all_data_", booted, ".RData"))
  } else load(paste0("02_surface_climate_response/048_stdanom_all_data_", booted, ".RData"))
} 

for (var in c("Surface Temperature", "Precipitation", "Sea Level Pressure", 
              "Southerly Winds", "Westerly Winds"))
{
  all_data[[var]][["std.anom"]] <- list()
  for (run_type in run_types)
  {
    all_data[[var]][["std.anom"]][[run_type]] <- apply(simplify2array(lapply(
      names(run_dict)[sapply(run_dict, function(x) {x == run_type})], 
      function(run) 
      {(all_data[[var]][["mean.active"]][[run]]-all_data[[var]][["mean"]][[run]])/all_data[[var]][["sd"]][[run]]}))
      , c(1,2), mean)
  }
  
  for (run_type in run_types)
  {
    all_data[[var]][["mean"]][[run_type]] <- apply(simplify2array(lapply(
      names(run_dict)[sapply(run_dict, function(x) {x == run_type})], 
      function(run)
      {all_data[[var]][["mean"]][[run]]}))
      , c(1,2), mean)
  }
  
  for (run_type in run_types)
  {
    all_data[[var]][["sd"]][[run_type]] <- sqrt(apply(simplify2array(lapply(
      names(run_dict)[sapply(run_dict, function(x) {x == run_type})], 
      function(run)
      {all_data[[var]][["sd"]][[run]]}))
      , c(1,2), sum))
  }
  
  if(T & !grepl("Winds",var)) #boot
  {
    for (conf in c(0.99))
    {
      for (bounds in c("lower","upper"))
      {
        all_data[[var]][[paste("anom-",conf,sep="")]][[bounds]] <- list()
        for (run_type in run_types)
        {
          all_data[[var]][[paste("anom-",conf,sep="")]][[bounds]] <- apply(simplify2array(lapply(
            names(run_dict)[sapply(run_dict, function(x) {x == run_type})],
            function(run)
            {(all_data[[var]][[paste(bounds,conf,sep="-")]][[run]]-all_data[[var]][["mean"]][[run]])/all_data[[var]][["sd"]][[run]]}))
            , c(1,2), mean)
        }
      }
    }
  }
}
