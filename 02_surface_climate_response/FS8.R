source("colours.R")
source("functions.R")
load("run_dictionary.RData")
data.dir <- "data/HadCM3/"
load("data/categorized_eruptions.RData")
source("colours.R")

volc_mean <- apply(D,1,mean)

savefile = "02_surface_climate_response/volc_samalas_lgm_pi.RData"

applyOverTSList <- function(ts_list, func, ...)
{
  latest_start <- max(sapply(ts_list, function(x) {x$tstart}))
  
  ts_list <- lapply(ts_list, FUN=
                      function(x)
                      {
                        stopifnot(x$deltat == 30)
                        time_diff <- (latest_start - x$tstart)/30
                        if(time_diff > 0) x$data <- x$data[(time_diff+1):length(x$data)]
                        x$tstart <- latest_start
                        return(x)
                      }
  )
  
  max_length <- min(sapply(ts_list, function(x) {length(x$data)}))
  
  ts_list <- lapply(ts_list, FUN=
                      function(x)
                      {
                        x$data <- x$data[1:max_length]
                        return(x)
                      }
  )
  
  ts.apply <- list(tstart=latest_start, deltat=30)
  ts.apply$data <- apply(simplify2array(lapply(ts_list, function(x) {x$data})), 1, func, ...)
  return(ts.apply)
}

makeAnomaly <- function(ts, start.volc)
{
  t <- ((1:length(ts$data))-1)*30+ts$tstart
  
  start.ts <- which.min(abs(t - ADToDaysSince(tAD$volc[start.volc])))
  
  if(start.ts < 13) return(NA)
  
  ts$data <- ts$data - mean(ts$data[-11:0+start.ts])
  
  return(ts)
}

makeRelative <- function(ts, start.volc)
{
  t <- ((1:length(ts$data))-1)*30+ts$tstart
  
  start.ts <- which.min(abs(t - ADToDaysSince(tAD$volc[start.volc])))
  
  if(start.ts < 13) return(NA)
  
  ts$data <- ts$data/mean(ts$data[-11:0+start.ts])
  
  return(ts)
}

if (!file.exists(savefile))
{
    if(!exists("gmst"))
    {
    gmst <- list()
    for (run_type in c("LGM*","PI*"))
    {
        for (run in names(run_dict)[sapply(run_dict, function(x) {x == run_type})])
        {
        nc <- nc_open(paste(data.dir, "surface_temperature/stationary/",run,".nc",sep=""))
        gmst[[run]] <- nc_var_to_TS(nc, "temp_1")
        gmst[[run]]$data <- fldmean(gmst[[run]]$data, na.rm=T)
        nc_close(nc)
        }
    }
    }

    if(!exists("gmpr"))
    {
    gmpr <- list()
    for (run_type in c("LGM*","PI*"))
    {
        for (run in names(run_dict)[sapply(run_dict, function(x) {x == run_type})])
        {
        nc <- nc_open(paste(data.dir, "precipitation/stationary/",run,".nc",sep=""))
        gmpr[[run]] <- nc_var_to_TS(nc, "precip")
        gmpr[[run]]$data <- gmpr[[run]]$data*31104000.0  #31,104,000 seconds (=1.0 years)
        gmpr[[run]]$data <- fldmean(gmpr[[run]]$data, na.rm=T)
        nc_close(nc)
        }
    }
    }

    if(!exists("sea_ice"))
    {
    sea_ice <- list()
    for (run_type in c("LGM*","PI*"))
    {
        for (run in names(run_dict)[sapply(run_dict, function(x) {x == run_type})])
        {
        nc <- nc_open(paste(data.dir, "sea_ice/stationary/",run,".nc",sep=""))
        sea_ice[[run]] <- nc_var_to_TS(nc, "iceconc")
        if (run_type == "LGM*") sea_ice[[run]]$data <- fldmean(sea_ice[[run]]$data, na.rm=T, mask=!ancil$lsm$`LGM_grid-ocean`)
        if (run_type == "PI*") sea_ice[[run]]$data <- fldmean(sea_ice[[run]]$data, na.rm=T, mask=!ancil$lsm$`PI_grid-ocean`)
        nc_close(nc)
        }
    }
    }

    i.samalas <- which.max(eruptions$peakVal)
    i.samalas.gmpr <- which.min(abs(DaysSinceToAD(((1:length(gmpr$`LGM*`$data))-1)*30+gmpr$`LGM*`$tstart) - tAD$volc[eruptions$start[i.samalas]]))
    i.start <- eruptions$start[i.samalas]
    ind <- max(0,i.start-150):min(length(volc_mean),i.start+350)

    print("temp")
    gmst.mean <- list()
    gmst.mean[["LGM*"]] <- applyOverTSList(lapply(gmst[1:3], 
                                                function(x) {x <- makeAnomaly(x,eruptions$start[i.samalas]); return(x)})
                                        , mean, na.rm=T)
    gmst.mean[["PI*"]] <- applyOverTSList(lapply(gmst[4:6], 
                                                function(x) {x <- makeAnomaly(x,eruptions$start[i.samalas]); return(x)})
                                        , mean, na.rm=T)
    gmst.max <- list()
    gmst.max[["LGM*"]] <- applyOverTSList(lapply(gmst[1:3], 
                                                function(x) {x <- makeAnomaly(x,eruptions$start[i.samalas]); return(x)})
                                        , max, na.rm=T)
    gmst.max[["PI*"]] <- applyOverTSList(lapply(gmst[4:6], 
                                                function(x) {x <- makeAnomaly(x,eruptions$start[i.samalas]); return(x)})
                                        , max, na.rm=T)
    gmst.min <- list()
    gmst.min[["LGM*"]] <- applyOverTSList(lapply(gmst[1:3], 
                                                function(x) {x <- makeAnomaly(x,eruptions$start[i.samalas]); return(x)})
                                        , min, na.rm=T)
    gmst.min[["PI*"]] <- applyOverTSList(lapply(gmst[4:6], 
                                                function(x) {x <- makeAnomaly(x,eruptions$start[i.samalas]); return(x)})
                                        , min, na.rm=T)
    rm(gmst)
    gc()
    print("precip")
    gmpr.mean <- list()
    gmpr.mean[["LGM*"]] <- applyOverTSList(lapply(gmpr[1:3], 
                                                function(x) {x <- makeAnomaly(x,eruptions$start[i.samalas]); return(x)})
                                        , mean, na.rm=T)
    gmpr.mean[["PI*"]] <- applyOverTSList(lapply(gmpr[4:6], 
                                                function(x) {x <- makeAnomaly(x,eruptions$start[i.samalas]); return(x)})
                                        , mean, na.rm=T)
    gmpr.max <- list()
    gmpr.max[["LGM*"]] <- applyOverTSList(lapply(gmpr[1:3], 
                                                function(x) {x <- makeAnomaly(x,eruptions$start[i.samalas]); return(x)})
                                        , max, na.rm=T)
    gmpr.max[["PI*"]] <- applyOverTSList(lapply(gmpr[4:6], 
                                                function(x) {x <- makeAnomaly(x,eruptions$start[i.samalas]); return(x)})
                                        , max, na.rm=T)
    gmpr.min <- list()
    gmpr.min[["LGM*"]] <- applyOverTSList(lapply(gmpr[1:3], 
                                                function(x) {x <- makeAnomaly(x,eruptions$start[i.samalas]); return(x)})
                                        , min, na.rm=T)
    gmpr.min[["PI*"]] <- applyOverTSList(lapply(gmpr[4:6], 
                                                function(x) {x <- makeAnomaly(x,eruptions$start[i.samalas]); return(x)})
                                        , min, na.rm=T)
    rm(gmpr)
    gc()
    print("sea ice")
    sea_ice.mean <- list()
    sea_ice.mean[["LGM*"]] <- applyOverTSList(lapply(sea_ice[1:3], 
                                                    function(x) {x <- makeAnomaly(x,eruptions$start[i.samalas]); return(x)})
                                            , mean, na.rm=T)
    sea_ice.mean[["PI*"]] <- applyOverTSList(lapply(sea_ice[4:6], 
                                                    function(x) {x <- makeAnomaly(x,eruptions$start[i.samalas]); return(x)})
                                            , mean, na.rm=T)
    sea_ice.max <- list()
    sea_ice.max[["LGM*"]] <- applyOverTSList(lapply(sea_ice[1:3], 
                                                    function(x) {x <- makeAnomaly(x,eruptions$start[i.samalas]); return(x)})
                                            , max, na.rm=T)
    sea_ice.max[["PI*"]] <- applyOverTSList(lapply(sea_ice[4:6], 
                                                function(x) {x <- makeAnomaly(x,eruptions$start[i.samalas]); return(x)})
                                            , max, na.rm=T)
    sea_ice.min <- list()
    sea_ice.min[["LGM*"]] <- applyOverTSList(lapply(sea_ice[1:3], 
                                                    function(x) {x <- makeAnomaly(x,eruptions$start[i.samalas]); return(x)})
                                            , min, na.rm=T)
    sea_ice.min[["PI*"]] <- applyOverTSList(lapply(sea_ice[4:6], 
                                                function(x) {x <- makeAnomaly(x,eruptions$start[i.samalas]); return(x)})
                                            , min, na.rm=T)
    rm(sea_ice)
    gc()
    save("gmst.mean","gmst.max","gmst.min","gmpr.mean","gmpr.max","gmpr.min", "sea_ice.mean", "sea_ice.max", "sea_ice.min", "i.samalas", "i.samalas.gmpr", "i.start", "ind", file=savefile)
    } else {
        load(savefile)
    }

COL <- c(COLS[["LGM"]], COLS[["PI"]])
LTYS <- c(1,1)

plot <- T

if(plot)
{
  filenam <- "02_surface_climate_response/047_volc_samalas_global_anomaly.pdf"
  cairo_pdf(file=filenam,width=11,height=3)
}
layout(t(1:3),width=rep(1,3),heights=1)
par(mar=c(2.8,1.8,1.8,0), oma=c(0,0.8,0,2.5), cex=0.8)

zero <- 0.45
zero.volc <- 0.0125
yscale.volc <- 0.625
volcticks <- c(0,0.2,0.4,0.6)
ygrid <- seq(0,1,0.1)

plot(ADToDaysSince(tAD$volc[ind]), volc_mean[ind]*yscale.volc+zero.volc, 
     type="l",xlab="",ylab="", ylim=range(volc_mean),axes=F,col=adjustcolor("black", 0.7))
ylim <- par("usr")[3:4]
yscale <- 0.25
for (i in 2:1)
{
  ind.gmst <- 10:(length(gmst.mean[[i]]$data)-10)
  plot_mean_conf((ind.gmst-1)*30+gmst.mean[[i]]$tstart
               ,filter(gmst.mean[[i]]$data,rep(1/12,12))[ind.gmst]*yscale+zero
               ,filter(gmst.max[[i]]$data,rep(1/12,12))[ind.gmst]*yscale+zero
               ,filter(gmst.min[[i]]$data,rep(1/12,12))[ind.gmst]*yscale+zero
               ,COL[i]
               ,LTYS[i])
}
plot_axis(1, ticks=round(DaysSinceToAD(axTicks(1))), label="Time CE")#, lineVal = 0.75, lineLab = 2
yticks <- seq(-1.5,0.5,0.5)
mtext("GMST Anomaly [K]", font=2, side=3, adj=0)
plot_axis(2, at=yscale*yticks+zero, ticks=yticks)#, lineVal = 0.75, lineLab = 2
abline(h=ygrid*yscale.volc+zero.volc, lty=3, col="darkgray")
abline(v=ADToDaysSince(c(1257.66667,1262.33333)),h=zero, lty=1, col="darkgray")
box()
plot_letter("a")
legend("topright"
       ,lty=LTYS
       ,col=COL
       ,lwd=2
       ,legend = c("LGM*","PI*"))

plot(ADToDaysSince(tAD$volc[ind]), volc_mean[ind]*yscale.volc+zero.volc
     , type="l",xlab="",ylab="", ylim=range(volc_mean),axes=F,col=adjustcolor("black", 0.7))
ylim <- par("usr")[3:4]

yscale <- 0.005
for (i in 2:1)
{
  ind.gmpr <- 10:(length(gmpr.mean[[i]]$data)-10)
  plot_mean_conf((ind.gmpr-1)*30+gmpr.mean[[i]]$tstart
                 ,filter(gmpr.mean[[i]]$data,rep(1/12,12))[ind.gmpr]*yscale+zero
                 ,filter(gmpr.max[[i]]$data,rep(1/12,12))[ind.gmpr]*yscale+zero
                 ,filter(gmpr.min[[i]]$data,rep(1/12,12))[ind.gmpr]*yscale+zero
                 ,COL[i]
                 ,LTYS[i])
}
plot_axis(1, ticks=round(DaysSinceToAD(axTicks(1))), label="Time CE")#, lineVal = 0.75, lineLab = 2
yticks <- seq(-50,25,25)
mtext("GMPR Anomaly [mm/year]", font=2, side=3, adj=0)
plot_axis(2, at=yscale*yticks+zero, ticks=yticks)#, lineVal = 0.75, lineLab = 2
abline(h=ygrid*yscale.volc+zero.volc, lty=3, col="darkgray")
abline(v=ADToDaysSince(c(1257.66667,1262.33333)),h=zero, lty=1, col="darkgray")
box()

plot_letter("b")

plot(ADToDaysSince(tAD$volc[ind]), volc_mean[ind]*yscale.volc+zero.volc,
     type="l",xlab="",ylab="", ylim=range(volc_mean),axes=F,col=adjustcolor("black", 0.7))
ylim <- par("usr")[3:4]

yscale <- 17.5
for (i in 2:1)
{
  ind.sea_ice <- 10:(length(sea_ice.mean[[i]]$data)-10)
  plot_mean_conf((ind.sea_ice-1)*30+sea_ice.mean[[i]]$tstart
                 ,filter(sea_ice.mean[[i]]$data,rep(1/12,12))[ind.sea_ice]*yscale+zero#-1*yscale
                 ,filter(sea_ice.max[[i]]$data,rep(1/12,12))[ind.sea_ice]*yscale+zero#-1*yscale
                 ,filter(sea_ice.min[[i]]$data,rep(1/12,12))[ind.sea_ice]*yscale+zero#-1*yscale
                 ,COL[i]
                 ,LTYS[i])
}
plot_axis(1, ticks=round(DaysSinceToAD(axTicks(1))), label="Time CE")#, lineVal = 0.75, lineLab = 2
yticks <- seq(-5,10,5)*0.1*0.01
mtext("GMICE Anomaly [%]", font=2, side=3, adj=0)
plot_axis(2, at=yscale*yticks+zero, ticks=yticks*100)#, lineVal = 0.75, lineLab = 2
plot_axis(4, at=volcticks*yscale.volc+zero.volc, ticks=volcticks, label="Aerosol Optical Depth",col=adjustcolor("black", 0.7), adj=0)#, lineVal = 0.75, lineLab = 2
abline(h=ygrid*yscale.volc+zero.volc, lty=3, col="darkgray")
abline(v=ADToDaysSince(c(1257.66667,1262.33333)),h=zero, lty=1, col="darkgray")
box()
plot_letter("c")

if(plot) dev.off()