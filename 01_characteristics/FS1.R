source("colours.R")
source("functions.R")
load("run_dictionary.RData")
load("ancil.RData")
data.dir <- "data/HadCM3/" #Generate fldmean prior to this analysis by using the `getFldmean.sh` script in `/HadCM3/`

#---FUNCTIONS---#
des <- c(0.5,rep(1,11),0.5)/12

#---DATA---#
filename <- "01_characteristics/057_state_scatter"

runs <- names(run_dict)
LGM <- unlist(run_dict, use.names = FALSE) == "LGM" | unlist(run_dict, use.names = FALSE) == "LGM*"
forced <- unlist(run_dict, use.names = FALSE) == "LGM*" | unlist(run_dict, use.names = FALSE) == "PI*"

solar <- read.csv("data/hadcm3_solar_forcing.csv")

load("data/crowley2012.RData")

if(!file.exists(paste(filename,".RData",sep="")))
{
  gmst <- list()
  for (run in runs) {
    nc <- nc_open(paste(data.dir, "surface_temperature/",run,".nc",sep="")) #alternatively compute fldmean using `./data/getFldmean.sh`
    gmst[[run]] <- nc_var_to_TS(nc, "temp_1")
    gmst[[run]]$data <- fldmean(gmst[[run]]$data, na.rm=T)
    nc_close(nc)
  }
  
  timesteps <- sort(as.integer(unique(unlist(sapply(gmst, function(x) x$time)))))
  n_timesteps <- length(timesteps)
  timesteps.range <- range(timesteps)
  gmst_sync <- matrix(NA,nrow=n_timesteps,ncol=length(runs))
  
  for (i in 1:length(runs)) {
    ind <- which(timesteps %in% gmst[[i]]$time)
    # We can do this since timesteps and thus the first dimension has _all_ timesteps
    gmst_sync[ind,i] <- gmst[[i]]$data
  }
  
  save("data","timesteps","gmst_sync",file=paste(filename,".RData",sep=""))
} else {
  load(paste(filename,".RData",sep=""))
}

#---PLOT---#
saveToPDF <- F

COL <- c(COLS[["LGM"]], COLS[["PI"]], "grey")
xlimz <- c(6,20)
ylimz <- c(6,20)

if (saveToPDF) cairo_pdf(filename=paste0(filename, ".pdf"), width=9, height=5)

# Set up plot
layout(matrix(c(1,1,2,3),2,2),widths=c(1,0.4),heights=1)
lwd.t <- 2
par(mar=c(3,3,1,5), oma=c(0,0,0,0), cex=1)

xlim <- range(DaysSinceToAD(timesteps))
plot(0,xlim=xlim,ylim=c(0,1), 
     type="n", xlab="", ylab="", yaxt='n', xaxt='n', bty="n",yaxs="i",
     panel.first={abline(v=seq(800,2000,100),col="lightgrey",lty=3)})

# PI GMST
y <- filter(apply(gmst_sync[,!LGM&forced],1,mean,na.rm=T),des)
ymean <- mean(y, na.rm=T)
scale <- 0.15
zero <- 0.92
lines(DaysSinceToAD(timesteps),(y-ymean)*scale+zero,col=COL[2])

ticks <- c(14.5,15.5)
plot_axis(2, at=(ticks-ymean)*scale+zero, ticks=ticks,col=COL[2],
          "PI* GMST [°C]", adj=1)

# LGM GMST
y <- filter(apply(gmst_sync[,LGM&forced],1,mean,na.rm=T),des)
ymean <- mean(y, na.rm=T)
zero <- 0.66
lines(DaysSinceToAD(timesteps),(y-ymean)*scale+zero,col=COL[1])

ticks <- seq(9,10,1)
plot_axis(4, at=(ticks-ymean)*scale+zero, ticks=ticks,col=COL[1],
          "LGM* GMST [°C]", adj=0.75)

# TSI
col.tsi <- pal[[1]]
solar.x <- DaysSinceToAD(solar$daysSince_1.1.year_hadcm3)
y <- solar$TSI[solar.x<xlim[2]]
ymean <- mean(y, na.rm=T)
scale <- 0.1
zero <- 0.35
lines(solar.x[solar.x<xlim[2]],(y-ymean)*scale+zero,col=col.tsi)

ticks <- seq(1365,1366,1)
plot_axis(2, at=(ticks-ymean)*scale+zero, ticks=ticks, col=col.tsi, 
          label=expression("Total Solar Irradiance [W/m"^2*"]"), adj=0.2)

# AOD
col.aod <- pal[[4]]
y <- apply(forcing, 1, mean)
ymean <- mean(y, na.rm=T)
scale <- 0.4
zero <- 0.02
lines(as.double(names(y)),(y-ymean)*scale+zero, col=col.aod)

ticks <- seq(0,0.6,0.2)
plot_axis(4, at=(ticks-ymean)*scale+zero, ticks=ticks, col=col.aod,
          label="Aerosol Optical Depth", adj=0)

box()

plot_axis(1, label="Time CE")
plot_letter("a")
par(mar=c(3,3,1,1))
plot(detrend(filter(apply(gmst_sync[,LGM&forced],1,mean,na.rm=T),des)), 
     detrend(filter(apply(gmst_sync[,LGM&!forced],1,mean,na.rm=T),des)),
     xlim=c(-1.3,0.3), ylim=c(-0.8,0.8), 
     xlab="", ylab="", yaxt='n', xaxt='n', type="l", bty="n", pch=19, cex=0.1,
     panel.first={grid();abline(h=0,v=0,col="lightgrey",lty=1,lwd=2)})
box()

plot_axis(1, label="LGM* GMST anom. [K]")
yticks <- c(-0.5,0,0.5)
plot_axis(2, at=yticks, ticks=yticks, label="LGM GMST anom. [K]")
plot_letter("b")

plot(detrend(filter(apply(gmst_sync[,LGM&forced],1,mean,na.rm=T),des)), 
     detrend(filter(apply(gmst_sync[,!LGM&forced],1,mean,na.rm=T),des)), 
     xlim=c(-1.3,0.3), ylim=c(-1.3,0.3),
     xlab="", ylab="", yaxt='n', xaxt='n', type="l", bty="n", pch=19, cex=0.1,
     panel.first= {
       grid();
       abline(h=0,v=0,col="lightgrey",lty=1,lwd=2);
       lines(c(-2,2),c(-2,2), lwd=1, col="red", lty=1);
     })
box()

plot_axis(1, label="LGM* GMST anom. [K]")
plot_axis(2, label="PI* GMST anom. [K]")
plot_letter("c")

if (saveToPDF) dev.off()
