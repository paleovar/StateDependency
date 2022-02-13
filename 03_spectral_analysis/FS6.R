source("colours.R")
source("functions.R")
load("run_dictionary.RData")

#---DATA---#
load("data/TransEBM_cell_weights.RData")
load("03_spectral_analysis/ebm_regional_vs_global_sat.RData")
smooth <- 0.02
regional.spectra1 <- lapply(regional.spectra, function(x) {LogSmooth(x, df.log=smooth, removeLast=0)})
global.spectra1 <- lapply(global.spectra, function(x) {LogSmooth(x, df.log=smooth, removeLast=0)})
zonal.spectra1 <- list()
for(i in names(zonal.spectra)){zonal.spectra1[[i]] <- lapply(zonal.spectra[[i]], function(x) {LogSmooth(x, df.log=smooth, removeLast=0)})}
rm(SP.2d, global.spectra, regional.spectra, zonal.spectra)
gc()

load("03_spectral_analysis/ebm_timmean_regional_vs_global_sat.RData")
smooth <- 0.02
regional.spectra2 <- lapply(regional.spectra, function(x) {LogSmooth(x, df.log=smooth, removeLast=0)})
global.spectra2 <- lapply(global.spectra, function(x) {LogSmooth(x, df.log=smooth, removeLast=0)})
zonal.spectra2 <- list()
for(i in names(zonal.spectra)){zonal.spectra2[[i]] <- lapply(zonal.spectra[[i]], function(x) {LogSmooth(x, df.log=smooth, removeLast=0)})}
rm(SP.2d, global.spectra, regional.spectra)
gc()

global.spectra1$lgm_uf <- NULL
global.spectra1$pi_uf <- NULL

global.spectra <- c(global.spectra1, global.spectra2)
names(global.spectra) <- c("lgm_f_sea_ice", "pi_f_sea_ice", names(global.spectra2))

regional.spectra1$lgm_uf <- NULL
regional.spectra1$pi_uf <- NULL

regional.spectra <- c(regional.spectra1, regional.spectra2)
names(regional.spectra) <- c("lgm_f_sea_ice", "pi_f_sea_ice", names(regional.spectra2))

rm(global.spectra1, global.spectra2, regional.spectra1, regional.spectra2)
gc()


#---PLOT---#
saveToPDF <- F

order <- c(1,2, 3, 4)
  # lgm_f, pi_f
  COL <- c(COLS[["LGM"]], COLS[["PI"]], c(adjustcolor(COLS[["LGM"]],0.8), "grey40"))
  
  LTYS <- c(1,1,3,3)
  cex.t <- 1
  lwd.t <- 2
  layout(1,1,1)
  
  filenam<-"03_spectral_analysis/058_ebm_sat_sea_ice_variability.pdf"
  if(saveToPDF){pdf(file=filenam,width=5,height=4)}

  par(mar=c(2.8,2.8,0.5,0.5), cex=cex.t, oma=c(0,0,0,0))

  
  ylimz <- c(0.00008, 1.1)#range(sapply(c(global.spectra,regional.spectra),function(x){x$spec}))#
  xlimz <- 1/c(500,2)
  
  plot(xlimz, ylimz, log="xy", axes=FALSE, xlab="", ylab="", type="n")
  grid()
  for (i in c(4,3, 2, 1)){ 
    plot_mean_conf(
      regional.spectra[[i]]$freq, regional.spectra[[i]]$spec,
      regional.spectra[[i]]$lim.2, regional.spectra[[i]]$lim.1,
      col=COL[i],
      lty=LTYS[i])
    plot_mean_conf(
      global.spectra[[i]]$freq, global.spectra[[i]]$spec,
      global.spectra[[i]]$lim.2, global.spectra[[i]]$lim.1,
      col=COL[i],
      lty=LTYS[i])
    }
  text(1/3, 0.12, "Local")
  text(1/3, 0.004, "Global")
  text(1/2.18, 0.0005, "Local")
  text(1/5, 0.0003, "Global")
  
  axis(1,labels=FALSE)
  mtext(side=1,at=axTicks(1),1/axTicks(1),line=0.5,cex=cex.t,lwd=lwd.t)
  mtext(side=1,"Period [years]",line=1.5,cex=cex.t,lwd=lwd.t)
  
  axis(2,labels=FALSE)
  mtext(side=2,at=axTicks(2)[seq(1,length(axTicks(2)),by=2)]
        , c(format(0.0001, scientific=F), 0.01, 1),line=0.5,cex=cex.t,lwd=lwd.t)
  mtext(side=2,expression("PSD [K"^2*" yr]"),line=1.5,cex=cex.t,lwd=lwd.t)
 
  legend("bottomleft"
           ,lty=LTYS[order]
           ,col=COL[order]
           ,lwd=lwd.t
           ,legend = c("LGM* (SID)","PI* (SID)", "LGM*","PI*")[order])
    
  box(lwd=lwd.t)
  
if(saveToPDF) dev.off()
