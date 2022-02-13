source("colours.R")
source("functions.R")
load("run_dictionary.RData")
data.dir <- "data/HadCM3/"

#---DATA---#
source("03_spectral_analysis/regional_vs_global.R")

smooth <- 0.02
regional.spectra <- lapply(regional.spectra, function(x) {LogSmooth(x, df.log=smooth, removeLast=0)})
global.spectra <- lapply(global.spectra, function(x) {LogSmooth(x, df.log=smooth, removeLast=0)})

zonal.spectra_smooth <- list()
for(i in names(zonal.spectra)){
  zonal.spectra_smooth[[i]] <- lapply(zonal.spectra[[i]], function(x) {LogSmooth(x, df.log=smooth, removeLast=0)})
}

rm(SP.2d,zonal.spectra)
gc()

#---PLOT---#
saveToPDF <- F
filenam <- "03_spectral_analysis/f3_sat_hadcm3_lats.pdf"

if(saveToPDF) pdf(file=filenam,width=10,height=5)
  
layout(matrix(c(1,1,2,3),2,2),widths=c(1,0.75),heights=c(1.,1))

par(mar=c(2.5,2.7,0.3,0.3), oma=c(0,0,0,0), cex=1)

order <- c(2,1,4,3)
# lgm_uf, lgm_f, pi_uf, pi_f
COL <- c(rep(COLS[["LGM"]],2), rep(adjustcolor(COLS[["PI"]], alpha.f = 0.8),2))
LTYS <- c(2,1,2,1,3,4)
cex.t <- 1
lwd.t <- 2

ylimz <- c(0.003, 2.8)
xlimz <- 1/c(500,2)
  
plot(xlimz, ylimz, log="xy", axes=FALSE, xlab="", ylab="", type="n")
grid()
for (i in 4:1) {
  plot_mean_conf(
    global.spectra[[i]]$freq, global.spectra[[i]]$spec,
    global.spectra[[i]]$lim.2, global.spectra[[i]]$lim.1,
    col=COL[i],
    lty=LTYS[i])
  plot_mean_conf(
    regional.spectra[[i]]$freq, regional.spectra[[i]]$spec,
    regional.spectra[[i]]$lim.2, regional.spectra[[i]]$lim.1,
    col=COL[i],
    lty=LTYS[i])
}
text(1/3, 8/9, "Local")
text(1/4, 0.03, "Global")
  
axis(1,labels=FALSE)
mtext(side=1,at=axTicks(1),1/axTicks(1),line=0.5,cex=cex.t,lwd=lwd.t)
mtext(side=1,"Period [years]",line=1.5,cex=cex.t,lwd=lwd.t)
  
axis(2, at=axTicks(2)[c(1,2,4,5,7,8)], tick=axTicks(2)[c(1,2,4,5,7,8)], labels=FALSE)
mtext(side=2,at=axTicks(2)[seq(2,length(axTicks(2)),by=3)]
        , axTicks(2)[seq(2,length(axTicks(2)),by=3)],line=0.5,cex=cex.t,lwd=lwd.t)
mtext(side=2,expression("PSD [K"^2*" yr]"),line=1.5,cex=cex.t,lwd=lwd.t)
  
legend("bottomleft"
         ,lty=LTYS[order]
         ,col=COL[order]
         ,lwd=lwd.t
         ,legend = c("LGM","LGM*","PI","PI*")[order])
  
plot_letter("a")
box(lwd=lwd.t)

add_unforced = T
type = "lgm"
    
order <- c(4,5, 2,3,1)
    # "tropics", "mid_lats_n", "mid_lats_s", "hight_lats_n", "hight_lats_s"
    COL <- pal
    LTYS <- c(rep(1, 5), rep(2, 5))
    cex.t <- 1
    lwd.t <- 2
    
    ylimz <- c(0.09, 120)
    xlimz <- 1/c(500,2)
    
    plot(xlimz, ylimz, log="xy", axes=FALSE, xlab="", ylab="", type="n")
    grid()
    
    for(i in 1:length(zonal.spectra_smooth[[1]])){ 
      plot_mean_conf(
        zonal.spectra_smooth[[paste0(type,"_f")]][[i]]$freq, zonal.spectra_smooth[[paste0(type,"_f")]][[i]]$spec,
        zonal.spectra_smooth[[paste0(type,"_f")]][[i]]$lim.2, zonal.spectra_smooth[[paste0(type,"_f")]][[i]]$lim.1,
        col=COL[i],
        lty=LTYS[i])
      if(add_unforced){
        plot_mean_conf(
          zonal.spectra_smooth[[paste0(type,"_uf")]][[i]]$freq, zonal.spectra_smooth[[paste0(type,"_uf")]][[i]]$spec,
          zonal.spectra_smooth[[paste0(type,"_uf")]][[i]]$lim.2, zonal.spectra_smooth[[paste0(type,"_uf")]][[i]]$lim.1,
          col=COL[i],
          lty=LTYS[i+5])
      }
    }
    
    axis(1,labels=FALSE)
    mtext(side=1,at=axTicks(1),1/axTicks(1),line=0.5,cex=cex.t,lwd=lwd.t)
    
    axis(2,labels=FALSE)
    mtext(side=2,at=axTicks(2)[seq(1,length(axTicks(2)),by=2)]
          , axTicks(2)[seq(1,length(axTicks(2)),by=2)],line=0.5,cex=cex.t,lwd=lwd.t)

    text(1/310, 95, "LGM")
    plot_letter("b")
    
    legend("topright"
           ,lty=LTYS
           ,col=COL[order]
           ,lwd=lwd.t
           ,legend =  c("high (N)", "high (S)",  "mid (N)", "mid (S)", "trop"),
           ncol=3, cex=cex.t-0.2, 
           inset=c(-0.05,0))
    
    legend("topright"
           ,lty=c(1,2)
           ,col=rep("black",2)
           ,lwd=lwd.t
           ,legend =  c("forced", "unforced"),
           ncol=1, cex=cex.t-0.2,
           inset=c(-0.,0.25))
    
    box(lwd=lwd.t)

type = "pi"
  
  plot(xlimz, ylimz, log="xy", axes=FALSE, xlab="", ylab="", type="n")
  grid()
  
  for(i in 1:length(zonal.spectra_smooth[[1]])){ 
    plot_mean_conf(
      zonal.spectra_smooth[[paste0(type,"_f")]][[i]]$freq, zonal.spectra_smooth[[paste0(type,"_f")]][[i]]$spec,
      zonal.spectra_smooth[[paste0(type,"_f")]][[i]]$lim.2, zonal.spectra_smooth[[paste0(type,"_f")]][[i]]$lim.1,
      col=COL[i],
      lty=LTYS[i])
    if(add_unforced){
      plot_mean_conf(
        zonal.spectra_smooth[[paste0(type,"_uf")]][[i]]$freq, zonal.spectra_smooth[[paste0(type,"_uf")]][[i]]$spec,
        zonal.spectra_smooth[[paste0(type,"_uf")]][[i]]$lim.2, zonal.spectra_smooth[[paste0(type,"_uf")]][[i]]$lim.1,
        col=COL[i],
        lty=LTYS[i+5])
    }
  }
  
  axis(1,labels=FALSE)
  mtext(side=1,at=axTicks(1),1/axTicks(1),line=0.5,cex=cex.t,lwd=lwd.t)
  mtext(side=1,"Period [years]",line=1.5,cex=cex.t,lwd=lwd.t)
  
  axis(2,labels=FALSE)
  mtext(side=2,at=axTicks(2)[seq(1,length(axTicks(2)),by=2)]
        , axTicks(2)[seq(1,length(axTicks(2)),by=2)],line=0.5,cex=cex.t,lwd=lwd.t)
  text(1/380, 95, "PI")
  plot_letter("c")
  

  
  box(lwd=lwd.t)

if(saveToPDF) dev.off()


