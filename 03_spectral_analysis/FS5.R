source("colours.R")
source("functions.R")
load("run_dictionary.RData")
data.dir <- "data/HadCM3/"

savefile <- "03_spectral_analysis/sea_ice_glob.RData"

#---DATA---#
if(!file.exists(savefile))
{
  sea_ice <- list()
  for (run_type in c("LGM*", "LGM", "PI*", "PI"))
  {
    for (run in names(run_dict)[sapply(run_dict, function(x) {x == run_type})])
    {
      print(run)
      nc <- nc_open(paste(data.dir, "/sea_ice/",run,".nc",sep=""))
      tmp <- nc_var_to_TS(nc, "iceconc")
      if (run_type == "LGM*" | run_type == "LGM") sea_ice[[run]]$data <- fldmean(tmp$data, na.rm=T, mask=!ancil$lsm$`LGM_grid-ocean`)
      if (run_type == "PI*" | run_type == "PI") sea_ice[[run]]$data <- fldmean(tmp$data, na.rm=T, mask=!ancil$lsm$`PI_grid-ocean`)
      rm(tmp)
      nc_close(nc)
      gc()
    }
  }
  sea_ice <- lapply(sea_ice, function(x) MakeEquidistant(
   # t.x=seq(1, length(x$data)+12)/12 - 1/12, by=1/12), t.y=x$data, time.target=seq(1, (length(x$data)+12)/12 - 1/12, by=1/12)) #for monthly resolution
    t.x=seq(1, length(x$data), by=1), t.y=x$data, time.target=seq(1, length(x$data), by=1)) #for monthly resolution
  )
  
  for(i in c("xmzkg", "xmzkh", "xmzki", "xnagh", "xnagd", "xnage")){
    sea_ice[[i]] <- tseries::na.remove(sea_ice[[i]]) #cut the end if there are NA values
  }
  save(sea_ice, file=savefile)
} else {load(savefile)}

lgm <- c(rep(T,6), rep(F,6))
forced <- c(rep(c(rep(T,3), rep(F,3)),2))

spec <- lapply(sea_ice, function(x){SpecMTM(x, detrend=T)})

spec.mean <- list()
spec.mean[["lgm_uf"]] <- MeanSpec(specList = spec[lgm & !forced])$spec
spec.mean[["lgm_f"]] <- MeanSpec(specList = spec[lgm & forced])$spec
spec.mean[["pi_uf"]] <- MeanSpec(specList = spec[!lgm & !forced])$spec
spec.mean[["pi_f"]] <- MeanSpec(specList = spec[!lgm & forced])$spec

smooth <- 0.04
spec.mean <- lapply(spec.mean, LogSmooth, df.log=smooth, removeLast=10)

#---PLOT---#
layout(1,1)
par(mar=c(3,3,0,0), cex=cex.t, oma=c(0,0,0,0))

# lgm_uf, lgm_f, pi_uf, pi_f
COLS <- c(rep(COLS$LGM, 2), rep(COLS$PI, 2))
order  <- c(2,1,4,3)
LTYS <- c(2,1,2,1)
cex.t <- 1
lwd.t <- 2

saveToPDF <- F

if(saveToPDF)
{
  filenam<-"03_spectral_analysis/064_sea_ice_variability.pdf"
  pdf(file=filenam,width=5,height=4)
}

par(mar=c(3,3,2,1), cex=cex.t, oma=c(0,0,0,0))

xlimz <- 1/c(500,2.3)
ylimz <- c(0.001,0.1)#range(sapply(spec.mean, function(x) {range(c(x$lim.1,x$lim.2))}))

plot(xlimz, ylimz, log="xy", axes=FALSE, xlab="", ylab="", type="n")
grid()
for (i in 1:4) {
  #SPECTRA RESCALED BY FACTOR OF 1000 TO GET UNIT %^2 yr 
  plot_mean_conf(
    spec.mean[[i]]$freq, spec.mean[[i]]$spec*1000,
    spec.mean[[i]]$lim.2*1000, spec.mean[[i]]$lim.1*1000,
    col=COLS[i],
    lty=LTYS[i])
}

axis(1,labels=FALSE)
mtext(side=1,at=axTicks(1),1/(axTicks(1)),line=0.5,cex=cex.t,lwd=lwd.t)
mtext(side=1,"Period [years]",line=1.5,cex=cex.t,lwd=lwd.t)

axis(2,labels=FALSE)
mtext(side=2,at=axTicks(2)[seq(2,length(axTicks(2)),by=2)]
      , axTicks(2)[seq(2,length(axTicks(2)),by=2)],line=0.5,cex=cex.t,lwd=lwd.t)
mtext(side=2,expression("PSD [%"^2*" yr]"),line=1.5,cex=cex.t,lwd=lwd.t)

legend("topright"
       ,lty=LTYS[order]
       ,col=COLS[order]
       ,legend = c("LGM","LGM*","PI","PI*")[order]
       ,lwd=lwd.t
       ,cex=1)

box(cex=cex.t)

if(saveToPDF) dev.off()
