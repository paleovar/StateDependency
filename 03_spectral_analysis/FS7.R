source("colours.R")
source("functions.R")
load("run_dictionary.RData")

#---FUNCTIONS---#
apply_transfer <- function(output_spec, input_spec, smooth=0.01){
  s <- transferSpec(list(output_spec, input_spec), input=2, output=1)
  s <- lapply(s$spec, function(x) x[2:which(s$spec$freq >=1)[1]])
  s_smoothed <- LogSmooth(s, df.log=0.01, removeLast = 0)
  s_smoothed$lim.1 <- s$lim.1
  s_smoothed$lim.2 <- s$lim.2
  return(s_smoothed) 
}

#---DATA---#
load("03_spectral_analysis/ebm_regional_vs_global_sat.RData")
regional.spectra1 <- regional.spectra 
global.spectra1 <- global.spectra
rm(SP.2d, global.spectra, regional.spectra)
gc()

load("03_spectral_analysis/ebm_timmean_regional_vs_global_sat.RData")
regional.spectra2 <- regional.spectra
global.spectra2 <- global.spectra 
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

ebm_spectra <- list(global.spectra, regional.spectra, zonal.spectra)
names(ebm_spectra) <- c("global.spectra", "regional.spectra", "zonal.spectra")
rm(global.spectra, regional.spectra, zonal.spectra)
gc()

#LOAD HADCM3 spectra
savefile <- "03_spectral_analysis/regional_vs_global_sat.RData"
load(savefile)
hadcm3_spectra <- list(global.spectra, regional.spectra, SP.2d)
names(hadcm3_spectra) <- c("global.spectra", "regional.spectra", "SP.2d")
rm(global.spectra, regional.spectra, SP.2d, savefile)
gc()

#EBM Spectra: forc + sea ice, forc, unforc + sea ice
#HadCM3 Spectra: forc + sea ice, unforc + sea ice

#HadCM3: forc + sea ice / unforc + sea ice = forc + driven sea ice
#EBM: forc + sea ice / unforc + sea ice = forc + driven sea ice
#EBM: forc + sea ice / forc = sea ice (driven + internal)
#gain: HadCM3 forc + sea ice / EBM forc = HadCM3 sea ice + internal dynamics + nonlinear response
####gain: HadCM3 forc + sea ice / EBM forc + sea ice = internal dynamics + nonlinear response

ratios_f <- list(
  global = list(lgm = apply_transfer(hadcm3_spectra$global.spectra$lgm_f, ebm_spectra$global.spectra$lgm_f), 
                 pi = apply_transfer(hadcm3_spectra$global.spectra$pi_f, ebm_spectra$global.spectra$pi_f)),
  regional = list(lgm = apply_transfer(hadcm3_spectra$regional.spectra$lgm_f, ebm_spectra$regional.spectra$lgm_f),
                  pi = apply_transfer(hadcm3_spectra$regional.spectra$pi_f, ebm_spectra$regional.spectra$pi_f))
)

ratios_f_sea_ice <- list(
  global = list(lgm = apply_transfer(hadcm3_spectra$global.spectra$lgm_f, ebm_spectra$global.spectra$lgm_f_sea_ice), 
                pi = apply_transfer(hadcm3_spectra$global.spectra$pi_f, ebm_spectra$global.spectra$pi_f_sea_ice)),
  regional = list(lgm = apply_transfer(hadcm3_spectra$regional.spectra$lgm_f, ebm_spectra$regional.spectra$lgm_f_sea_ice),
                  pi = apply_transfer(hadcm3_spectra$regional.spectra$pi_f, ebm_spectra$regional.spectra$pi_f_sea_ice))
)

#---PLOT---#
saveToPDF <- F

order <- c(3,4,1,2)
j=2
# global_f_lgm, global_f_pi, global_f_sea_ice_lgm, global_f_sea_ice_pi
COL <- rep(c(COLS[["LGM"]], COLS[["PI"]]),2)
LTYS <- c(3, 3, 1, 1)
cex.t <- 1
lwd.t <- 2

filenam<-"03_spectral_analysis/058_spectral_ratio_regional.pdf"
if(saveToPDF){pdf(file=filenam,width=5,height=4)}

layout(1,1,1)
par(mar=c(3,3,1,1), cex=cex.t, oma=c(0,0,0,0))

ylimz <- c(1.9, 1050)
xlimz <- 1/c(300,5)

plot(xlimz, ylimz, log="xy", axes=FALSE, xlab="", ylab="", type="n", yaxs="i")
grid()
ratios <- c(ratios_f$regional, ratios_f_sea_ice$regional)

for (i in c(2,1,4,3)) { 
  plot_mean_conf(
    ratios[[i]]$freq, ratios[[i]]$spec,
    ratios[[i]]$lim.2, ratios[[i]]$lim.1,
    col=COL[i],
    lty=LTYS[i])
}

text(1/140, 820 , "local: HadCM3 / TransEBM", cex=0.8)

axis(1,labels=FALSE)
mtext(side=1,at=axTicks(1),1/axTicks(1),line=0.5,cex=cex.t,lwd=lwd.t)
mtext(side=1,"Period [years]",line=1.5,cex=cex.t,lwd=lwd.t)

axis(2,labels=FALSE)
mtext(side=2,at=axTicks(2)[seq(1,length(axTicks(2)),by=2)]
      , axTicks(2)[seq(1,length(axTicks(2)),by=2)],line=0.5,cex=cex.t,lwd=lwd.t)
mtext(side=2,expression("Ratio of local mean PSD"),line=1.5,cex=cex.t,lwd=lwd.t)

legend("topleft",
        inset= c(0, 0.08), 
       ,lty=LTYS[order]
       ,col=COL[order]
       ,lwd=lwd.t
       ,legend = c("LGM* / LGM*", "PI* / PI*",  "LGM* / LGM* (SID)", "PI* / PI* (SID)")[order], cex=.8) 

box(lwd=lwd.t)

if(saveToPDF) dev.off()
