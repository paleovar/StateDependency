require(openxlsx)
source("colours.R")
source("functions.R")
load("run_dictionary.RData")
data.dir <- "data/HadCM3/"

#---DATA---#
# Get amoc
if(!file.exists("data/amoc_hadcm3.xlsx")){
    source("03_spectral_analysis/amoc.R")
    } else {amoc_raw <- read.xlsx("data/amoc_hadcm3.xlsx")}

runs <- names(amoc_raw[grepl("x", colnames(amoc_raw))])
states <- list("LGM"     =c(F,F,F,T,T,T,rep(F,8)),
               "LGM*"    =c(T,T,T,rep(F,11)),
               "PI"      =c(rep(F,9),T,T,T,F,F),
               "PI*"     =c(rep(F,6),T,T,T,rep(F,5)))

order <- c(2,1,4,3)

filename <- "03_spectral_analysis/053_amoc_spectra_decorr.pdf"
amoc.detrended <- lapply(amoc_raw[runs], detrend)

amoc.decorr <- lapply(amoc.detrended, get_corr_length)
corr_length <- array(NA,length(states))
corr_length.sd <- array(NA,length(states))

spec <- lapply(lapply(amoc.detrended, na.contiguous), SpecMTM) #warnings are fine, timestep is 1

spec.mean <- list()

for (i in 1:length(states))
{
  corr_length[i] <- mean(unlist(amoc.decorr[states[[i]]]))
  corr_length.sd[i] <- sd(unlist(amoc.decorr[states[[i]]]))/sqrt(sum(states[[i]])*100)
  spec.mean[[names(states)[i]]] <- MeanSpec(specList = spec[states[[i]]])$spec
}

smooth = 0.04
spec.mean <- lapply(spec.mean, LogSmooth, df.log=smooth, removeLast=5)

#---PLOT---#
plot <- F
if(plot) cairo_pdf(file=filename,width=8,height=4)

# LGM LGM* PI PI*
COL <- c(rep(COLS$LGM, 2), rep(COLS$PI, 2))
LTYS <- c(2,1,2,1)
PCHS <- c(19,17,19,17)

layout(matrix(1:2,1,2), widths=c(1.2,0.8))
par(mar=c(3,3,1,1), oma=c(0,0,0,0))
xlimz <- 1/c(200,2)
ylimz <- range(sapply(spec.mean, function(x) {range(c(x$lim.1,x$lim.2))}))*c(1,1)
plot(xlimz, ylimz, log="xy", axes=FALSE, xlab="", ylab="", type="n", xaxs="i")
for (i in 1:length(states))
{
  plot_mean_conf(spec.mean[[i]]$freq,spec.mean[[i]]$spec,spec.mean[[i]]$lim.1,spec.mean[[i]]$lim.2,
                 col=COL[i],lty=LTYS[i])
}
plot_axis(1, ticks=1/axTicks(1), label="Period [years]")
plot_axis(2, label=expression("PSD [Sv"^2*" yr]"))
grid()
legend("topright"
       ,lty=LTYS[order]
       ,col=COL[order]
       ,legend = c("LGM","LGM*","PI","PI*")[order]
       ,lwd=2
       ,cex=1)
plot_letter("a")
box()

par(mar=c(3,3,1,1))
plot(index(corr_length), corr_length, axes=FALSE, xlim=c(1,length(corr_length))+c(-0.5,0.5), xlab="", ylab="",
     ylim = c(2,4.5), pch=PCHS, col=COL)
arrows(index(corr_length), corr_length-corr_length.sd,
       index(corr_length), corr_length+corr_length.sd, col=COL, length=0.05, angle=90, code=3, lwd=2)
plot_axis(2, label="Correlation length [years]")
abline(h=axTicks(2), col = "lightgray", lty = "dotted")
plot_letter("b")
box()

par(new=T)
layout(1,1,1)
par(mar=c(1,1,1,4))
plot(0, axes=FALSE, xlab="", ylab="", type="n")
legend("bottomright", 
       legend=names(states)[order], 
       col=COL[order], 
       lty=rep(1,4), 
       pch=PCHS[order], lwd=2, y.intersp=1, x.intersp=1, inset=c(0.,0.21))

if(plot) dev.off()
