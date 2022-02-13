source("colours.R")
source("functions.R")
load("run_dictionary.RData")
load("04_variance_ratios/proxylist.RData")
data.dir = "data/HadCM3/" #The directory were the HadCM3 data is located. Please run the corresponding `.sh` scripts first.

#---DATA---#
source("04_variance_ratios/proxy_var_ratio_average.R")

#---FUNCTIONS---#
llines2 <- function (x, conf = TRUE, col = "black", alpha = 0.3, removeFirst = 0, 
          removeLast = 0, ...) 
{
  index <- (removeFirst + 1):(length(x$freq) - removeLast)
  x$freq <- x$freq[index]
  x$spec <- x$spec[index]
  x$lim.1 <- x$lim.1[index]
  x$lim.2 <- x$lim.2[index]
  if (conf) 
    polygon(c(x$freq, rev(x$freq)), c(x$lim.1, rev(x$lim.2)), 
            col = ColTransparent(col, alpha), border = NA)
}

log10.axis <- function(side, at, ...) {
  at.minor <- log10(outer(1:2, 10^(min(at):max(at))))
  #lab <- sapply(at, function(i) as.expression(bquote(10^ .(i))))
  axis(side=side, at=at.minor, labels=NA, tcl=par("tcl")*0.5, ...)
  axis(side=side, at=at, labels=F, ...)
}

minor.ticks.axis <- function(ax,n,t.ratio=0.5,mn,mx,...){
  
  lims <- par("usr")
  if(ax %in%c(1,3)) lims <- lims[1:2] else lims[3:4]
  
  major.ticks <- pretty(lims,n=5)
  if(missing(mn)) mn <- min(major.ticks)
  if(missing(mx)) mx <- max(major.ticks)
  
  major.ticks <- major.ticks[major.ticks >= mn & major.ticks <= mx]
  
  labels <- sapply(major.ticks,function(i)
    as.expression(bquote(10^ .(i)))
  )
  axis(ax,at=major.ticks,labels=F,...)
  
  n <- n+2
  minors <- log10(pretty(10^major.ticks[1:2],n))-major.ticks[1]
  minors <- minors[-c(1,n)]
  
  minor.ticks = c(outer(minors,major.ticks,`+`))
  minor.ticks <- minor.ticks[minor.ticks > mn & minor.ticks < mx]
  
  
  axis(ax,at=minor.ticks,tcl=par("tcl")*t.ratio,labels=FALSE)
}

plot_all_spec_incl <- function(allvals1,allvals2,allvals3,ids, df.log=0.01, mar=0.5,
                               Col.ts,Lty.ts,xlimz,tsc.in,varratiomean=NULL,var.ratio.CI=NULL){
  
  f <- 1/tsc.in
  spec1 <- lapply(allvals1, function(x) list(spec=x$spec, f=x$freq, dof=x$dof))
  spec2 <- lapply(allvals2, function(x) list(spec=x$spec, f=x$freq, dof=x$dof))
  spec3 <- lapply(allvals3, function(x) list(spec=x$spec, f=x$freq, dof=x$dof))
  
  for (i in 1:length(spec2)) {
    ylimz<-c(min(min(spec1[[i]]$spec[which(spec1[[i]]$f <= f[[1]] & spec1[[i]]$f >= f[[2]])]),
                 min(spec2[[i]]$spec[which(spec2[[i]]$f <= f[[1]] & spec2[[i]]$f >= f[[2]])]),
                 min(spec3[[i]]$spec[which(spec3[[i]]$f <= f[[1]] & spec3[[i]]$f >= f[[2]])])),
             max(max(spec1[[i]]$spec[which(spec1[[i]]$f <= f[[1]] & spec1[[i]]$f >= f[[2]])]),
                 max(spec2[[i]]$spec[which(spec2[[i]]$f <= f[[1]] & spec2[[i]]$f >= f[[2]])]),
                 max(spec3[[i]]$spec[which(spec3[[i]]$f <= f[[1]] & spec3[[i]]$f >= f[[2]])])))
    plot(xlimz,c(ylimz[[1]]-mar*ylimz[[1]],ylimz[[2]]+mar*ylimz[[2]]),axes=FALSE,log="xy",type="n",xlab="",ylab="")
    
    if (!is.null(spec1[[i]])) LLines((AddConfInterval(LogSmooth(allvals1[[i]], df.log=df.log, removeLast = 0))), col=Col.ts[1],lty=Lty.ts[1],lwd=1)
    if (!is.null(spec2[[i]])){llines2((AddConfInterval(LogSmooth(allvals2[[i]], df.log=df.log, removeLast = 0))), lwd=0)
                              LLines((AddConfInterval(LogSmooth(allvals2[[i]], df.log=df.log,  removeLast = 0))), col=Col.ts[2],lty=Lty.ts[2],lwd=2)}
    if (!is.null(spec3[[i]])){llines2((AddConfInterval(LogSmooth(allvals3[[i]], df.log=df.log, removeLast = 0))), lwd=0)
                              LLines((AddConfInterval(LogSmooth(allvals3[[i]], df.log=df.log, removeLast = 0))), col=Col.ts[3],lty=Lty.ts[3],lwd=2)}
    
    axis(2,line=-1,tick=F)
    #log10.axis(2, at=seq(0, 3000, 500))
    #minor.ticks.axis(2,1, mn=300, mx=10000)
    axis(1,line=-1,tick=F,labels=tsc.in,at=f)
    #minor.ticks.axis(1,9)
    #log10.axis(1, at=seq(0, 3, 1))
    box()
    abline(v=f,col="grey")
    
    title(ids[i])
    
    rect(axTicks(1)[1]/10,ylimz[1]/10,f[2],ylimz[2]*10,col=adjustcolor("grey",0.25),border=NA)
    rect(f[1],ylimz[1]/10,xlimz[2]*10,ylimz[2]*10,col=adjustcolor("grey",0.25),border=NA)
    }
}

#---PLOT---#

col.blue <- COLS[["LGM"]]
col.orange <- COLS[["PI"]]

LTYS <- c(1,1,2)

saveToPDF <- T
if(saveToPDF) 
{
  cairo_pdf("04_variance_ratios/063_all_spectra_average_short_new.pdf",
            width=12,height=7, bg="transparent")
}

df.log <- 0.02
par(mfrow=c(6,11),mar=c(1,1,1,.5),oma=c(1,1,1,1),las=0,lwd=1)

COL <- c(pal[[5]],col.orange,col.orange)
state <- "PI"
tsc <- names(tsc.in)[1]
plot_all_spec_incl(spec.proxy[[state]][[tsc]],
                     meanspec.hadcm[[paste0(state,"*")]][[tsc]],
                     meanspec.hadcm[[state]][[tsc]],
                     QC.all[[state]][[tsc]]$metafilt$ID,
                     Col.ts=COL,Lty.ts=LTYS,
                     mar=0.2,
                     df.log=df.log,
                     xlimz=c(1/(1.2*max(tsc.in[[1]])), 1/(0.8*min(tsc.in[[1]]))),tsc.in[[tsc]])
tsc <- names(tsc.in)[2]
plot_all_spec_incl(spec.proxy[[state]][[tsc]],
                   meanspec.hadcm[[paste0(state,"*")]][[tsc]],
                   meanspec.hadcm[[state]][[tsc]],
                   QC.all[[state]][[tsc]]$metafilt$ID,
                   Col.ts=COL,Lty.ts=LTYS,
                   mar=0.2,
                   df.log=df.log,
                   xlimz=c(1/(1.2*max(tsc.in[[2]])), 1/(0.8*min(tsc.in[[2]]))),tsc.in[[tsc]])

if(saveToPDF) dev.off()

if(saveToPDF) 
{
  cairo_pdf("04_variance_ratios/063_all_spectra_average_long_new.pdf",
            width=12,height=14, bg="transparent")
}

par(mfrow=c(13,12),mar=c(1,1,1,.5),oma=c(1,1,1,1),las=0,lwd=1)
df.log <- 0.02

COL <- c(pal[[5]],col.orange,col.orange)
state <- "PI"
tsc <- names(tsc.in)[3]
plot_all_spec_incl(spec.proxy[[state]][[tsc]],
                   meanspec.hadcm[[paste0(state,"*")]][[tsc]],
                   meanspec.hadcm[[state]][[tsc]],
                   QC.all[[state]][[tsc]]$metafilt$ID,
                   Col.ts=COL,Lty.ts=LTYS,
                   mar=0.8,
                   df.log=df.log,
                   xlimz=c(1/(1.2*max(tsc.in[[3]])), 1/(0.8*min(tsc.in[[3]]))),tsc.in[[tsc]])
tsc <- names(tsc.in)[4]
plot_all_spec_incl(spec.proxy[[state]][[tsc]],
                   meanspec.hadcm[[paste0(state,"*")]][[tsc]],
                   meanspec.hadcm[[state]][[tsc]],
                   QC.all[[state]][[tsc]]$metafilt$ID,
                   Col.ts=COL,Lty.ts=LTYS,
                   mar=0.8,
                   df.log=df.log,
                   xlimz=c(1/(1.2*max(tsc.in[[4]])), 1/(0.8*min(tsc.in[[4]]))),tsc.in[[tsc]])
state <- "LGM"
COL <- c(pal[[5]],col.blue,col.blue)
plot_all_spec_incl(spec.proxy[[state]][[tsc]],
                   meanspec.hadcm[[paste0(state,"*")]][[tsc]],
                   meanspec.hadcm[[state]][[tsc]],
                   QC.all[[state]][[tsc]]$metafilt$ID,
                   Col.ts=COL,Lty.ts=LTYS,
                   mar=.8,
                   df.log=df.log,
                   xlimz=c(1/(1.2*max(tsc.in[[4]])), 1/(0.8*min(tsc.in[[4]]))),tsc.in[[tsc]])

if(saveToPDF) dev.off()

all_IDs <- unique(unlist(lapply(QC.all, function(x) lapply(x, function(y) y$metafilt$ID))))
obs_length <- length(QC.all[["PI"]][["2-5"]][["metafilt"]][["ID"]])
proxy_length <- length(all_IDs) - obs_length
print(paste0(proxy_length, " proxies and ", obs_length, " observations used."))
