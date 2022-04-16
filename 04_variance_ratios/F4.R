source("colours.R")
source("functions.R")
load("run_dictionary.RData")
load("04_variance_ratios/proxylist.RData")
data.dir = "data/HadCM3/" #The directory were the HadCM3 data is located. Please run the corresponding `.sh` scripts first.

source("04_variance_ratios/proxy_var_ratio_average.R")

saveToPDF <- F

if(saveToPDF) 
  {
   cairo_pdf("04_variance_ratios/060_proxy_var_ratio_bilint_rev_revised.pdf",
              width=11,height=4)
  }
  
  xticks <- seq(-90,90,30)
  xlabels <- c("","60°S","","0°","","60°N","")
  yticks <- 10^(-3:3)
  letters <- c("a","b","c","d","e")
  letters <- c("a","b","c","d")
  i.letter <- 0
  
  layout(matrix(c(1:4,rep(5,4)),2,4,byrow=T), widths=c(1,1,1,1), heights=c(1,0.1))
  par(mar=c(2,0,0,0), oma=c(1,3,.5,0),cex=1) #bottom, left, top, right
  
  archives <- unique(unlist(sapply(QC.all, function(x) {unique(sapply(x, function(y) y$metafilt$Archive))}))) 
  PCH <- seq_len(length.out=length(archives))
  names(PCH) <- c("ice core","marine sediment","observation",
                  "tree","lake sediment","hybrid","documents")
  order.legend <- c(3,2,6,4,1,5)
  
  #plot.state <- c(rep("PI",length(tsc.in)),"LGM")
  plot.state <- rep("PI",length(tsc.in))
  plot.tsc <- c(names(tsc.in), names(tsc.in)[length(names(tsc.in))])
 # for (i in 1:(length(tsc.in)+1))
  for (i in 1:(length(tsc.in)))
  {
    state <- plot.state[i]
    tsc <- plot.tsc[i]
    state.f <- paste0(state, "*")
    
    plot(1,type="n",xlim=c(-100,100),ylim=c(0.001,1000),
         axes=F,xlab="",ylab="",log="y", xaxs="i",
         panel.first = 
           {
             abline(h=yticks,v=xticks,col="lightgrey");
             abline(h=1,lwd=2);
           })
    if(!is.null(var.ratio[[state]][[tsc]]) & !is.null(var.ratio[[state.f]][[tsc]]))
    {
      PCH.archive <- PCH[QC.all[[state]][[tsc]]$metafilt$Archive]
      lats <- QC.all[[state]][[tsc]]$metafilt$Lat
      points(lats, var.ratio[[state]][[tsc]],pch=PCH.archive,col=adjustcolor("black",0.5),cex=.8)
      arrows(lats, var.ratio[[state]][[tsc]]+var.ratio.CI[[state]][[tsc]][,2],
             lats, var.ratio[[state]][[tsc]]-var.ratio.CI[[state]][[tsc]][,1],
             col=adjustcolor("black",0.4), length=0, code=3,lwd=2)
      points(lats, var.ratio[[state.f]][[tsc]],pch=PCH.archive,col=pal[[1]],cex=.8)
      arrows(lats, var.ratio[[state.f]][[tsc]]+var.ratio.CI[[state.f]][[tsc]][,2],
             lats, var.ratio[[state.f]][[tsc]]-var.ratio.CI[[state.f]][[tsc]][,1],
             col=adjustcolor(pal[[1]], 0.7), length=0, code=3)
    }
    if(tsc == names(tsc.in)[1]) 
    {
      plot_axis(2, at=yticks, ticks=yticks, label="Variance ratio r(sim./obs.)")
      text(-90,300,"unforced simulations",adj=0,cex=.95)
      text(-90,130,"naturally forced simulations",adj=0,cex=.95, col=pal[[1]])
    }
    plot_axis(1, at=xticks, ticks=xlabels, label="Latitude")
    box()
    text(-80,740,paste0(tsc,"yrs"),adj=0,cex=0.8)
    #if(state=="LGM"){text(95,780,state,adj=1,cex=0.8, col=COLS[[state]])} else {text(95,780,state,adj=1,cex=0.8)}
    text(-94, 0.001,substitute(paste("f=", c, "(", u, ",", l, ")"), list(c=format(round(changes[[state]][[tsc]], 2), nsmall=2),
                                                                        l=format(round(lower.bounds[[state]][[tsc]],2), nsmall = 2),
                                                                        u=format(round(upper.bounds[[state]][[tsc]],2), nsmall = 2)
                                                                          )), pos=4, cex=0.8)
    plot_letter(letters[i.letter <- i.letter+1])
  }
  
  par(mar=c(0,0,0,0))
  plot(c(0,0),c(1,1),type="n",axes=F,xlab="",ylab="")
  legend("top"
         ,pch=PCH[archives][order.legend]
         ,legend=archives[order.legend]
         ,xpd=NA
         ,inset=.75,
         ,xjust=0.5
         ,horiz=T
         ,seg.len=0
         ,cex=0.6)
  
if(saveToPDF) dev.off()
