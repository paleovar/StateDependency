source("colours.R")
source("functions.R")
load("run_dictionary.RData")
load("ancil.RData")
data.dir <- "data/HadCM3/"

#---FUNCTIONS---#
ncvarnames <- c("temp_1")
vars <- c("Surface Temperature")
i.var = 1
tscale1 <- 500 #the bigger one
tscale2 <- 50 #the smaller one

simpleawmean<-function(fld, lats = seq(from = -90, to = 90, length.out = 73)){
  #https://github.com/paleovar/iHadCM3LastMill/blob/master/Functions/aw_mean.R
  zonmean<-apply(fld,2,mean,na.rm=TRUE)
  w.lats<-cos(lats*pi/180)/sum(cos(lats*pi/180))
  tavg<-sum(w.lats*zonmean,na.rm=TRUE)
  return(tavg)
}

simpleawsd<-function(fld, lats = seq(from = -90, to = 90, length.out = 73)){
  #https://github.com/paleovar/iHadCM3LastMill/blob/master/Functions/aw_mean.R
  zonmean<-apply(fld,2,sd,na.rm=TRUE)
  w.lats<-cos(lats*pi/180)/sum(cos(lats*pi/180))
  tavg<-sum(w.lats*zonmean,na.rm=TRUE)
  return(tavg)
}

#---DATA---#
# Set global zlim
zlim <- array(NA,1) #temp only
mini <- array(NA,1) #temp only
source("03_spectral_analysis/data_world.R")
zlim[i.var] <- max(abs(sapply(scaling, function(x) range(-1*x, na.rm=T))))
mini[i.var] <- min(sapply(scaling, function(x) range(-1*x, na.rm=T)))

means <- list()
sd <- list()
for(run_type in c("LGM", "LGM*", "PI", "PI*")){
  means[[run_type]] <- simpleawmean(scaling[[run_type]]*(-1))
  sd[[run_type]] <- simpleawsd(scaling[[run_type]]*(-1))
}

#---PLOT---#
# Set up plot
zlim <- round(max(zlim)*c(-1,1),1)
mini <-  round(min(mini),1)
zlim <- c(-2.2, 2.2)
col <- colorRampPalette(c(rev(brewer.pal(8,"Blues")),"#FFFFFF",brewer.pal(8,"Reds")))
nbreakpoints=11 #or 9 depending on timescale

saveToPDF <- F

if(saveToPDF){
    cairo_pdf(paste0("03_spectral_analysis/043_spectral_analysis_world", "_", tscale1, "_", tscale2,".pdf"),
              width=5,height=4)
}

layout(matrix(c(1:4,5,5),3,2,byrow=T), widths=c(1,1), heights=c(1,1,0.2))
par(oma=c(2,1.5,1.5,1), cex=0.6)
letters <- c("a","b","c","d")
i.letter <- 0

for (i.var in c(1)) 
{
  var <- vars[i.var]
  
  par(mar=c(1,.5,0.9,.5))
  
  for (run_type in c("LGM", "LGM*", "PI", "PI*"))
  {
    plot_matrix(scaling[[run_type]]*(-1), zlim, col=col(nbreakpoints))
    add_map(lgm=c(F,T)[grepl("LGM",run_type)+1])
    box()
    plot_letter(letters[i.letter <- i.letter+1],cex = 0.8)
    text(1.07, 0.6, substitute(paste(hat(beta), "=", m, ''%+-%'', s), list(m=format(round(means[[run_type]], 2), nsmall=2),
                                                                       s=format(round(sd[[run_type]],2), nsmall = 2))
                               ), adj=0, pos=4, cex=1.2)
    title(run_type, line=0.15)
  }
  par(mar=c(1.2,10,0.7,10))
  offset <- 3
  plot_colorbar(seq(zlim[1],zlim[2], length.out = nbreakpoints+1)[c((1+offset),(nbreakpoints+1))], 
                  col=col(nbreakpoints)[(1+offset):(nbreakpoints)], axis.pos=1, add.axis=F, cex=0.15)
  plot_axis(1, lineVal=0.8, lineLab=2.5, 
              lab=substitute(paste(nn, beta), list(nn="Spectral exponent ")), font=0.3)
}

if(saveToPDF) dev.off()
