require(boot)
require(parallel)
require(graphics)

source("colours.R")
source("functions.R")
load("run_dictionary.RData")
load("ancil.RData")

#---DATA---#
booting <- T
data.dir <- "data/HadCM3/"
source("02_surface_climate_response/data_anomaly.R")

seaicePolys <- list()
for (run_type in run_types)
{
  all_data$`Sea Ice`[[run_type]][,143:144] <- 1.0
  seaicePolys[[run_type]] <- generate_polygon_area(all_data$`Sea Ice`[[run_type]] >= 0.5)
}

iceshieldPolys <- list()
iceshieldPolys[["LGM*"]] <- generate_polygon_area(ancil$iceshield$lgm == 1)
iceshieldPolys[["PI*"]] <- generate_polygon_area(ancil$iceshield$pi == 1)

#---PLOT---#
nbreakpoints <- 17
latgrid <- c(-91.25,-60,-30,0,30,60,91.25)
zlim <- list("Precipitation"=c(-3.4,3.4),
             "Sea Level Pressure"=c(-3.4,3.4), 
             "500mbar Wind Fields"=c(-3.4,3.4))
col <- colorRampPalette(rev(c("black","#962D26","#de2d26","#fc9272","#fee0d2","white",
                              "#deebf7","#3182bd","#7fcdbb","#edf8b1",
                              "#ffc200")))

plot <- F

if(plot)
{
  pdf(file=paste("02_surface_climate_response/048_stdanom_matrix_new.pdf",sep=""),
      width=11,height=11)
}

layout(matrix(c(1:9,c(10,10,0)),4,3,byrow=T),widths=c(1,1,0.5),heights=c(1,1,1,0.15))
par(oma=c(3,4,0,1))

lineCol <- list("LGM*"=COLS[["LGM"]], "PI*"=COLS[["PI"]])

letters <- c("a","b","c","d","e","f","g","h","i")
i.letter <- 0

upper_mar <- 4
lower_mar <- 2

vars <- c("Temperature", "Precipitation","Sea Level Pressure", "500mbar Wind Fields")

for (i in 2:4)
{
  var <- vars[i]
  
  par(mar=c(lower_mar,0,upper_mar,0.2))
  for (run_type in run_types)
  {
    if(run_type=="PI*") par(mar=c(lower_mar,0.2,upper_mar,0))
    
    if(i == 4)
    {
      fld_u <- all_data$`Westerly Winds`$std.anom[[run_type]]
      fld_v <- all_data$`Southerly Winds`$std.anom[[run_type]]
      plot_windfield(fld_u,fld_v,
                     stretch=7,#0.2
                     strength.cutoff=.1,#2
                     zlim=zlim[[var]],
                     col=col(nbreakpoints),
                     by=3,
                     arrowhead.size=0.8)
    } else {
      plot_matrix(all_data[[var]]$std.anom[[run_type]], zlim[[var]], col=col(nbreakpoints))
    }
    add_map(lgm=(run_type=="LGM*"))
    
    if(run_type=="LGM*")
    {
      plot_axis(2,at=latgrid, ticks=c("90°S","60°S","30°S","0","30°N","60°N","90°N"), 
                lineVal = 0.75, lineLab = 2, lab="Latitude")
    }
    abline(h=latgrid, col="grey", lty=2, lwd=1)
    
    add_polygons(iceshieldPolys[[run_type]],T,col=rgb(0, 0, 0, alpha = 1, maxColorValue = 8))
    add_polygons(seaicePolys[[run_type]],F,col=rgb(0, 0, 0, alpha = 1, maxColorValue = 8))
    if(i!=4)
    {
      add_dots(cond=((all_data[[var]]$std.anom[[run_type]] > all_data[[var]]$`anom-0.99`$lower)
                     & (all_data[[var]]$std.anom[[run_type]] < all_data[[var]]$`anom-0.99`$upper)),
               col=rgb(t(col2rgb("black")), alpha = 100, maxColorValue = 255))
    }
    
    legend("topright", legend=paste(run_type,"    ",sep=""), text.font=2, lwd=2, col=lineCol[[run_type]], bg="white")
    
    box()
    par(cex=1)
    plot_letter(letters[i.letter <- i.letter+1])
    par(cex=0.66)
  }
  
  mtext(vars[i], side=3, line=0.2, at=par("usr")[1], font=2)
  
  par(mar=c(lower_mar,0.4,upper_mar,1))
  if(i==1){xlim <- c(-2,0)}else{xlim <- range(-0.8,0.8)}
  if(i == 4)
  {
    var <- "Westerly Winds"
    nlats <- 72
  } else {
    nlats <- 73
  }
  
  plot(rev(apply(all_data[[var]]$std.anom[[run_types[1]]], 2, mean, na.rm=T)), seq(-91.25,91.25,length.out=nlats),
       type="l",axes=F,yaxs="i",xlim=xlim,col=lineCol[[run_types[1]]],lwd=2)
  lines(rev(apply(all_data[[var]]$std.anom[[run_types[2]]], 2, mean, na.rm=T)), seq(-91.25,91.25,length.out=nlats),
        type="l",col=lineCol[[run_types[2]]],lwd=2)
  lines(rev(apply(all_data[[var]]$std.anom[[run_types[2]]], 2, mean, na.rm=T)), seq(-91.25,91.25,length.out=nlats),
        type="l",col=lineCol[[run_types[2]]],lwd=2)
  
  abline(v=0, col="grey")
  abline(v=c(-2,-1.5,-1,-0.5,0.5,1), lty=2, col="grey")
  abline(h=latgrid, col="grey", lty=2, lwd=1)
  
  plot_axis(1, lineVal=1, lineLab=2.2, lab="MSA")
  
  if (i==2)
  {
    plot_axis(3, lab="AOD", at=seq(-0.5,0.5,0.5), ticks=c(0.28, 0.23, 0.18), col="black", lineLab=1.7)
    par(new=T)
    AOD <-  (apply(D[which(volc_mean > 0.13),],2,mean)-0.23)*(-10)
    AOD <- rep(c(rep(AOD[1],2),AOD[2],AOD[3],rep(AOD[4],2)),each=200)
    plot(AOD,1:1200, type="l", lty=2, lwd=1.5, xlim=xlim,axes=F,yaxs="i",col="black")
  }
  
  box()
  par(cex=1)
  plot_letter(letters[i.letter <- i.letter+1])
  par(cex=0.66)
}

par(mar=c(1,1,1,1))
plot_colorbar(c(-3,3), col=col(nbreakpoints)[2:(length(col(nbreakpoints))-1)], axis.pos=1, add.axis=F)
plot_axis(1,at=axTicks(1), lineVal=1, lineLab=2, lab="MSA")

if(plot) dev.off()
