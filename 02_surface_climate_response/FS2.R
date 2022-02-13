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
zlim <- list("Surface Temperature"=c(-80,35),
             "Precipitation"=c(0,7000),
             "dO18"=c(-100,0), 
             "Sea Level Pressure"=c(970,1040), 
             "500mbar Wind Fields"=c(0,30))
col <- list("Surface Temperature"=colorRampPalette(rev(brewer.pal(9,"YlGnBu"))),
            "Precipitation"=colorRampPalette(c(rep("white",0),brewer.pal(6,"GnBu"),"#6C3CE8","#6C3CE8")),
            "Sea Level Pressure"=colorRampPalette(IPCCColPal$BlPu),
            "500mbar Wind Fields"=colorRampPalette(IPCCColPal$YlRd))
colorbarlab <- list("Surface Temperature"="Mean surface temperature [°C]",
                    "Precipitation"="Mean precipitation [mm/year]",
                    "Sea Level Pressure"="Mean sea level pressure [hPa]",
                    "500mbar Wind Fields"="Mean wind velocity (500mbar) [m/s]")
zonalmeanlab <- list("Surface Temperature"="Zon. mean surface temperature [°C]",
                    "Precipitation"="Zon. mean precipitation [mm/year]",
                    "Sea Level Pressure"="Zon. mean sea level pressure [hPa]",
                    "500mbar Wind Fields"="Zon. mean wind vel. (500mbar) [m/s]")

plot <- F

if(plot)
{
  cairo_pdf(file="02_surface_climate_response/048_mean_matrix.pdf",
            width=11,height=15)
}

layout(matrix(c(1:3,c(4,4,0),5:7,c(8,8,0),9:11,c(12,12,0),13:15,c(16,16,0)),8,3,byrow=T),widths=c(1,1,0.5),heights=c(rep(c(1,0.1),4)))
par(oma=c(3,4,0,2))

lineCol <-  list("LGM*"=COLS[["LGM"]], "PI*"=COLS[["PI"]])

upper_mar <- 4
lower_mar <- 0.5

letters <- c("a","b","c","d","e","f","g","h","i","j","k","l")
i.letter <- 0

vars <- c("Surface Temperature","Precipitation","Sea Level Pressure", "500mbar Wind Fields")
units <- c("°C","mm/year","hPa","m/s")
sign_digits <- c(1,0,1,1)
for (i in 1:4)
{
  var <- vars[i]
  
  par(mar=c(lower_mar,0,upper_mar,0.2))
  for (run_type in run_types)
  {
    if(run_type=="PI*") par(mar=c(lower_mar,0.2,upper_mar,0))
    
    if(grepl("Wind",var))
    {
      fld_u <- all_data$`Westerly Winds`$mean[[run_type]]
      fld_v <- all_data$`Southerly Winds`$mean[[run_type]]
      plot_windfield(fld_u,fld_v,
                     stretch=0.4,
                     strength.cutoff=1.5,
                     zlim=zlim[[var]],
                     col=col$`500mbar Wind Fields`(nbreakpoints),
                     by=3,
                     arrowhead.size=0.8)
      fld_u.sd <- all_data$`Westerly Winds`$sd[[run_type]]
      fld_v.sd <- all_data$`Southerly Winds`$sd[[run_type]]
      fld.sd <- sqrt(fld_u.sd*fld_u.sd+fld_v.sd*fld_v.sd)[,c(1:72,72)]
      fld.sd[fld.sd == Inf | fld.sd == -Inf] <- NA
      fld <- sqrt(fld_u*fld_u+fld_v*fld_v)[,c(1:72,72)]
      fld[fld == Inf | fld == -Inf] <- NA
      meantext <- paste("Global mean: (",round(fldmean(fld, na.rm=T),sign_digits[i]),
                        "±",round(fldmean(fld.sd, na.rm=T),sign_digits[i]),") ",units[i],sep="")
    } else {
      plot_matrix(all_data[[var]]$mean[[run_type]], zlim[[var]], col=col[[var]](nbreakpoints))
      meantext <- paste("Global mean: (",round(fldmean(all_data[[var]]$mean[[run_type]], na.rm=T),sign_digits[i]),
                        "±",round(fldmean(all_data[[var]]$sd[[run_type]], na.rm=T),sign_digits[i]),") ",units[i],
                        sep="")
    }
    
    if(run_type == "LGM*")
    {
      mtext(meantext, side=3, at=par("usr")[1],adj=0, cex=0.8)
    } else {
      mtext(meantext, side=3, at=par("usr")[2],adj=1, cex=0.8)
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
   
    legend("topright", legend=paste(run_type,"    ",sep=""), text.font=2, lwd=2, col=lineCol[[run_type]], bg="white")
    
    box()
    par(cex=1)
    plot_letter(letters[i.letter <- i.letter+1])
    par(cex=0.66)
  }
  
  mtext(vars[i], side=3, line=0.2, at=par("usr")[1], font=2)
  
  par(mar=c(lower_mar,0.4,upper_mar,1))
  if(grepl("Wind",var))
  {
    var <- "Westerly Winds"
    nlats <- 72
  } else {
    nlats <- 73
  }
  fld.lgm <- rev(apply(all_data[[var]]$mean[[run_types[1]]], 2, mean, na.rm=T))
  fld.pi <- rev(apply(all_data[[var]]$mean[[run_types[2]]], 2, mean, na.rm=T))
  xlim <- range(c(fld.lgm[fld.lgm != Inf],fld.pi[fld.pi != Inf]), na.rm=T)
  plot(fld.lgm, seq(-91.25,91.25,length.out=nlats),
       type="l",axes=F,yaxs="i",xlim=xlim,col=lineCol[[run_types[1]]],lwd=2)
  lines(fld.pi, seq(-91.25,91.25,length.out=nlats),
        type="l",col=lineCol[[run_types[2]]],lwd=2)
  abline(v=0, col="grey")
  abline(v=axTicks(1), lty=2, col="grey")
  abline(h=latgrid, col="grey", lty=2, lwd=1)
  
  
  if (grepl("Wind",var)) var <- "500mbar Wind Fields"
  
  if (i %in% c(2,3)) xTicks <- axTicks(1)[seq(1,length(axTicks(1)),2)]
  else xTicks <- axTicks(1)
  plot_axis(1, at=xTicks, ticks=xTicks, lineVal=1, lineLab=2.2, lab=zonalmeanlab[[var]])
  box()
  
  par(cex=1)
  plot_letter(letters[i.letter <- i.letter+1])
  par(cex=0.66)
  
  par(mar=c(1.2,6,0,6))
  plot_colorbar(zlim[[var]], col=col[[var]](nbreakpoints), axis.pos=1, add.axis=F)
  plot_axis(1,at=axTicks(1), lineVal=1, lineLab=2.1, lab=colorbarlab[[var]])
}

if(plot) dev.off()