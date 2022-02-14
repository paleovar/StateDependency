source("colours.R")
source("functions.R")
load("run_dictionary.RData")
load("ancil.RData")
data.dir <- "data/HadCM3/" #Generate fldmean prior to this analysis by using the `getFldmean.sh` script in `/HadCM3/`

#---FUNCTIONS---#
des <- c(0.5,rep(1,11),0.5)/12

dist_vals <- function(mat){
  t <- c(mat)
  t <- t[!is.na(t)]
  return(t)
}

#---DATA---#
filename <- "01_characteristics/057_state_dist"

runs <- names(run_dict)
LGM <- unlist(run_dict, use.names = FALSE) == "LGM" | unlist(run_dict, use.names = FALSE) == "LGM*"
forced <- unlist(run_dict, use.names = FALSE) == "LGM*" | unlist(run_dict, use.names = FALSE) == "PI*"

if(!file.exists(paste(filename,".RData",sep="")))
{
  gmst <- list()
  for (run in runs) {
    nc <- nc_open(paste(data.dir, "surface_temperature/",run,".nc",sep=""))  #alternatively compute fldmean using `./data/getFldmean.sh`
    gmst[[run]] <- nc_var_to_TS(nc, "temp_1")
    gmst[[run]]$data <- fldmean(gmst[[run]]$data, na.rm=T)
    nc_close(nc)
  }
  
  timesteps <- sort(as.integer(unique(unlist(sapply(gmst, function(x) x$time)))))
  n_timesteps <- length(timesteps)
  timesteps.range <- range(timesteps)
  gmst_sync <- matrix(NA,nrow=n_timesteps,ncol=length(runs))
  
  for (i in 1:length(runs)) {
    ind <- which(timesteps %in% gmst[[i]]$time)
    # We can do this since timesteps and thus the first dimension has _all_ timesteps
    gmst_sync[ind,i] <- gmst[[i]]$data
  }
  
  save("data","timesteps","gmst_sync",file=paste(filename,".RData",sep=""))
} else {
  load(paste(filename,".RData",sep=""))
}

#---PLOT---#
saveToPDF <- T

gmst_sync_f_d <- apply(gmst_sync, 2, function(x){detrend(filter(x, des))})

var <- list()
sd <- list()
variances <- apply(gmst_sync_f_d, 2, function(x){var(x, na.rm = T)})
variances <- sqrt(variances)
var$lgm_uf <- mean(variances[LGM&!forced])
var$lgm_f <- mean(variances[LGM&forced])
var$pi_uf <- mean(variances[!LGM&!forced])
var$pi_f <- mean(variances[!LGM&forced])

sd$lgm_uf <- sd(variances[LGM&!forced])
sd$lgm_f <- sd(variances[LGM&forced])
sd$pi_uf <- sd(variances[!LGM&!forced])
sd$pi_f <- sd(variances[!LGM&forced])

if(saveToPDF) cairo_pdf(filename=paste0(filename, ".pdf"), width=4, height=3)

par(mai = c(.5, .5, .1, .1), cex=1)
xlimz <- range(gmst_sync_f_d, na.rm=T)
print(xlimz)
xlimz=round(xlimz,1)
ylimz <- c(0,4)

dens <- 12
greys_lgm <- COLS_a[["LGM"]]
greys_pi <- COLS_a[["PI"]]

hist(dist_vals(gmst_sync_f_d[,LGM&!forced]), density = dens, angle = 30, col = greys_lgm, border = greys_lgm, 
     breaks=seq(-1.4, 0.6, by=0.05), prob=T, ylim=ylimz, xlab="", ylab="", yaxt='n', xaxt='n',
     xaxs = "i", yaxs = "i", main="", panel.first={grid(ny=4, nx=0)}) 
abline(v = seq(-1.5, 0.5, by=0.5),
       lty = 3, col = "grey80")

hist(dist_vals(gmst_sync_f_d[,LGM&forced]), breaks=seq(-1.4, 0.6, by=0.05), col = greys_lgm, border = greys_lgm,  prob=T, add=T)
hist(dist_vals(gmst_sync_f_d[,LGM&forced]),  density = dens, angle = -30, breaks=seq(-1.4, 0.6, by=0.05), col = greys_lgm, border = greys_lgm,  prob=T, add=T)

hist(dist_vals(gmst_sync_f_d[,!LGM&!forced]),  density = dens, angle = 30, breaks=seq(-1.4, 0.6, by=0.05), col=greys_pi, border=greys_pi, prob=T, add=T)

hist(dist_vals(gmst_sync_f_d[,!LGM&forced]), breaks=seq(-1.4, 0.6, by=0.05), col=greys_pi, border=greys_pi, prob=T, add=T)
hist(dist_vals(gmst_sync_f_d[,!LGM&forced]),  density = dens, angle = -30, breaks=seq(-1.4, 0.6, by=0.05), col=greys_pi, border=greys_pi, prob=T, add=T)

lines(density(dist_vals(gmst_sync_f_d[,!LGM&!forced])), col = COLS[["PI"]], lwd=2, lty=2)   
lines(density(dist_vals(gmst_sync_f_d[,!LGM&forced])), col=COLS[["PI"]], lwd=2) 
lines(density(dist_vals(gmst_sync_f_d[,LGM&!forced])),  xlim=c(-1.4, 0.6), col = COLS[["LGM"]], lwd=2, lty=2)  
lines(density(dist_vals(gmst_sync_f_d[,LGM&forced])), col=COLS[["LGM"]], lwd=2)  

start=1.9
shift=0.29
start.x=-1.4

#Uncertainty
delta_var_lgm_f_lgm_uf = (var$lgm_f/var$lgm_uf)*sqrt((sd$lgm_f/var$lgm_f)**2 + (sd$lgm_uf/var$lgm_uf)**2)
delta_var_pi_f_pi_uf = (var$pi_f/var$pi_uf)*sqrt((sd$pi_f/var$pi_f)**2 + (sd$pi_uf/var$pi_uf)**2)
delta_var_lgm_f_pi_f = (var$lgm_f/var$pi_f)*sqrt((sd$lgm_f/var$lgm_f)**2 + (sd$pi_f/var$pi_f)**2)
delta_var_lgm_uf_pi_uf = (var$lgm_uf/var$pi_uf)*sqrt((sd$lgm_uf/var$lgm_uf)**2 + (sd$pi_uf/var$pi_uf)**2)

text(start.x, start, "forcing dependence", pos=4, cex=.8)
text(start.x, start -shift, substitute(paste(r[sigma], "(LGM* / LGM) =", a, "±", b), list(a=round(var$lgm_f/var$lgm_uf,2), b=round(delta_var_lgm_f_lgm_uf,2))), cex=.8, pos=4)
text(start.x, start -2*shift , substitute(paste(r[sigma], "(PI* / PI) =", a, "±", b), list(a=round(var$pi_f/var$pi_uf,2),  b=round(delta_var_pi_f_pi_uf,2))), cex=.8, pos=4)
text(start.x, start -4*shift, "state dependence", pos=4, cex=.8)
text(start.x, start -5*shift, substitute(paste(r[sigma], "(LGM* / PI*) =", a, "±", b), list(a=round(var$lgm_f/var$pi_f,2), b=round(delta_var_lgm_f_pi_f,2))), cex=.8, pos=4)
text(start.x, start -6*shift, substitute(paste(r[sigma], "(LGM / PI) =", a,  "±", b), list(a=round(var$lgm_uf/var$pi_uf,2), b=round(delta_var_lgm_uf_pi_uf,2))), cex=.8, pos=4)

legend(-1.35, 3.75, legend=c("LGM*", "LGM", "PI*", "PI"), col=c(rep(COLS[["LGM"]],2),rep(COLS[["PI"]],2)), lty=c(1,2, 1,2), lwd=rep(2,4), cex=.9)
plot_axis(1, label="GMST anomaly [K]")
plot_axis(2,  label="Density")
box()

dev.off()

