#THIS SCRIPT IS CALLED BY FS9.R

#---DATA---#
save <- F

runs <- c("xmzke","xmzkg","xmzkh","xmzka","xmzkd","xmzki","xnagd","xnagf","xnagg","xnagb","xnage","xnagh")
path <-  "data/HadCM3/"
savefile <- paste0("data/max-amoc.RData")
  
nc <- nc_open("data/amoc_grid.nc")
lats <- ncvar_get(nc, "latitude")
lons <- ncvar_get(nc, "longitude")
nc_close(nc)
gc()

load("data/categorized_eruptions.RData")
eruptions <- NULL
D <- NULL

zonal_width <- 6.3781e6*tan(1.25/360*2*pi)*cos(lats[1:(length(lats)-1)]/180*pi)

amoc.max <- list()
daysSince <- list()

if(!file.exists(savefile)){
  for (run in runs)
    {
    nc <- nc_open(paste(path,run,".nc",sep="")) 
    vel <- ncvar_get(nc, "field704")
    daysSince[[run]] <- ncvar_get(nc, "t")
    d <- dim(vel)
    zonal_width <- array(zonal_width, d[1:3])
    depth <- ncvar_get(nc, "depth")
    height <- t(array(depth-c(0,depth[1:(length(depth)-1)]),d[2:1]))
    height <- array(height,c(dim(height),d[3]))
    nc_close(nc)
    
    flow <- 0.01*vel*zonal_width*height/1e6
    inf.mask <- which(abs(flow) == Inf,arr.ind=T)
    flow[inf.mask] <- NA
    flow <- flow[,d[2]:1,]
    amoc <- array(NA,d[1:3])
    for (z in 1:d[2])
    {
      amoc[,z,] <- -apply(array(flow[,1:z,], c(d[1],z,d[3])), c(1,3), sum, na.rm=T)
    }
    # Assert that at every timestep there is at least one not-NA value
    stopifnot(sum(apply(!is.na(amoc[33:67,21-(11:15),]),3,sum) == 0) == 0)
    
    amoc.max[[run]] <- apply(amoc[33:67,21-(11:15),],c(3),max,na.rm=T)
  }
  save("amoc.max", "tAD", "daysSince",file=savefile)
} else {
  load(savefile)
}

timesteps <- sort(as.integer(unique(unlist(daysSince))))
n_timesteps <- length(timesteps)
timesteps.range <- range(timesteps)
export_amoc <- matrix(NA,nrow=n_timesteps,ncol=length(runs)+1)
rownames(export_amoc) <- DaysSinceToAD(timesteps)
colnames(export_amoc) <- c("daysSince",runs)

export_amoc[,1] <- sort(as.integer(unique(unlist(daysSince))))

for (i in 2:(length(runs)+1))
{
  ind <- which(timesteps %in% daysSince[[i-1]])
  export_amoc[ind,i] <- amoc.max[[i-1]]
}

if(save) write.xlsx(export_amoc, "data/amoc_hadcm3.xlsx", colNames=T, rowNames=T, sheetName="amoc",append=FALSE)
