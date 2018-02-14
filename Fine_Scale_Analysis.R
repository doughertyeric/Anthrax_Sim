library(lubridate)
library(data.table)
library(dplyr)
library(tlocoh)
library(adehabitatHR)
library(rgeos)
library(rgdal)
library(raster)
library(maptools)

#Determine extent of the site, add buffer, and create raster to which others snap
site <- data.frame(matrix(c(-25000, 25000, 25000, -25000, 25000, 25000, -25000, -25000),4,2))
colnames(site) <- c("Longitude", "Latitude")
for (i in 1:nrow(site)) {
  site$Longitude[i] <- site$Longitude[i] + 292000
  site$Latitude[i] <- site$Latitude[i] + 8025000  
} 
site <- SpatialPoints(site[ , c("Longitude","Latitude")], proj4string=CRS("+proj=utm +south +zone=35 +ellps=WGS84"))

poly <- Polygon(site)
ps <- Polygons(list(poly),1)
site.poly <- SpatialPolygons(list(ps))

extent.all <- extent(SpatialPoints(site, CRS("+proj=utm +south +zone=35 +ellps=WGS84")))
snap <- raster(resolution = c(10,10), 
               xmn = (extent.all@xmin - 2500), xmx = (extent.all@xmax + 2500), 
               ymn = (extent.all@ymin - 2500), ymx = (extent.all@ymax + 2500), 
               crs = "+proj=utm +south +zone=35 +ellps=WGS84")

#Read in the raster of the LIZ locations and sizes
LIZs <- read.csv('/Users/ericdougherty/Box Sync/Dissertation/Anthrax_Sim/Ecol_Letters_ABM/Landscape/Run12/LIZs_Buffers_Run12.csv')
LIZs <- LIZs[,-1]

LIZ.layer <- function(LIZs) {
  
  # Increase attractivity of hex with LIZ
  coords <- coordinates(LIZs[,2:3])
  coords <- SpatialPoints(coords, proj4string=CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
  
  #buffers <- c()
  #for (i in 1:nrow(LIZs)) {
  #  age <- LIZs[i,4]
  #  size.factor <- log(LIZs[i,6])/log(350)
  #  buffers[i] <- size.factor * (rnorm(1, (3 * (4 - age)), ((3 * (4 - age))/8)))
  #}
  
  buffers <- LIZs$buffers
  buffers.sp <- gBuffer(coords, buffers, byid=TRUE, id = seq(1,length(buffers),1))
  df <- data.frame(matrix(seq(1,length(buffers),1), ncol=1))
  df[,2] <- rep(1,length(buffers.sp@polygons))
  colnames(df) <- c("ID", "value")
  buffers.sp <- SpatialPolygonsDataFrame(buffers.sp, data = df)
  
  site <- data.frame(matrix(c(-25000, 25000, 25000, -25000, 25000, 25000, -25000, -25000),4,2))
  colnames(site) <- c("Longitude", "Latitude")
  for (i in 1:nrow(site)) {
    site$Longitude[i] <- site$Longitude[i] + 292000
    site$Latitude[i] <- site$Latitude[i] + 8025000  
  } 
  site <- SpatialPoints(site[ , c("Longitude","Latitude")], proj4string=CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
  
  poly <- Polygon(site)
  ps <- Polygons(list(poly),1)
  site.poly <- SpatialPolygons(list(ps))
  
  extent.all <- extent(SpatialPoints(site, CRS("+proj=utm +south +zone=35 +ellps=WGS84")))
  snap <- raster(resolution = c(10,10), 
                 xmn = (extent.all@xmin - 2500), xmx = (extent.all@xmax + 2500), 
                 ymn = (extent.all@ymin - 2500), ymx = (extent.all@ymax + 2500), 
                 crs = "+proj=utm +south +zone=35 +ellps=WGS84")
  
  #strt <- Sys.time()
  LIZ.layer <- rasterize(buffers.sp, snap, background=0, getCover=TRUE)
  #print(Sys.time() - strt)
  
  return(LIZ.layer)
}

LIZ.layer <- LIZ.layer(LIZs)

################################################################################

#Loop over individuals
all.probs <- data.frame(matrix(0,20,2))
for (k in 11:20) {
  
  extent.all <- extent(SpatialPoints(site, CRS("+proj=utm +south +zone=35 +ellps=WGS84")))
  snap <- raster(resolution = c(10,10), 
                 xmn = (extent.all@xmin - 2500), xmx = (extent.all@xmax + 2500), 
                 ymn = (extent.all@ymin - 2500), ymx = (extent.all@ymax + 2500), 
                 crs = "+proj=utm +south +zone=35 +ellps=WGS84")
  
  #int.pts <- read.csv(paste0("/Users/ericdougherty/Box Sync/Dissertation/Anthrax_Sim/Ecol_Letters_ABM/Paths/Run12/ID",k,"_IntermediatePts_Run12.csv"))
  #int.pts <- read.csv(paste0("/global/home/users/edougher/ID",k,"_IntermediatePts_Run12.csv"))
  int.pts <- read.csv(paste0('/Users/ericdougherty/Desktop/Ecol_Letters_ABM/Final_Outputs/Run12/ID',k,'_IntermediatePts_Run12.csv'))
  int.sp.pts <- SpatialPoints(int.pts[,2:3], CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
  dens <- rasterize(int.sp.pts, snap, fun='count', background=0)
  extent(dens) <- snap
  dens.non.zero <- na.omit(sum(dens[dens > 0]))
  dens <- dens / dens.non.zero
  
  extent.all <- extent(SpatialPoints(site, CRS("+proj=utm +south +zone=35 +ellps=WGS84")))
  snap <- raster(resolution = c(10,10), 
                 xmn = (extent.all@xmin - 2500), xmx = (extent.all@xmax + 2500), 
                 ymn = (extent.all@ymin - 2500), ymx = (extent.all@ymax + 2500), 
                 crs = "+proj=utm +south +zone=35 +ellps=WGS84")
  
  reduced <- int.pts[int.pts[,4] == 2 | int.pts[,5] == 2,]
  reduced.sp.pts <- SpatialPoints(reduced[,2:3], CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
  reduced.dens <- rasterize(reduced.sp.pts, snap, fun='count', background=0)
  extent(reduced.dens) <- snap
  #reduced.non.zero <- na.omit(sum(reduced.dens[reduced.dens > 0]))
  reduced.dens <- reduced.dens / dens.non.zero
  
  #########
  
  dens.prob <- LIZ.layer * dens
  reduced.dens.prob <- LIZ.layer * reduced.dens
  
  all.probs[k,1] <- cellStats(dens.prob, stat='sum')
  all.probs[k,2] <- cellStats(reduced.dens.prob, stat='sum')
  print(k)
}

write.csv(all.probs, "Fine_Scale_Probabilities.csv")
