library(sp)
library(rgeos)
library(raster)
library(foreach)
library(doParallel)
library(truncnorm)

determine.vmax <- function(xy, window) {
  vmax <- 0
  for (i in 1:(nrow(xy) - (window))) {
    points <- data.frame(xy[i:(i+window),])
    total.dist <- 0
    for (j in 1:(nrow(points)-1)) {
      dist.x <- as.numeric(points[j+1,1] - points[j,1])
      dist.y <- as.numeric(points[j+1,2] - points[j,2])
      dist <- sqrt(dist.x^2 + dist.y^2)
      total.dist <- total.dist + dist
    }
  if (total.dist > vmax) {
    vmax <- total.dist
  }
  }
  return(vmax)
}

determine.maxstep <- function(path.length, gamma.mean, gamma.sd) {
  temp.steps <- rgamma(5*path.length, shape=((gamma.mean^2)/(gamma.sd^2)), rate=(gamma.mean/(gamma.sd^2)))
  temp.steps <- temp.steps[order(temp.steps, decreasing = TRUE)]
  max.step <- temp.steps[1]
  #max.step <- mean(top.steps)
  return(max.step)
}

sim.paths <- function(pt1, pt2, vmax, time.steps, n.sim) {
  sub.v <- vmax/time.steps
  pt2 <- SpatialPoints(pt2, CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
  intermediate.pts <- data.frame(matrix(0,(time.steps*n.sim),2))
  for (j in 1:n.sim) {
    new.pos <- SpatialPoints(pt1, CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
    intermediate.pts[(((j-1)*time.steps)+1),] <- pt1
    for (n in 2:(time.steps-1)) {
      origin.circle.rad <- sub.v*n
      origin.circle <- gBuffer(spgeom=new.pos, width=origin.circle.rad)
      destination.circle.rad <- sub.v*(time.steps-n)
      destination.circle <- gBuffer(spgeom=pt2, width=destination.circle.rad)
      new.pos.region <- gIntersection(origin.circle, destination.circle)
      new.pos <- spsample(new.pos.region, n=1, type='random', iter=10)
      intermediate.pts[(((j-1)*time.steps)+n),] <- new.pos@coords
    }
    intermediate.pts[(((j-1)*time.steps)+time.steps),] <- pt2@coords
  }
  return(intermediate.pts)
}

sim.path.density <- function(pt1, pt2, vmax, time.steps, n.sim) {
  sub.v <- vmax/time.steps
  pt2 <- SpatialPoints(pt2, CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
  intermediate.pts <- data.frame(matrix(0,(time.steps*n.sim),2))
  for (j in 1:n.sim) {
    new.pos <- SpatialPoints(pt1, CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
    intermediate.pts[(((j-1)*time.steps)+1),] <- pt1
    for (n in 2:(time.steps-1)) {
      origin.circle.rad <- sub.v*n
      origin.circle <- gBuffer(spgeom=new.pos, width=origin.circle.rad)
      destination.circle.rad <- sub.v*(time.steps-n)
      destination.circle <- gBuffer(spgeom=pt2, width=destination.circle.rad)
      new.pos.region <- gIntersection(origin.circle, destination.circle)
      new.pos <- spsample(new.pos.region, n=1, type='random', iter=10)
      intermediate.pts[(((j-1)*time.steps)+n),] <- new.pos@coords
    }
    intermediate.pts[(((j-1)*time.steps)+time.steps),] <- pt2@coords
  }
  intermediate.sp.pts <- SpatialPoints(intermediate.pts, CRS("+proj=utm +south +zone=35 +ellps=WGS84"))

  max.extent <- extent(intermediate.sp.pts)
  r <- raster(max.extent)
  projection(r) <- CRS("+proj=utm +south +zone=35 +ellps=WGS84")
  res(r) <- c(10,10)
  dens <- rasterize(intermediate.sp.pts, r, fun='count', background=0)
  dens@data@values <- dens@data@values/n.sim
  dens@data@values[dens@data@values > 0.99] <- 0.5
  out.list <- list(intermediate.pts, dens)
  return(out.list)
}

mosaicList <- function(raster.list, snap){
  
  for (i in 1:length(raster.list)) {
    raster.list[[i]] <- projectRaster(raster.list[[i]], snap, method = "ngb",
                                      crs = "+proj=utm +south +zone=35 +ellps=WGS84")
  }
  
  # edit settings of the raster list for use in do.call and mosaic
  names(raster.list) <- NULL
  #####This function deals with overlapping areas
  raster.list$fun <- sum
  
  #run do call to implement mosaic over the list of raster objects.
  mos <- do.call(raster::mosaic, raster.list)
  
  #set crs of output
  crs(mos) <- crs(x = "+proj=utm +south +zone=35 +ellps=WGS84")
  return(mos)
}

#########

test.zeb <- read.csv("~/Box Sync/Dissertation/Anthrax_Sim/Ecol_Letters_ABM/Paths/Run13/ID1_Steps_Run13.csv")

xy <- test.zeb[,2:3]
state.vector <- test.zeb[,6]
n.sim <- 1
max.step.state1 <- determine.maxstep(nrow(xy), 8.619, 7.840)
max.step.state2 <- determine.maxstep(nrow(xy), 106, 13)
max.step.state3 <- determine.maxstep(nrow(xy), 628, 403)

extent.all <- extent(SpatialPoints(xy, CRS("+proj=utm +south +zone=35 +ellps=WGS84")))
snap <- raster(resolution = c(10,10), 
               xmn = (extent.all@xmin - max.step.state3), xmx = (extent.all@xmax + max.step.state3), 
               ymn = (extent.all@ymin - max.step.state3), ymx = (extent.all@ymax + max.step.state3), 
               crs = "+proj=utm +south +zone=35 +ellps=WGS84")

#dens@data@values <- log(dens@data@values)
#plot(dens)
num_cores <- detectCores()
cl<-makeCluster(num_cores)
registerDoParallel(cl)

#### Non-parallelized ####
#all.intermediate.pts <- list()
#density.rasters <- list()
#for (i in 1:(nrow(xy) - 1)) {
#  state.val <- state.vector[i]
#  if (state.val == 1) {
#    path.dens <- sim.path.density(xy[i,], xy[(i+1),], max.step.state1, 21, n.sim)
#  } else if (state.val == 2) {
#    path.dens <- sim.path.density(xy[i,], xy[(i+1),], max.step.state2, 21, n.sim)
#  } else if (state.val == 3) {
#    path.dens <- sim.path.density(xy[i,], xy[(i+1),], max.step.state3, 21, n.sim)
#  }
#  all.intermediate.pts[[i]] <- path.dens[[1]]
#  density.rasters[[i]] <- path.dens[[2]]
#}

#### Parallelized Point Output ####
strt <- Sys.time()
#simulated.paths <- foreach(z = 1:20, .combine='rbind', .packages=c("sp", "rgeos", "raster"), .export=c("sim.paths")) %dopar% {
path.dens <- data.frame(matrix(0,0,2))
path.dens.temp <- data.frame(matrix(0,0,2))
for (i in 1:(nrow(xy) - 1)) {
  state.val1 <- as.numeric(state.vector[i])
  state.val2 <- as.numeric(state.vector[i+1])
  dist.x <- abs(xy[i,1] - xy[(i+1),1])
  dist.y <- abs(xy[i,2] - xy[(i+1),2])
  current.dist <- sqrt(dist.x^2 + dist.y^2)
  
  if (state.val2 == 2) {
    if (current.dist < max.step.state2) {
      #path.dens.temp <- sim.paths(xy[i,], xy[(i+1),], max.step.state2, 21, n.sim)
      path.dens.temp <- sim.paths(xy[i,], xy[(i+1),], (rtruncnorm(n=1,mean=1.2,sd=0.05,a=1.05)*current.dist), 21, n.sim)
    } else {
      path.dens.temp <- sim.paths(xy[i,], xy[(i+1),], (1.2*current.dist), 21, n.sim)
    }
  } else if (state.val2 == 3) {
    if (current.dist < max.step.state3) {
      #path.dens.temp <- sim.paths(xy[i,], xy[(i+1),], max.step.state3, 21, n.sim)
      path.dens.temp <- sim.paths(xy[i,], xy[(i+1),], (rtruncnorm(n=1,mean=1.2,sd=0.05,a=1.05)*current.dist), 21, n.sim)
    } else {
      path.dens.temp <- sim.paths(xy[i,], xy[(i+1),], (1.2*current.dist), 21, n.sim)
    }
  } else if (state.val2 == 1) {
    if (current.dist < max.step.state1) {
      #path.dens.temp <- sim.paths(xy[i,], xy[(i+1),], max.step.state1, 21, n.sim)
      path.dens.temp <- sim.paths(xy[i,], xy[(i+1),], (rtruncnorm(n=1,mean=1.2,sd=0.05,a=1.05)*current.dist), 21, n.sim)
    } else {
      path.dens.temp <- sim.paths(xy[i,], xy[(i+1),], (1.2*current.dist), 21, n.sim)
    }
  }
  path.dens <- data.frame(rbind(path.dens, path.dens.temp))
  print(i)
}
print(Sys.time() - strt)

#### Parallelized Point and Raster Output ####
#strt <- Sys.time()
#simulated.paths <- foreach(i = 1:9, .packages=c("sp", "rgeos", "raster"), .export=c("sim.path.density")) %dopar% {
#  state.val <- state.vector[i]
#  if (state.val == 1) {
#    path.dens <- sim.path.density(xy[i,], xy[(i+1),], max.step.state1, 21, n.sim)
#  } else if (state.val == 2) {
#    path.dens <- sim.path.density(xy[i,], xy[(i+1),], max.step.state2, 21, n.sim)
#  } else if (state.val == 3) {
#    path.dens <- sim.path.density(xy[i,], xy[(i+1),], max.step.state3, 21, n.sim)
#  }
#}
#print(Sys.time() - strt)
