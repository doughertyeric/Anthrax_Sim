library(lubridate)
library(data.table)
library(dplyr)
library(tlocoh)
library(adehabitatHR)
library(rgeos)
library(rgdal)
library(raster)
library(sf)
library(maptools)


##########################################################################################################

train.test <- function (data, seed) {
  df <- data.frame(matrix(TRUE,nrow(data),100))
  set.seed(seed)
  
  for (i in 1:ncol(df)) {
    samp <- sample(seq(1,nrow(data),1), round(0.002222*(nrow(data))))
    for (j in 1:length(samp)) {
      for (k in 1:nrow(df)) {
        if (k == samp[j]) {
          df[k,i] <- FALSE
        }
      }
    }
  }
  
  count <- 0
  for (i in 1:ncol(df)) {
    new <- table(df[,i])[1]
    count <- count + new
  }
  
  for (i in 3:nrow(df)) {
    for (j in 1:ncol(df)) {
      if (df[i-1,j] == FALSE && df[i-2,j] != FALSE) {
        for (k in 0:98) {
          df[i+k,j] <- FALSE
        }
      }
    }
  }
  
  df <- df[1:nrow(data),]
  return(df)
}

algo <- function(df, k.vals, data.lxy, data) {
  
  count <- 0
  for (i in 1:ncol(df)) {
    new <- table(df[,i])[1]
    count <- count + new
  }
  
  trace <- data.frame(matrix(0,length(k.vals),6))
  
  for (k in 1:length(k.vals)) {
    
    current.k_val <- k.vals[k]
    
    # Calculate the nearest neighbors and create lhs object for full dataset
    full.lxy <- lxy.nn.add(data.lxy, k=current.k_val, s=0, status=F)
    full.lhs <- lxy.lhs(full.lxy, k=current.k_val, s=0, status=F)
    coords <- full.lhs[[1]]$pts@coords
    
    # Create list for the negative log likelihood values for each test/train split in df
    likelihood <- list()
    
    for (j in 1:ncol(df)) {
      
      # Create a one-column data frame from df
      df1 <- df[1:nrow(data),j]
      
      # Create selection of hulls based on Boolean
      hulls.sel.idx <- which(df1)
      full.hulls <- hulls(full.lhs)
      selected.hulls <- full.hulls[[1]] [ full.hulls[[1]]@data$pts.idx %in% hulls.sel.idx , ]
      
      # Determine the number of points in training dataset
      df.temp <- data.frame(as.numeric(df[,j]))
      colnames(df.temp) <- "subset"
      total.pts <- sum(df.temp)
      
      # Find middle points of testing data and define as -1
      for (i in 2:nrow(df.temp)) {
        if (df.temp[i,1] == 0 && df.temp[i-1,1] != 0 && df.temp[i-1,1] != -1) {
          df.temp[i+50,1] <- -1
        }
      }
      
      df.temp <- df.temp[1:nrow(data),]
      df.temp <- cbind(coords, df.temp)
      
      # Extract the middle points for testing and record associated coordinates
      test.pts <- data.frame()
      q = 1
      for (i in 1:nrow(df.temp)) {
        if (df.temp[i,3] == -1) {
          test.pts[q,1] <- df.temp[i,1]
          test.pts[q,2] <- df.temp[i,2]
          q = q + 1
        }
      }
      
      colnames(test.pts) <- c("x", "y")
      test.pts <- SpatialPoints(test.pts[ , c("x","y")], proj4string=CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
      
      poly <- SpatialPolygons(selected.hulls@polygons, proj4string = CRS('+proj=utm +south +zone=35 +ellps=WGS84'))
      
      # Determine the number of hulls under each test point,
      overlay <- data.frame(matrix(0,length(test.pts@coords[,1]),1))
      
      for (i in 1:length(test.pts@coords[,1])) {
        overlay.list <- over(test.pts[i], poly, returnList=TRUE)
        overlay[i,1] <- length(overlay.list[[1]])
      }
      
      hull.mean <- mean(full.lhs[[1]]$hulls@data$area)
      
      # calculate likelihood
      for (i in 1:nrow(overlay)) {
        overlay[i,2] <- overlay[i,1]/length(selected.hulls[[1]])
        overlay[i,3] <- -log(overlay[i,2])
        overlay[i,4] <- log(overlay[i,2],exp(1))
        if (overlay[i,1] == 0) {
          overlay[i,3] <- -log((1/nrow(data))/100)
          overlay[i,4] <- log((1/nrow(data))/100, exp(1))
          overlay[i,4] <- log((1/nrow(data))/100, exp(1))
        }
      }
      
      colnames(overlay) <- c("over", "prob", "loglike", "ln.like")
      
      # Add values to likelihood list
      likelihood[[j]] <- as.list(overlay)
    }
    
    log.like <- data.frame(matrix(0,length(likelihood),4))
    for (i in 1:length(likelihood)) {
      log.like[i,1] <- sum(likelihood[[i]]$loglike, na.rm=TRUE)
      log.like[i,2] <- sum(likelihood[[i]]$ln.like, na.rm=TRUE)
      log.like[i,3] <- -2*(log.like[i,2]) + 2*k
      log.like[i,4] <- -2*(log.like[i,2]) + k*log(as.numeric(count),exp(1))
    }
    #Sys.time()
    
    # Assign mean of the means of negative log likelihoods as the current posterior probability for Metropolis-Hastings
    new.postLike <- sum(log.like[,1])
    new.lnLike <- sum(log.like[,2])
    new.AIC <- sum(log.like[,3])
    new.BIC <- sum(log.like[,4])
    
    #cat("iteration:", counter, "chain:", current.k_val, current.s_val, "likelihood:", new.postLike, "\n")
    
    trace[k,1] <- current.k_val
    trace[k,2] <- hull.mean
    trace[k,3] <- new.postLike
    trace[k,4] <- new.lnLike
    trace[k,5] <- new.AIC
    trace[k,6] <- new.BIC
    
  } # End of k loop
  
  return(trace) 
}

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

#########################################################################################################################

#Import simulated paths and create a list object
temp = list.files(pattern="*_Steps_Run13.csv")
ids <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
movementpaths <- lapply(temp, read.csv)
movementpaths <- lapply(movementpaths, setNames, c("X.1", "stepID", "X", "Y", "stepsize", 
                                                   "angle", "mode", "ID", "Time"))
movementpaths <- lapply(movementpaths, mutate, Time = ymd_hms(Time))

temp = list.files(pattern="*_kmeans_Run13.csv")
movement.paths2 <- lapply(temp, read.csv)

all.paths <- list()
for (i in 1:length(movementpaths)) {
  move1 <- movementpaths[[i]]
  move2 <- movement.paths2[[i]]
  move.temp <- data.frame(cbind(move1[,2:9], move2[,10]))
  colnames(move.temp) <- c("stepID", "X", "Y", "stepsize", "angle", "mode", "ID", "Time", "kmeans")
  all.paths[[i]] <- move.temp
}

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
#LIZs <- raster('/Users/ericdougherty/Box Sync/Dissertation/Anthrax_Sim/Ecol_Letters_ABM/Landscape/Run12/LIZs_Run12.tif')
LIZs <- read.csv('LIZs_Buffers_Run13.csv')
LIZs <- LIZs[,-1]
LIZ.layer <- LIZ.layer(LIZs)

##########################################################################################################

#Loop over the individuals to determine the MCP, locoh home range, then output the probability of contact
all.probs.broad <- data.frame(matrix(0,20,2))
for (i in 1:length(all.paths)) {
  
  #Remove each individual and project points
  temp <- all.paths[[i]]
  temp.sp <- SpatialPoints(temp[ , c("X","Y")], CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
  temp.pts <- coordinates(temp.sp)
  colnames(temp.pts) <- c("x","y")
  
  #Create 95% MCP from the SpatialPoints object, assign probability values to each cell
  temp.mcp <- mcp(temp.sp, percent=95)
  MCP <- rasterize(temp.mcp, snap, background=0, fun='first')#getCover=TRUE)
  MCP.non.zero <- na.omit(length(MCP[MCP > 0]))
  MCP <- MCP * (0.95 / MCP.non.zero)
  
  for (j in 1:(nrow(temp) - 1)) {
    state.val1 <- as.numeric(temp$kmeans[j])
    state.val2 <- as.numeric(temp$kmeans[j+1])
    temp[j,10] <- state.val1
    temp[j,11] <- state.val2
  }
  reduced <- temp[temp[,10] == 2 | temp[,11] == 2,]
  reduced <- reduced[!is.na(reduced$X),]
  reduced.sp <- SpatialPoints(reduced[ , c("X","Y")], CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
  reduced.mcp <- mcp(reduced.sp, percent=95)
  reduced.MCP <- rasterize(reduced.mcp, snap, background=0, fun='first')#getCover=TRUE)
  reduced.MCP <- reduced.MCP * (0.95 / MCP.non.zero)
  
  #Use the LIZ layer to determine the overall probability of contact (and the number of encounter cells)
  MCP.prob <- LIZ.layer * MCP
  reduced.MCP.prob <- LIZ.layer * reduced.MCP

  all.probs.broad[i,1] <- cellStats(MCP.prob, stat='sum')
  all.probs.broad[i,2] <- cellStats(reduced.MCP.prob, stat='sum')

  print(i)
  
}

write.csv(all.probs.broad, "Broad_Scale_Probabilities_MCP.csv")

#############################################################################################

opt.k <- data.frame(matrix(0,20,2))
for (i in 1:length(all.paths)) {
  #Remove each individual and project points
  temp <- all.paths[[i]]
  temp.sp <- SpatialPoints(temp[ , c("X","Y")], CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
  temp.pts <- coordinates(temp.sp)
  colnames(temp.pts) <- c("x","y")
  
  locoh.temp <- temp[seq(1, nrow(temp), 3),]
  df <- train.test(locoh.temp, seed=1234)
  locoh.temp.sp <- SpatialPoints(locoh.temp[ , c("X","Y")], CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
  locoh.temp.pts <- coordinates(locoh.temp.sp)
  colnames(locoh.temp.pts) <- c("x","y")
  
  #Create LoCoH home range, assign probability values to each cell
  locoh.temp.lxy <- xyt.lxy(xy=locoh.temp.pts, dt=locoh.temp$Time, id=locoh.temp$ID[1], proj4string=CRS("+proj=utm +south +zone=35 +ellps=WGS84"), dup.dt.check=FALSE)
  strt <- Sys.time()
  trace <- algo(df, k.vals=seq(2,30,2), locoh.temp.lxy, locoh.temp)
  print(Sys.time() - strt)
  min.val <- which.min(trace[,5]) #AIC in column 5, BIC in column 6
  opt.k[i,1] <- trace[min.val,1]
  
  for (j in 1:(nrow(temp) - 1)) {
    state.val1 <- as.numeric(temp$kmeans[j])
    state.val2 <- as.numeric(temp$kmeans[j+1])
    temp[j,10] <- state.val1
    temp[j,11] <- state.val2
  }
  reduced <- temp[temp[,10] == 2 | temp[,11] == 2,]
  reduced <- reduced[!is.na(reduced$X),]
  reduced.temp <- reduced[seq(1, nrow(reduced), 3),]
  df2 <- train.test(reduced.temp, seed=1234)
  reduced.sp <- SpatialPoints(reduced.temp[ , c("X","Y")], CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
  reduced.temp.pts <- coordinates(reduced.sp)
  colnames(reduced.temp.pts) <- c("x","y")

  reduced.lxy <- xyt.lxy(xy=reduced.temp.pts, dt=reduced.temp$Time, id=reduced.temp$ID[1], proj4string=CRS("+proj=utm +south +zone=35 +ellps=WGS84"), dup.dt.check=FALSE)
  strt <- Sys.time()
  trace2 <- algo(df2, k.vals=seq(2,30,2), reduced.lxy, reduced.temp)
  print(Sys.time() - strt)
  min.val2 <- which.min(trace2[,5])
  opt.k[i,2] <- trace2[min.val2,1]
  
  print(i)
}

#write.csv(opt.k, "LoCoH_Optimal_k.csv")
opt.k <- read.csv("LoCoH_Optimal_k.csv")
opt.k <- opt.k[,-1]

all.probs.broad2 <- data.frame(matrix(0,20,2))
for (i in c(1,2,4,5,7,8,9,10,11,12,13,14,17,18)) {
#for (i in 1:length(all.paths)) {
  temp <- all.paths[[i]]
  temp.sp <- SpatialPoints(temp[ , c("X","Y")], CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
  temp.pts <- coordinates(temp.sp)
  colnames(temp.pts) <- c("x","y")
  
  locoh.temp <- temp[seq(1, nrow(temp), 3),]
  #df <- train.test(locoh.temp, seed=1234)
  locoh.temp.sp <- SpatialPoints(locoh.temp[ , c("X","Y")], CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
  locoh.temp.pts <- coordinates(locoh.temp.sp)
  colnames(locoh.temp.pts) <- c("x","y")
  
  #Create LoCoH home range, assign probability values to each cell
  locoh.temp.lxy <- xyt.lxy(xy=locoh.temp.pts, dt=locoh.temp$Time, id=locoh.temp$ID[1], proj4string=CRS("+proj=utm +south +zone=35 +ellps=WGS84"), dup.dt.check=FALSE)
  temp.lxy <- lxy.nn.add(locoh.temp.lxy, k=opt.k[i,1], s=0, status=F)
  temp.lhs <- lxy.lhs(temp.lxy, k=opt.k[i,1], s=0, status=F)
  temp.lhs <- lhs.iso.add(temp.lhs)
  
  lhs.exp.shp(temp.lhs, iso=TRUE, dir=".", file.base = "temp", status=F, show.time=F)
  temp.shp <- readShapePoly(paste0("temp.iso.srt-area.iso-q.h",as.character(nrow(locoh.temp)),".i5.00.iso.shp"))
  temp.shp <- SpatialPolygons(temp.shp@polygons, proj4string = CRS('+proj=utm +south +zone=35 +ellps=WGS84'))
  locoh <- rasterize(temp.shp, snap, background=0, fun='first')#, getCover=TRUE)
  locoh.non.zero <- na.omit(length(locoh[locoh > 0]))
  locoh <- locoh * (0.95 / locoh.non.zero)
  
  locoh.prob <- LIZ.layer * locoh
  all.probs.broad2[i,1] <- cellStats(locoh.prob, stat='sum')
  print(all.probs.broad2[i,1])
  
  #Remove temporary .shp files from locoh construction process
  file.remove(paste("temp.iso.srt-area.iso-q.h",as.character(nrow(locoh.temp)),".i5.00.iso.dbf", sep = ""))
  file.remove(paste("temp.iso.srt-area.iso-q.h",as.character(nrow(locoh.temp)),".i5.00.iso.prj", sep = ""))
  file.remove(paste("temp.iso.srt-area.iso-q.h",as.character(nrow(locoh.temp)),".i5.00.iso.shp", sep = ""))
  file.remove(paste("temp.iso.srt-area.iso-q.h",as.character(nrow(locoh.temp)),".i5.00.iso.shx", sep = ""))
  
  for (j in 1:(nrow(temp) - 1)) {
    state.val1 <- as.numeric(temp$kmeans[j])
    state.val2 <- as.numeric(temp$kmeans[j+1])
    temp[j,10] <- state.val1
    temp[j,11] <- state.val2
  }
  reduced <- temp[temp[,10] == 2 | temp[,11] == 2,]
  reduced <- reduced[!is.na(reduced$X),]
  reduced.temp <- reduced[seq(1, nrow(reduced), 3),]
  reduced.sp <- SpatialPoints(reduced.temp[ , c("X","Y")], CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
  reduced.temp.pts <- coordinates(reduced.sp)
  colnames(reduced.temp.pts) <- c("x","y")
  
  reduced.lxy <- xyt.lxy(xy=reduced.temp.pts, dt=reduced.temp$Time, id=reduced.temp$ID[1], proj4string=CRS("+proj=utm +south +zone=35 +ellps=WGS84"), dup.dt.check=FALSE)
  reduced.lxy <- lxy.nn.add(reduced.lxy, k=opt.k[i,2], s=0, status=F)
  reduced.lhs <- lxy.lhs(reduced.lxy, k=opt.k[i,2], s=0, status=F)
  reduced.lhs <- lhs.iso.add(reduced.lhs)
  
  lhs.exp.shp(reduced.lhs, iso=TRUE, dir=".", file.base = "reduced", status=F, show.time=F)
  reduced.shp <- readShapePoly(paste0("reduced.iso.srt-area.iso-q.h",as.character(nrow(reduced.temp.pts)),".i5.00.iso.shp"))
  reduced.shp <- SpatialPolygons(reduced.shp@polygons, proj4string = CRS('+proj=utm +south +zone=35 +ellps=WGS84'))
  reduced <- rasterize(reduced.shp, snap, background=0, fun='first')#, getCover=TRUE)
  reduced <- reduced * (0.95 / locoh.non.zero)
  
  #Use the LIZ layer to determine the overall probability of contact (and the number of encounter cells)
  reduced.prob <- LIZ.layer * reduced
  all.probs.broad2[i,2] <- cellStats(reduced.prob, stat='sum')
  print(all.probs.broad2[i,2])
  
  #Remove temporary .shp files from locoh construction process
  file.remove(paste("reduced.iso.srt-area.iso-q.h",as.character(nrow(reduced.temp.pts)),".i5.00.iso.dbf", sep = ""))
  file.remove(paste("reduced.iso.srt-area.iso-q.h",as.character(nrow(reduced.temp.pts)),".i5.00.iso.prj", sep = ""))
  file.remove(paste("reduced.iso.srt-area.iso-q.h",as.character(nrow(reduced.temp.pts)),".i5.00.iso.shp", sep = ""))
  file.remove(paste("reduced.iso.srt-area.iso-q.h",as.character(nrow(reduced.temp.pts)),".i5.00.iso.shx", sep = ""))
  
  print(i)
}

write.csv(all.probs.broad2, "Broad_Scale_Probabilities_LoCoH.csv")
