library(sp)
library(rgeos)
library(prevR)
library(foreach)
library(doParallel)
library(fdrtool)
library(plyr)
library(stringr)
library(latticeExtra)
library(truncnorm)
library(CircStats)
library(raster)

###### Basic Functions ######

# Basic function for points within polygon
within.func <- function(pts,ply){
  within.temp <- gIntersects(ply, pts, byid=TRUE)
  within <- which(within.temp == TRUE)
  return(within)
}

outside.func <- function(pts,ply){
  outside.temp <- gIntersects(ply, pts, byid=TRUE)
  outside <- which(outside.temp == FALSE)
  return(outside)
}

# Calculates angle between two points
anglefun <- function(xx, yy, bearing=TRUE) {
  ## calculates the compass bearing of the line between two points
  ## xx and yy are the differences in x and y coordinates between two points
  ## Options:
  ## bearing = FALSE returns +/- pi instead of 0:2*pi
  
  b<-sign(yy)
  b[b==0]<-1  #corrects for the fact that sign(0) == 0
  tempangle = b*(xx<0)*pi+atan(yy/xx)
  if(bearing){
    #return a compass bearing 0 to 2pi
    #if bearing==FALSE then a heading (+/- pi) is returned
    tempangle[tempangle<0]<-tempangle[tempangle<0]+2*pi
  }
  return(tempangle)
}

# Assign positions at each next step for one individual
movement <- function(xy, step, heading, TimeSinceSunrise) {
  
  pi = 3.141593
  x_init <- xy[1,1]
  y_init <- xy[1,2]
  
  if (heading < 0) {
    heading <- abs(heading) + pi
  }
  
  #rad_y <- angle*0.0174532925
  y_change <- sin(heading)*step
  y_new <- y_init + y_change
  
  # Use cosine to determine the movement in x (longitude)
  #rad_x <- angle*0.0174532925
  x_change <- cos(heading)*step
  x_new <- x_init + x_change
  
  x_init <- x_new
  y_init <- y_new
  
  move.temp <- as.data.frame(matrix(0,1,3))
  move.temp[1,1] <- x_new
  move.temp[1,2] <- y_new
  move.temp[1,3] <- TimeSinceSunrise
  
  return(move.temp)
}

###### Initialization ######

# Create hetergeneous landscape based on forage quality
landscape.hetero <- function(site, cell.size, prop.graze, mean.qual) {
  
  for (i in 1:nrow(site)) {
    site$Longitude[i] <- site$Longitude[i] + 292000
    site$Latitude[i] <- site$Latitude[i] + 8025000  
  } 
  site <- SpatialPoints(site[ , c("Longitude","Latitude")], proj4string=CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
  
  poly <- Polygon(site)
  ps <- Polygons(list(poly),1)
  site.poly <- SpatialPolygons(list(ps))
  
  hex_points <- spsample(site.poly, type = "hexagonal", cellsize = cell.size, offset = c(0,0))
  hex_points@proj4string <- CRS("+proj=utm +south +zone=35 +ellps=WGS84")
  hex_grid <- HexPoints2SpatialPolygons(hex_points, dx = cell.size)
  
  hex_area <- hex_grid@polygons[[1]]@Polygons[[1]]@area
  mean.cell.qual <- (mean.qual / 1000000)*hex_area
  
  # Randomly assign values for quality of browse and forage in hex
  veg <- data.frame(matrix(0,length(hex_grid@polygons),6))
  for (i in 1:nrow(veg)){
    veg[i,1] <- i
    veg[i,2] <- rbinom(1,1,prop.graze) #whether this is a max cell or not
    veg[i,3] <- rtruncnorm(n=1, mean=mean.cell.qual, sd=(mean.cell.qual/4), a=0)
    veg[i,4] <- 0 #carrying capacity for graze based on initial value
    veg[i,5] <- 0 #density of LIZs
    veg[i,6] <- 1 
    veg[i,7] <- rbinom(1,1,prop.graze) #whether this is a min graze cell or not
    #veg[i,8] <- runif(1,1,7) # browse value
    #veg[i,9] <- 0 # NDVI (mean)
  }
  
  for (k in 1:nrow(veg)) {
    if (veg[k,2] == 1) {
      id <- veg[k,1]
      pos <- hex_points[id]
      quality.radius <- (1000/cell.size)*cell.size*abs(rnorm(1,1.5,0.5))
      nearby <- gBuffer(pos, width=quality.radius, byid=TRUE)
      within <- within.func(hex_points, nearby)
      for (j in 1:length(within)) {
        id2 <- within[j]
        veg[id2,3] <- rtruncnorm(n=1, mean=(4*(mean.cell.qual/3)), sd=(mean.cell.qual/4), a=0)
      }
    } else if (veg[k,7] == 1) {
      id <- veg[k,1]
      pos <- hex_points[id]
      quality.radius <- (1000/cell.size)*cell.size*abs(rnorm(1,1.5,0.5))
      nearby <- gBuffer(pos, width=quality.radius, byid=TRUE)
      within <- within.func(hex_points, nearby)
      for (j in 1:length(within)) {
        id2 <- within[j]
        veg[id2,3] <- rtruncnorm(n=1,mean=(2*(mean.cell.qual/3)),sd=(mean.cell.qual/4),a=0)
      }
    }
  }
  
  veg[,4] <- veg[,3] 
  #veg[,9] <- (veg[,3] + veg[,8])/2
  
  colnames(veg) <- c("hex_ID", "High_graze", "graze", "CC.graze", "LIZ.density", "region", "Low_graze")
  landscape <- SpatialPolygonsDataFrame(hex_grid, veg, match.ID=FALSE)
  landscape@proj4string <- CRS("+proj=utm +south +zone=35 +ellps=WGS84")
  
  landscape.list <- list(landscape, hex_points, hex_grid)
  return(landscape.list)
}

edge.adjustment <- function(site, landscape, cell.size, hex_points, hex_grid, mean.qual) {
  site.buff <- data.frame(matrix(0,4,2)) 
  colnames(site.buff) <- c("Longitude", "Latitude")
  for (i in 1:nrow(site)) {
    site.buff[i,1] <- site$Longitude[i] + 292000
    site.buff[i,2] <- site$Latitude[i] + 8025000  
  }
  
  buffer <- (1000/cell.size)*cell.size * 1.5
  
  for (i in 1:nrow(site.buff)) {
    if (site.buff[i,1] < 292000) {
      site.buff[i,1] <- site.buff[i,1] + buffer
    } else {
      site.buff[i,1] <- site.buff[i,1] - buffer
    }
  }
  
  for (i in 1:nrow(site.buff)) {
    if (site.buff[i,2] < 8025000) {
      site.buff[i,2] <- site.buff[i,2] + buffer
    } else {
      site.buff[i,2] <- site.buff[i,2] - buffer
    }
  }
  
  site.buff <- SpatialPoints(site.buff[ , c("Longitude","Latitude")], proj4string=CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
  
  poly <- Polygon(site.buff)
  ps <- Polygons(list(poly),1)
  non.edge.poly <- SpatialPolygons(list(ps), proj4string=CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
  
  hex_area <- hex_grid@polygons[[1]]@Polygons[[1]]@area
  mean.cell.qual <- (mean.qual / 1000000)*hex_area
  
  edge.cells <- outside.func(hex_points, non.edge.poly)
  
  for (j in 1:length(edge.cells)) {
    hex_ID <- as.numeric(edge.cells[j])
    landscape@data$graze[hex_ID] <- rtruncnorm(n=1,mean=(mean.cell.qual/3),sd=(mean.cell.qual/4),a=0)
    landscape@data$CC.graze[hex_ID] <- landscape@data$graze[hex_ID]
  }
  forage <- data.frame(landscape@data)
  return(forage)
}

assign.LIZ.locations <- function(LIZs, landscape, hex_points, cell.size) {
  mean.graze <- mean(landscape@data$graze)
  sd.graze <- sd(landscape@data$graze)
  high.threshold <- mean.graze + sd.graze 
  high.quality <- landscape@data[landscape@data$graze > high.threshold,]
  high.rows <- nrow(high.quality)
  
  for (i in 1:nrow(LIZs)) {
    rand <- round(runif(1,0.5,(high.rows+0.5)))
    id <- high.quality$hex_ID[rand]
    pos <- hex_points@coords[id,]
    radius <- (1000/cell.size)*cell.size/2
    x_shift <- runif(1,-radius,radius)
    y_shift <- runif(1,-radius,radius)
    LIZs[i,2] <- pos[1] + x_shift
    LIZs[i,3] <- pos[2] + y_shift
  }
  return(LIZs)
}

# Assign LIZs preferentially (non-randomly) to cells with high forage quality
assign.positions <- function(LIZs, landscape, hex_points, cell.size, inds, site, state.vector) {
  mean.graze <- mean(landscape@data$graze)
  sd.graze <- sd(landscape@data$graze)
  high.threshold <- mean.graze + sd.graze 
  high.quality <- landscape@data[landscape@data$graze > high.threshold,]
  high.rows <- nrow(high.quality)
  width <- max(site[,1])-min(site[,1])
  height <- max(site[,2])-min(site[,2])
  
  for (i in 1:nrow(LIZs)) {
    rand.forage <- runif(1,0,1)
    if (rand.forage < 0.7) {
      rand <- round(runif(1,0.5,(high.rows+0.5)))
      id <- high.quality$hex_ID[rand]
      pos <- hex_points@coords[id,]
      radius <- (1000/cell.size)*cell.size/2
      x_shift <- runif(1,-radius,radius)
      y_shift <- runif(1,-radius,radius)
      LIZs[i,2] <- pos[1] + x_shift
      LIZs[i,3] <- pos[2] + y_shift
    } else {
      LIZs[i,2] <- runif(1, (min(site[,1])+0.2*width), (max(site[,1])-0.2*width)) + 292000
      LIZs[i,3] <- runif(1, (min(site[,2])+0.2*height), (max(site[,2])-0.2*height)) + 8025000
    }
  }
  
  xy <- data.frame(matrix(0,N,2))
  for (j in 1:nrow(inds)) {
    if (state.vector[j,1] == 2) {
      rand <- round(runif(1,0.5,(high.rows+0.5)))
      id <- high.quality$hex_ID[rand]
      pos <- hex_points@coords[id,]
      radius <- (1000/cell.size)*cell.size/2
      x_shift <- runif(1,-radius,radius)
      y_shift <- runif(1,-radius,radius)
      xy[j,1] <- pos[1] + x_shift
      xy[j,2] <- pos[2] + y_shift
    } else {
      xy[j,1] <- runif(1, (min(site[,1])+0.2*width), (max(site[,1])-0.2*width)) + 292000
      xy[j,2] <- runif(1, (min(site[,2])+0.2*height), (max(site[,2])-0.2*height)) + 8025000
    }
  }
  xy <- cbind(seq(1,N,1), xy)
  
  position.list <- list(LIZs, xy)
  return(position.list)
}

# Assign values to agents
initialize.agents <- function(N, region) {
  
  sp.mean <- 350
  traits <- data.frame(matrix(0,N,6)) #traitframe
  m = 1
    
  # Loop over individuals of a given species
  for (j in 1:N) {
    size <- rnorm(1,sp.mean,sp.mean/8) 
    alive <- 1
    infected <- 0
    p.range <- (1000/cell.size)*cell.size + abs(rnorm(1,0.5,0.2)) * size
      
    traits[m,1] <- paste(m, 'A', sep='')
    traits[m,2] <- size
    traits[m,3] <- alive
    traits[m,4] <- infected
    traits[m,5] <- p.range
    traits[m,6] <- round(runif(1,0.5,10.5))
    m = m + 1
  }
  return(traits)
}

# Set behavioral states according to proportion observed in zebra
initialize.states <- function(N) {
  state.vector <- data.frame(matrix(0,N,1))
  for (i in 1:N) {
    rand <- runif(1,0,1)
    if (rand < 0.08007182) {
      state.vector[i,1] <- 1
    } else if (rand > 0.08007181 && rand < 0.6216637) {
      state.vector[i,1] <- 2
    } else {
      state.vector[i,1] <- 3
    }
  }
  return(state.vector)
}

# Randomly distribute LIZs
LIZ <- function(M, inds, site, t) {
  LIZs <- data.frame(matrix(0,1,5))
  m = 1
    
  for (k in 1:round(M)) { #proportion d of total individuals per year (with 2.39 being length of time until LIZ is no longer present)
    LIZs[m,1] <- paste(m, "C", sep='') #LIZ.ID
      
    LIZs[m,2] <- 0
    LIZs[m,3] <- 0
    LIZs[m,4] <- floor(runif(1,1,3.9999)) #years since death
    rand <- runif(1,0,1)
    if (rand < 0.826590) { #Probability of springbok
      LIZs[m,5] <- 1
      LIZs[m,6] <- rnorm(1,350,(350/8)) # body size
    } else if (rand > 0.826589 && rand < 0.974950) {
      LIZs[m,5] <- 2
      LIZs[m,6] <- rnorm(1,35,(35/8))
    } else {
      LIZs[m,5] <- 3
      LIZs[m,6] <- rnorm(1,3500,(3500/8))
    }
    m = m + 1
  }
    
  return(LIZs)
}

# Green up associated with LIZs
green.up <- function(landscape, LIZs) {
  
  # Increase attractivity of hex with LIZ
  coords <- coordinates(LIZs[,2:3])
  coords <- SpatialPoints(coords, proj4string=CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
    
  current_hex <- over(coords, landscape)
  for (k in 1:nrow(current_hex)) {
    idx <- current_hex[k,1]
    age <- LIZs[k,4]
    factor <- 1 + (0.4 - (age*(rnorm(1,0.1,0.02))))
    landscape@data[idx,4] <- landscape@data[idx,4] * factor
  }
  return(landscape)
}

# LIZ Layer for Probability of Contact
LIZ.layer <- function(LIZs, region) {
  
  # Increase attractivity of hex with LIZ
  coords <- coordinates(LIZs[,2:3])
  coords <- SpatialPoints(coords, proj4string=CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
  
  buffers <- c()
  for (i in 1:nrow(LIZs)) {
    age <- LIZs[i,4]
    size.factor <- log(LIZs[i,6])/log(350)
    buffers[i] <- size.factor * (rnorm(1, (3 * (4 - age)), ((3 * (4 - age))/8)))
  }

  buffers.sp <- gBuffer(coords, buffers, byid=TRUE, id = seq(1,length(buffers),1))
  df <- data.frame(matrix(seq(1,length(buffers),1), ncol=1))
  df[,2] <- rep(1,length(buffers.sp@polygons))
  colnames(df) <- c("ID", "value")
  buffers.sp <- SpatialPolygonsDataFrame(buffers.sp, data = df)
  
  area.extent <- extent(region)
  r <- raster(area.extent)
  projection(r) <- CRS("+proj=utm +south +zone=35 +ellps=WGS84")
  res(r) <- c(10,10)
  
  strt <- Sys.time()
  LIZ.layer <- rasterize(buffers.sp, r, background=0, getCover=TRUE)
  print(Sys.time() - strt)
  
  return(LIZ.layer)
}

#### Agent State Processes #####

# Movement of individuals based on behavioral states (parameterized using zebra data)
move.func.simp <- function(N, inds, landscape, hex_points, hex_grid, state.vector, TimeSinceSunrise, LIZs, site, prev.angles) {
  
  moves <- data.frame(matrix(0,N,3))
  tracker <- data.frame(matrix(0,N,3))
  
  points <- coordinates(inds[,7:8])
  coords <- SpatialPoints(points, proj4string=CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
  curr.cells <- as.numeric(over(coords, hex_grid))
  
  #moves <- foreach(j = 1:N, .packages=c("sp", "rgeos", "prevR", 'fdrtool'), .export=c("movement", "anglefun", "bounce", "within.func"), .combine=rbind) %dopar% {
  for (j in 1:N) {    
    
    if (state.vector[j,1] == 2) {
      coords.temp <- SpatialPoints(coords[j,], proj4string=CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
      perception <- gBuffer(coords.temp, width=inds[j,5], byid=TRUE)
    
      within <- within.func(hex_points, perception)
      curr.cell <- curr.cells[j]
      curr.quality <- landscape@data[curr.cell,3]
    
      if (length(within) >= 1) {
        nearby <- c()
        for (k in 1:length(within)) {
          hex_id <- within[k]
          nearby[k] <- landscape@data[hex_id,3]
        }
        max.attract <- max(nearby)
        max.cell <- within[which.max(nearby)]
      } else {
        max.cell <- curr.cell
      }
      
      if (curr.cell == max.cell) {
        cell_center <- data.frame(matrix(hex_points@coords[curr.cell,], nrow=1))
        x_diff <- as.numeric(points[j,1] - cell_center[1,1])
        y_diff <- as.numeric(points[j,2] - cell_center[1,2])
        heading <- as.numeric(anglefun(x_diff, y_diff))
      } else {
        max_center <- data.frame(matrix(hex_points@coords[max.cell,], nrow=1))
        x_diff <- as.numeric(points[j,1] - max_center[1,1])
        y_diff <- as.numeric(points[j,2] - max_center[1,2])
        heading <- as.numeric(anglefun(x_diff, y_diff))
      }
      
      stepsize <- rgamma(1, shape=((106^2)/(13^2)), rate=(106/(13^2))) #rtruncnorm(n=1, mean=106, sd=13, a=0)
      
      tracker[j,1:3] <- c(stepsize, heading, state.vector[j,1])
      temp.pts <- data.frame(matrix(points[j,], nrow=1))
      moves[j,1:3] <- movement(temp.pts, stepsize, heading, TimeSinceSunrise)
      
    } else if (state.vector[j,1] == 1) {
      stepsize <- rgamma(1, shape=((8.619^2)/(7.840^2)), rate=(8.619/(7.840^2))) #rtruncnorm(n=1, mean=8.619, sd=7.840, a=0)
      heading.adjust <- rvm(n=1, mean=3.119, k=0.489)
      old.heading <- as.numeric(prev.angles[j,1])
      heading <- old.heading + heading.adjust
      heading[heading < 0] <- heading[heading < 0] + 2 * pi
      
      tracker[j,1:3] <- c(stepsize, heading, state.vector[j,1])
      temp.pts <- data.frame(matrix(points[j,], nrow=1))
      moves[j,1:3] <- movement(temp.pts, stepsize, heading, TimeSinceSunrise)
      
    } else {
      stepsize <- rgamma(1, shape=((628^2)/(403^2)), rate=(628/(403^2))) #rtruncnorm(n=1, mean=628, sd=403, a=0)
      heading.adjust <- rvm(n=1, mean=0.0004, k=3.163)
      old.heading <- as.numeric(prev.angles[j,1])
      heading <- old.heading + heading.adjust
      heading[heading < 0] <- heading[heading < 0] + 2 * pi
        
      tracker[j,1:3] <- c(stepsize, heading, state.vector[j,1])
      temp.pts <- data.frame(matrix(points[j,], nrow=1))
      moves[j,1:3] <- movement(temp.pts, stepsize, heading, TimeSinceSunrise)
    }
    
    #plot(perception)
    #points(temp.pts[,1], temp.pts[,2], pch=19)
    #points(hex_points@coords[curr.cell,1], hex_points@coords[curr.cell,2], pch=19, col='blue')
    #points(move[,1], move[,2], pch=19, col='red')
    
  }
  
  for (r in 1:nrow(inds)) {
    while (point.in.SpatialPolygons(moves[r,1], moves[r,2], region) == FALSE) {
      old.pts <- data.frame(matrix(points[r,], nrow=1))
      stepsize <- (1.5 * rgamma(1, shape=((628^2)/(403^2)), rate=(628/(403^2))))
      old.heading <- tracker[r,2]
      if ((old.heading + pi) > (2*pi)) {
        heading <- (old.heading - pi)
      } else {
        heading <- (old.heading + pi)
      }
      moves[r,1:3] <- movement(old.pts, stepsize, heading, TimeSinceSunrise)
      tracker[r,1:3] <- c(stepsize, heading, state.vector[j,1])
    }
  }
  
  #for (q in 1:nrow(inds)) {
  #  moves[q,4:6] <- dist.to.LIZ(moves[q,1:2], LIZs, 1, site)
  #  moves[q,7] <- prob.contact(moves[q,4], 1/((moves[q,6]/317.2266)*rnorm(1,3,0.5)))
  #}
  
  moves.list <- list(moves, tracker)
  return(moves.list)
}

# Transition individuals between states according the transition probabilities from moveHMM
state.vector <- function(N, prev.vector, TimeSinceSunrise) {
  state.vector <- data.frame(matrix(0,N,1))
  if (TimeSinceSunrise == 23) {
    prop1 <- sum(prev.vector[,1] == 1)
    prop2 <- sum(prev.vector[,1] == 2)
    prop3 <- sum(prev.vector[,1] == 3)
    diff <- round(0.3783363*N) - (prop3)
    if ((diff/N) > 0) {
      needed.state3 <- (diff/N)
    } else {
      needed.state3 <- 0
    }
    shift.2to3 <- rbinom(prop2,1,needed.state3)
    shift.1to2 <- rbinom(prop1,1,needed.state3)
    counter1 <- 1
    counter2 <- 1
    for (j in 1:nrow(prev.vector)) {
      if (prev.vector[j,1] == 1) {
        if (shift.1to2[counter1] == 1) {
          state.vector[j,1] <- 2
        } else {
          state.vector[j,1] <- 1
        }
        counter1 <- counter1 + 1
      } else if (prev.vector[j,1] == 2) {
        if (shift.2to3[counter2] == 1) {
          state.vector[j,1] <- 3
        } else {
          state.vector[j,1] <- 2
        }
        counter2 <- counter2 + 1
      } else {
        state.vector[j,1] <- 3
      }
    }
  } else {
    for (i in 1:N) {
      if (prev.vector[i,1] == 1) {
        prob.1to2 <- exp(-1.48222568 + (0.029693 * TimeSinceSunrise))
        trans <- rbinom(1,1,prob.1to2)
        if (trans == 0) {
          state.vector[i,1] <- 1
        } else {
          state.vector[i,1] <- 2
        }
      } else if (prev.vector[i,1] == 3) {
        prob.3to2 <- exp(-1.4800534 + (0.023695 * TimeSinceSunrise))
        trans <- rbinom(1,1,prob.3to2)
        if (trans == 0) {
          state.vector[i,1] <- 3
        } else {
          state.vector[i,1] <- 2
        }
      } else {
        prob.2to1 <- exp(-2.66930802 + (0.03868 * TimeSinceSunrise))
        prob.2to3 <- exp(-2.317976802 + (0.001909 * TimeSinceSunrise))
        trans.2to1 <- rbinom(1,1,prob.2to1)
        trans.2to3 <- rbinom(1,1,prob.2to3)
        while (trans.2to1 == 1 && trans.2to3 == 1) {
          trans.2to1 <- rbinom(1,1,prob.2to1)
          trans.2to3 <- rbinom(1,1,prob.2to3)
        }
        if (trans.2to1 == 0 && trans.2to3 == 0) {
          state.vector[i,1] <- 2
        } else if (trans.2to1 == 1 && trans.2to3 == 0) {
          state.vector[i,1] <- 1
        } else if (trans.2to1 == 0 && trans.2to3 == 1) {
          state.vector[i,1] <- 3
        }
      }
    }
  }
  return(state.vector)
}
  
#### Landscape Processes #####

# Individuals feeding decreases the quality of the cell
feeding <- function(N, inds, landscape, state.vector) {
  coords <- inds[,7:8]
  coords <- coordinates(coords)
  colnames(coords) <- c("Longitude", "Latitude")
  coords <- SpatialPoints(coords[ , c("Longitude","Latitude")], proj4string=CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
  
  curr.pos <- over(coords, landscape)[,1:4]
  for (j in 1:N) {
    if (state.vector[j,1] == 2) {
      if (!is.na(curr.pos[j,1])) {
        hex <- curr.pos[j,1]
        if ((landscape@data[hex,3] - (inds[j,2]*12/(365*(24*0.5415919)))) > 0) {
          landscape@data[hex,3] <- landscape@data[hex,3] - (inds[j,2]*12/(365*(24*0.5415919)))
        } else {
          landscape@data[hex,3] <- landscape@data[hex,3]
        }
      }
    }
  }
  return(landscape)
}

disperse.inds <- function(balls, bins) {
  bin.df <- data.frame(matrix(0,bins,1))
  for (i in 1:balls) {
    rand <- round(runif(1,0.5,(bins+0.5)))
    bin.df[rand,1] <- bin.df[rand,1] + 1
  }
  return(bin.df)
}

general.feeding <- function(density, landscape, mean.qual, region) {
  K.high <- round(density * (region@polygons[[1]]@area/1000000) * 0.5415919 * 0.8)
  K.low <- round(density * (region@polygons[[1]]@area/1000000) * 0.5415919 * 0.2)
  if (max(land@data[,3]) > mean.qual) {
    high.quality <- which(land@data[,3] > mean.qual)
    high.decrease <- disperse.inds(K.high, length(high.quality))
    for (i in 1:length(high.quality)) {
      hex <- high.quality[i]
      if (high.decrease[i,1] > 0) {
        if ((landscape@data[hex,3] - (sum(rnorm(high.decrease[i,1],350,(350/8))*12)/(365*(24*0.5415919)))) > 0) {
          landscape@data[hex,3] <- landscape@data[hex,3] - (sum(rnorm(high.decrease[i,1],350,(350/8))*12)/(365*(24*0.5415919)))
        } else {
          landscape@data[hex,3] <- landscape@data[hex,3]
        }
      }
    }
  } else {
    K.low <- round(density * (region@polygons[[1]]@area/1000000) * 0.5415919)
  }
  low.quality <- which(land@data[,3] < mean.qual)
  low.decrease <- disperse.inds(K.low, length(low.quality))
  for (j in 1:length(low.quality)) {
    hex2 <- low.quality[j]
    if (low.decrease[j,1] > 0) {
      if ((landscape@data[hex2,3] - (sum(rnorm(low.decrease[j,1],350,(350/8))*12)/(365*24))) > 0) {
        landscape@data[hex2,3] <- landscape@data[hex2,3] - (sum(rnorm(low.decrease[j,1],350,(350/8))*12)/(365*(24*0.5415919)))
      } else {
        landscape@data[hex2,3] <- landscape@data[hex2,3]
      }
    }
  }
  return(landscape)
}

# Background (logistic) growth rate of vegetation
veg.growth <- function(landscape, veg.rate) {
  for (i in 1:nrow(landscape@data)) {
    if (landscape@data[i,3] < landscape@data[i,4]) {
      if (landscape@data[i,3] > 0) {
        landscape@data[i,3] <- landscape@data[i,3] + (veg.rate*landscape@data[i,3]) * (1 - landscape@data[i,3]/landscape@data[i,4])
      } else {
        landscape@data[i,3] <- 0.001 * landscape@data[i,4]
      }
    }
  }  
  return(landscape)
}

###### Post Run ######

save.paths <- function(all.steps, N, run.number, steps.per.day) {
  all.paths <- list()
  for (y in 1:N) {
    ind.steps <- data.frame(matrix(0,length(all.steps),7))
    start.day <- as.Date("2010/02/01", format="%Y/%m/%d")
    for (x in 1:length(all.steps)) {
      pos <- data.frame(matrix(all.steps[[x]][y,], nrow=1))
      tracker <- data.frame(matrix(tracker.list[[x]][y,], nrow=1))
      ind.steps[x,1] <- pos[,1]
      ind.steps[x,2] <- pos[,2]
      ind.steps[x,3] <- tracker[,1]
      ind.steps[x,4] <- tracker[,2]
      ind.steps[x,5] <- tracker[,3]
      ind.steps[x,6] <- y
      day <- x%/%steps.per.day
      date <- start.day + day
      hour <- ((x - ((day-1)*steps.per.day))%/%(steps.per.day/24)) - 24
      hour.char <- if(hour<10){paste("0",hour,sep='')}else{as.character(hour)}
      min <- (x - ((day-1)*steps.per.day))%%(steps.per.day/24)
      ind.steps[x,7] <- paste(date," ",hour.char,":",if(min==0){"00"}else{(min*20)}, sep='')
    }
    all.paths[[y]] <- data.frame(ind.steps)
    write.csv(ind.steps, paste("ID",y,"_Steps_Run",run.number,".csv", sep=''))
  }
  return(all.paths)
}

###################################################################

run.num <- 13
site <- data.frame(matrix(c(-25000, 25000, 25000, -25000, 25000, 25000, -25000, -25000),4,2)) #2500 km2
#site <- data.frame(matrix(c(-35400, 35400, 35400, -35400, 35400, 35400, -35400, -35400),4,2)) #~5000 km2
cell.size=200
mean.qual <- runif(1,900000,1100000)
colnames(site) <- c("Longitude", "Latitude")

# Initialize landscape
#land.temp <- landscape(site, cell.size)
land.temp <- landscape.hetero(site, cell.size, prop.graze=0.02, mean.qual=mean.qual)

hex_points <- land.temp[[2]]
hex_grid <- land.temp[[3]]
land <- land.temp[[1]]
forage <- edge.adjustment(site, land, cell.size, hex_points, hex_grid, mean.qual)
land@data <- forage
region <- gUnaryUnion(land, id = land@data$region)
#spplot(land[,"graze"])

N = 20
M = round(0.01353528 * (region@polygons[[1]]@area/1000000) * 3)
steps.per.day <- 72
days <- 90
agents <- initialize.agents(N, region)
init.state <- data.frame(matrix(0,N,1))
behav.states <- data.frame(matrix(0,N,days*(steps.per.day/3)))
init.state <- initialize.states(N)
LIZs <- LIZ(M, agents, site, 1)
position.list <- assign.positions(LIZs, land, hex_points, cell.size, agents, site, init.state)
LIZs <- position.list[[1]]
xy <- position.list[[2]]
land <- green.up(land, LIZs)
#p <- p + layer(panel.points(LIZs[,2], LIZs[,3], col="white", pch=19))
#p + layer(panel.points(xy[,2], xy[,3], col="black", pch=19))
#dist.mat <- dist.to.LIZ(xy[,2:3], LIZs, N, site)
bound <- cbind(agents, xy[,2:3])
colnames(bound) <- c("id", "size", "alive", "infected", "p.range", "num", "x", "y")

#LIZ.layer <- LIZ.layer(LIZs, region)
#writeRaster(LIZ.layer, paste0("LIZ_Layer_Run",run.num,".csv"))
starting.land <- land
starting.LIZs <- LIZs
starting.agents <- bound
save.image("~/Desktop/Ecol_Letters_ABM/Anthrax_Sim_Workspace_v9.RData")

#cumulative.contacts <- list()
all.steps <- list()
landscape.list <- list()
all.steps[[1]] <- xy[,2:3]
#daily.prob <- list()
tracker.list <- list()
init.angles <- data.frame(matrix(rep(0,N),ncol=1))

#plot(region)
#points(xy[,2], xy[,3], pch=19)
#points(LIZs[,2], LIZs[,3], pch=19, col='blue')

strt <- Sys.time()
for (z in 1:days) {
  for (t in 1:(steps.per.day + 1)) {
    TimeSinceSunrise <- (t-1)%/%3
    step.no <- t%%3
    if (step.no == 1) {
      land <- general.feeding(0.89, land, mean.qual, region)
      if (z == 1 && TimeSinceSunrise == 0) {
        prev.state <- init.state
        behav.states[,1] <- state.vector(N, prev.state, TimeSinceSunrise)
        prev.angles <- init.angles
      } else {
        prev.state <- data.frame(matrix(behav.states[,(((z-1)*24) + TimeSinceSunrise)]))
        behav.states[,(((z-1)*24) + (TimeSinceSunrise+1))] <- state.vector(N, prev.state, TimeSinceSunrise)
      } 
    }
    current.states <- data.frame(matrix(behav.states[,(((z-1)*24) + (TimeSinceSunrise+1))]))
    moves <- move.func.simp(N, bound, land, hex_points, hex_grid, current.states, TimeSinceSunrise, LIZs, site, prev.angles)
    step <- moves[[1]][,1:2]
    all.steps[[((z-1)*steps.per.day) + t]] <- step
    #daily.prob[[t]] <- moves[[1]][,7]
    bound <- cbind(agents, moves[[1]][,c(1:2)])
    tracker.list[[((z-1)*steps.per.day) + t]] <- moves[[2]]
    prev.angles <- data.frame(matrix(moves[[2]][,2], ncol=1))
    #points(step[,1], step[,2], pch=19, col='red', cex=0.3)
    if (step.no == 0) {
      land <- veg.growth(land, veg.rate=0.1)
    }
    if (step.no == 2) {
      land <- feeding(N, bound, land, current.states)
    }
  }
  landscape.list[[z]] <- data.frame(land@data)
  #cumulative.contacts[[z]] <- cumulative.contact(daily.prob, N, steps.per.day)
  #total.contact <- sum(cumulative.contacts[[z]][,1])
  cat("Day", z, "\n")
}
print(Sys.time() - strt)

all.paths <- save.paths(all.steps, N, run.num, steps.per.day)

for (i in 1:length(landscape.list)) {
  write.csv(landscape.list[[i]], paste("Landscape_Data_Day_",i,".csv",sep=''))
}
saveRDS(landscape.list, paste0("Landscape_Run",run.num,".rds"))

coords <- coordinates(LIZs[,2:3])
coords <- SpatialPoints(coords, proj4string=CRS("+proj=utm +south +zone=35 +ellps=WGS84"))

buffers <- c()
for (i in 1:nrow(LIZs)) {
  age <- LIZs[i,4]
  size.factor <- log(LIZs[i,6])/log(350)
  buffers[i] <- size.factor * (rnorm(1, (3 * (4 - age)), ((3 * (4 - age))/8)))
}

buffer.LIZ <- cbind(LIZs, buffers)
write.csv(buffer.LIZ, paste0("LIZs_Buffers_Run",run.num,".csv"))

buffers.sp <- gBuffer(coords, buffers, byid=TRUE, id = seq(1,length(buffers),1))
df <- data.frame(matrix(seq(1,length(buffers),1), ncol=1))
df[,2] <- rep(1,length(buffers.sp@polygons))
colnames(df) <- c("ID", "value")
buffers.sp <- SpatialPolygonsDataFrame(buffers.sp, data = df)