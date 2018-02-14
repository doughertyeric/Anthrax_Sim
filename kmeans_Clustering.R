library(lubridate)
library(data.table)
library(dplyr)

run.num = 13
temp = list.files(pattern="*_Steps_Run13.csv")
ids <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
movementpaths <- lapply(temp, read.csv)
movementpaths <- lapply(movementpaths, setNames, c("stepID", "X", "Y", "stepsize", 
                                                   "angle", "mode", "ID", "Time"))
movementpaths <- lapply(movementpaths, mutate, Time = ymd_hm(Time))

paths <- list()
for (i in 1:20){
  clust <- movementpaths[[i]]
  
  for (j in 1:nrow(clust)) {
    if (clust$mode[j] == 1 && clust$stepsize[j] > 100) {
      clust$mode[j] <- 3
    } 
  }
  write.csv(clust, paste0("ID",i,"_Steps_Run",run.num,".csv"))
  
  kdist <- kmeans(clust[, c("stepsize")], centers= c(0,150,1000))
  clust$kmeans <- kdist$cluster
  #sorttable <- summarise(group_by(clust, kmeans), med.dist = median(stepsize))
  #clust$modek <- ifelse(clust$kmeans == as.numeric(sorttable[sorttable$med.dist==min(sorttable$med.dist), "kmeans"]), 1,
  #                      ifelse(clust$kmeans == as.numeric(sorttable[sorttable$med.dist==max(sorttable$med.dist), "kmeans"]), 3, 2))
  paths[[i]] <- clust
  print(paste0("-----", i, " ----"))
  #print(sorttable)
  #print(table(clust$kmeans))
  #print(table(clust$modek))
  print(table(clust[clust$mode == 2, "kmeans"])/sum(table(clust[clust$mode == 2, "kmeans"])))
}

p <- rbindlist(paths)

for (i in 1:20) {
  temp <- subset(p, ID == i)
  write.csv(temp, paste0("ID",i,"_kmeans_Run",run.num,".csv"))
}

# Summary stats across all simulants

for (i in 1:3) {
  temp <- subset(p, mode == i)
  print(paste("Median",i,median(temp$stepsize)))
  print(paste("Mean",i,mean(temp$stepsize)))
  print(paste("SD",i,sd(temp$stepsize)))
}

####################

temp = list.files(pattern="*_kmeans_Run13.csv")
ids <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
movement.paths2 <- lapply(temp, read.csv)
movement.paths2 <- lapply(movement.paths2, setNames, c("X.1", "stepID", "X", "Y", "stepsize", 
                                                   "angle", "mode", "ID", "Time", "kmeans"))
movement.paths2 <- lapply(movement.paths2, mutate, Time = ymd_hms(Time))

all.paths.temp <- data.frame(matrix(0,1,9))
colnames(all.paths.temp) <- c("stepID", "X", "Y", "stepsize", "angle", "mode", "ID", "Time", "kmeans")
for (i in 1:length(movement.paths2)) {
  move1 <- movement.paths2[[i]]
  move1 <- move1[,-1]
  colnames(move1) <- c("stepID", "X", "Y", "stepsize", "angle", "mode", "ID", "Time", "kmeans")
  all.paths.temp <- data.frame(rbind(all.paths.temp, move1))
}
all.paths.temp <- all.paths.temp[-1,]

count <- 0
for (i in 1:nrow(all.paths.temp)) {
  if (all.paths.temp$mode[i] == all.paths.temp$kmeans[i]) {
    count <- count + 1
  }
}
print(count/nrow(all.paths.temp))

foraging <- all.paths.temp[all.paths.temp$mode == 2,]
count2 <- 0
for (i in 1:nrow(foraging)) {
  if (foraging$kmeans[i] == 2) {
    count2 <- count2 + 1
  }
}
print(count2/nrow(foraging))
