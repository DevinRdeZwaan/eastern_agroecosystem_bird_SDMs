##########################################################################################
### R code for:

### "Balancing conservation priorities for grassland and forest specialist bird
### communities in agriculturally dominated landscapes"

###  DR de Zwaan, N Alavi, GW Mitchell, DR Lapen, J Duffe, & S Wilson (2021)
###  Biological Conservation	

### Run with R version 3.6.3
##########################################################################################  


#################################
### Code description and notices:
#################################

### The following code:
### 1) Pre-processes BBS count data (already assembled)
### 2) Fits presence/absence ensemble models for species-at-risk (SAR) and
###     species of special concern (SOSC)
### 3) Fits stacked species distribution models (SSDM)
### 4) Evaluates model fit

### Important note:
### The original package used to create ensemble models (package 'SSDM') no longer
### properly calculates omission rate, sensitivity, specificity, or AUC in the
### former 'evaluation' process (at least at the time of this analysis).
### As a result, I have rewritten the code to calculate these by hand.


################################################################################
### Begin code
################################################################################

### Set working directory

setwd("")

### Load all relevant packages
library(tidyr)
library(ggmap)
library(maps)
library(mapdata)
library(raster)
library(rgdal)
library(gdalUtils)
library(sp)
library(raster)
library(plyr)
library(dplyr)
library(ggplot2)
library(landscapemetrics)
library(tidyverse)
library(vegan)
library(viridis)
library(SSDM)
library(spThin)
library(ROCR)
library(mgcv)
library(gratia)

### Increase raster memory to 16 gb
options(rasterMaxMemory = 16e9)
memory.limit(size = NA)
memory.limit(size = 16e9)


################################################################################
### Pre-process bird count data
################################################################################

### Read in bird data
bird_data <- read.csv("bird_data_2014_2018.csv")

head(bird_data)

### Subset 2018 only

bird_data_2018 <- subset(bird_data, year==2018)

### Select location and abundance data for SAR/SOSC

bird_data_use <- bird_data_2018 %>% 
                 select(longitude,latitude,BOBO,EAME,EWPE,WOTH,HOLA,SAVS,KILL,LEFL)


### Convert abundances to occurrences
bird_data_use$BOBO_pa <- ifelse(bird_data_use$BOBO >0, 1,0) 
bird_data_use$EAME_pa <- ifelse(bird_data_use$EAME >0, 1,0)
bird_data_use$EWPE_pa <- ifelse(bird_data_use$EWPE >0, 1,0)
bird_data_use$WOTH_pa <- ifelse(bird_data_use$WOTH >0, 1,0)
bird_data_use$HOLA_pa <- ifelse(bird_data_use$HOLA >0, 1,0)
bird_data_use$SAVS_pa <- ifelse(bird_data_use$SAVS >0, 1,0)
bird_data_use$KILL_pa <- ifelse(bird_data_use$KILL >0, 1,0)
bird_data_use$LEFL_pa <- ifelse(bird_data_use$LEFL >0, 1,0)

### Remove original count data

bird_data_use <- bird_data_use[,-c(3:10)] 

### Create presence only dataset for each species.

BOBO_presence_data <- subset(bird_data_use, BOBO_pa > 0)
EAME_presence_data <- subset(bird_data_use, EAME_pa > 0)
EWPE_presence_data <- subset(bird_data_use, EWPE_pa > 0)
WOTH_presence_data <- subset(bird_data_use, WOTH_pa > 0)
HOLA_presence_data <- subset(bird_data_use, HOLA_pa > 0)
SAVS_presence_data <- subset(bird_data_use, SAVS_pa > 0)
KILL_presence_data <- subset(bird_data_use, KILL_pa > 0)
LEFL_presence_data <- subset(bird_data_use, LEFL_pa > 0)


### Remove other species data

BOBO_presence_data <- BOBO_presence_data[,-c(4:10)]
EAME_presence_data <- EAME_presence_data[,-c(3,5:10)]
EWPE_presence_data <- EWPE_presence_data[,-c(3,4,6:10)]
WOTH_presence_data <- WOTH_presence_data[,-c(3:5,7:10)]
HOLA_presence_data <- HOLA_presence_data[,-c(3:6,8:10)]
SAVS_presence_data <- SAVS_presence_data[,-c(3:7,9:10)]
KILL_presence_data <- KILL_presence_data[,-c(3:8,10)]
LEFL_presence_data <- LEFL_presence_data[,-c(3:9)]


####################################
### Occurrence sample sizes for 2018
####################################

length(BOBO_presence_data$BOBO_pa) #421
length(EAME_presence_data$EAME_pa) #323
length(HOLA_presence_data$HOLA_pa) #111
length(SAVS_presence_data$SAVS_pa) #694
length(KILL_presence_data$KILL_pa) #307
length(WOTH_presence_data$WOTH_pa) #195
length(EWPE_presence_data$EWPE_pa) #202
length(LEFL_presence_data$LEFL_pa) #158


################################################
### Create absence data in 2018 for each species
################################################

### Subset data
BOBO_p_a <- bird_data_use %>% select(longitude,latitude,BOBO_pa)
EAME_p_a <- bird_data_use %>% select(longitude,latitude,EAME_pa)
HOLA_p_a <- bird_data_use %>% select(longitude,latitude,HOLA_pa)
SAVS_p_a <- bird_data_use %>% select(longitude,latitude,SAVS_pa)
KILL_p_a <- bird_data_use %>% select(longitude,latitude,KILL_pa)
WOTH_p_a <- bird_data_use %>% select(longitude,latitude,WOTH_pa)
EWPE_p_a <- bird_data_use %>% select(longitude,latitude,EWPE_pa)
LEFL_p_a <- bird_data_use %>% select(longitude,latitude,LEFL_pa)

### Add species column
BOBO_p_a$species <- rep("BOBO", length(BOBO_p_a$BOBO))
EAME_p_a$species <- rep("EAME", length(EAME_p_a$EAME))
HOLA_p_a$species <- rep("HOLA", length(HOLA_p_a$HOLA))
SAVS_p_a$species <- rep("SAVS", length(SAVS_p_a$SAVS))
KILL_p_a$species <- rep("KILL", length(KILL_p_a$KILL))
WOTH_p_a$species <- rep("WOTH", length(WOTH_p_a$WOTH))
EWPE_p_a$species <- rep("EWPE", length(EWPE_p_a$EWPE))
LEFL_p_a$species <- rep("LEFL", length(LEFL_p_a$LEFL))

### Select only 0s
BOBO_absences <- subset(BOBO_p_a, BOBO_pa == 0)
EAME_absences <- subset(EAME_p_a, EAME_pa == 0)
HOLA_absences <- subset(HOLA_p_a, HOLA_pa == 0)
SAVS_absences <- subset(SAVS_p_a, SAVS_pa == 0)
KILL_absences <- subset(KILL_p_a, KILL_pa == 0)
WOTH_absences <- subset(WOTH_p_a, WOTH_pa == 0)
EWPE_absences <- subset(EWPE_p_a, EWPE_pa == 0)
LEFL_absences <- subset(LEFL_p_a, LEFL_pa == 0)


########################
### Spatial thinning
########################

### Note:
### Thin presence and absence data separately for each species to reduce 
### spatial autocorrelation


#############
### Presences
#############

### Add species column to each dataset
BOBO_presence_data$species <- rep("BOBO", length(BOBO_presence_data$longitude))
EAME_presence_data$species <- rep("EAME", length(EAME_presence_data$longitude))
HOLA_presence_data$species <- rep("HOLA", length(HOLA_presence_data$longitude))
SAVS_presence_data$species <- rep("SAVS", length(SAVS_presence_data$longitude))
KILL_presence_data$species <- rep("KILL", length(KILL_presence_data$longitude))
WOTH_presence_data$species <- rep("WOTH", length(WOTH_presence_data$longitude))
EWPE_presence_data$species <- rep("EWPE", length(EWPE_presence_data$longitude))
LEFL_presence_data$species <- rep("LEFL", length(LEFL_presence_data$longitude))

### Make vector of datasets and empty list to fill

presence_dfs <- list(BOBO_presence_data,EAME_presence_data,HOLA_presence_data,
                     SAVS_presence_data,KILL_presence_data,WOTH_presence_data,
                     EWPE_presence_data,LEFL_presence_data)

presence_thinned_list <- list()

### Randomly thin by 1km (to avoid adjacent occurrences)

set.seed(656)

for(i in 1:8) {
  
  presence_thinned_list[[i]] <- thin(presence_dfs[[i]], 
                                   long.col = "longitude", lat.col = "latitude",
                                   spec.col = "species", thin.par = 1, reps = 10, #1km and 10 replication
                                   locs.thinned.list.return = TRUE,
                                   write.files = FALSE, write.log.file = FALSE, 
                                   verbose = TRUE)
}

### Extract species-specific presences

BOBO_presence_thinned <- data.frame(presence_thinned_list[[1]][[1]])
EAME_presence_thinned <- data.frame(presence_thinned_list[[2]][[1]])
HOLA_presence_thinned <- data.frame(presence_thinned_list[[3]][[1]])
SAVS_presence_thinned <- data.frame(presence_thinned_list[[4]][[1]])
KILL_presence_thinned <- data.frame(presence_thinned_list[[5]][[1]])
WOTH_presence_thinned <- data.frame(presence_thinned_list[[6]][[1]])
EWPE_presence_thinned <- data.frame(presence_thinned_list[[7]][[1]])
LEFL_presence_thinned <- data.frame(presence_thinned_list[[8]][[1]])

### Create presence column

BOBO_presence_thinned$BOBO_pa <- rep(1,length(BOBO_presence_thinned$Longitude))
EAME_presence_thinned$EAME_pa <- rep(1,length(EAME_presence_thinned$Longitude))
HOLA_presence_thinned$HOLA_pa <- rep(1,length(HOLA_presence_thinned$Longitude))
SAVS_presence_thinned$SAVS_pa <- rep(1,length(SAVS_presence_thinned$Longitude))
KILL_presence_thinned$KILL_pa <- rep(1,length(KILL_presence_thinned$Longitude))
WOTH_presence_thinned$KILL_pa <- rep(1,length(WOTH_presence_thinned$Longitude))
EWPE_presence_thinned$KILL_pa <- rep(1,length(EWPE_presence_thinned$Longitude))
LEFL_presence_thinned$KILL_pa <- rep(1,length(LEFL_presence_thinned$Longitude))

### Rename columns

colnames(BOBO_presence_thinned) <- c("longitude","latitude","BOBO_pa")
colnames(EAME_presence_thinned) <- c("longitude","latitude","EAME_pa")
colnames(HOLA_presence_thinned) <- c("longitude","latitude","HOLA_pa")
colnames(SAVS_presence_thinned) <- c("longitude","latitude","SAVS_pa")
colnames(KILL_presence_thinned) <- c("longitude","latitude","KILL_pa")
colnames(WOTH_presence_thinned) <- c("longitude","latitude","WOTH_pa")
colnames(EWPE_presence_thinned) <- c("longitude","latitude","EWPE_pa")
colnames(LEFL_presence_thinned) <- c("longitude","latitude","LEFL_pa")



#################
### Absences
#################

### Thin absences by looping through each species

absence_thinned_list <- list()
species_list <- list(BOBO_absences,EAME_absences,HOLA_absences,SAVS_absences,
                     KILL_absences,WOTH_absences,EWPE_absences,LEFL_absences)

set.seed(111)

for(i in 1:8) {
  
  absence_thinned_list[[i]] <- thin(species_list[[i]], long.col = "longitude", 
                                  lat.col = "latitude",spec.col = "species", 
                                  thin.par = 1, reps = 10, 
                                  locs.thinned.list.return = TRUE,
                                  write.files = FALSE, write.log.file = FALSE, verbose = TRUE)
}

### Create dataframes from lists

BOBO_absences_thin_df <- data.frame(absence_thinned_list[[1]][[1]])
EAME_absences_thin_df <- data.frame(absence_thinned_list[[2]][[1]])
HOLA_absences_thin_df <- data.frame(absence_thinned_list[[3]][[1]])
SAVS_absences_thin_df <- data.frame(absence_thinned_list[[4]][[1]])
KILL_absences_thin_df <- data.frame(absence_thinned_list[[5]][[1]])
WOTH_absences_thin_df <- data.frame(absence_thinned_list[[6]][[1]])
EWPE_absences_thin_df <- data.frame(absence_thinned_list[[7]][[1]])
LEFL_absences_thin_df <- data.frame(absence_thinned_list[[8]][[1]])

### Create column of zeros to match 1's in presence data
BOBO_absences_thin_df$BOBO_pa <- rep(0,length(BOBO_absences_thin_df$Longitude))
EAME_absences_thin_df$EAME_pa <- rep(0,length(EAME_absences_thin_df$Longitude))
HOLA_absences_thin_df$HOLA_pa <- rep(0,length(HOLA_absences_thin_df$Longitude))
SAVS_absences_thin_df$SAVS_pa <- rep(0,length(SAVS_absences_thin_df$Longitude))
KILL_absences_thin_df$KILL_pa <- rep(0,length(KILL_absences_thin_df$Longitude))
WOTH_absences_thin_df$WOTH_pa <- rep(0,length(WOTH_absences_thin_df$Longitude))
EWPE_absences_thin_df$EWPE_pa <- rep(0,length(EWPE_absences_thin_df$Longitude))
LEFL_absences_thin_df$LEFL_pa <- rep(0,length(LEFL_absences_thin_df$Longitude))

### Rename column names to match
colnames(BOBO_absences_thin_df) <- c("longitude","latitude","BOBO_pa")
colnames(EAME_absences_thin_df) <- c("longitude","latitude","EAME_pa")
colnames(HOLA_absences_thin_df) <- c("longitude","latitude","HOLA_pa")
colnames(SAVS_absences_thin_df) <- c("longitude","latitude","SAVS_pa")
colnames(KILL_absences_thin_df) <- c("longitude","latitude","KILL_pa")
colnames(WOTH_absences_thin_df) <- c("longitude","latitude","WOTH_pa")
colnames(EWPE_absences_thin_df) <- c("longitude","latitude","EWPE_pa")
colnames(LEFL_absences_thin_df) <- c("longitude","latitude","LEFL_pa")



#############################################################################
### Thin HOLA in a spatially structured way
#############################################################################

### Note:
### Because of the low sample size and bias towards southern Ontario, HOLA
### appears to over predict in the south and under predict in the north when
### thinning is done evenly across the range (preliminary attempts).

### To counter this, we divided HOLA occurrences into two halves: above and 
### below 42.75 deg latitude. This latitude was based on evenly distributing
### the density of occurrences.
### Those below this latitude were thinned, those above were not.

### Divide occurrences by latitude
HOLA_presence_sub_42 <- subset(HOLA_presence_data, latitude <= 42.75)
HOLA_presence_above_42 <- subset(HOLA_presence_data, latitude >= 42.75)

### Thin
set.seed(22)
HOLA_thinned_sub_42 <- thin(HOLA_presence_sub_42, 
                            long.col = "longitude", lat.col = "latitude",
                            spec.col = "species", thin.par = 1, reps = 10, 
                            locs.thinned.list.return = TRUE,
                            write.files = FALSE, write.log.file = FALSE, 
                            verbose = TRUE)

HOLA_presence_thinned_sub_42 <- data.frame(HOLA_thinned_sub_42[[1]])
HOLA_presence_thinned_sub_42$HOLA_pa <- rep(1,length(HOLA_presence_thinned_sub_42$Longitude))
colnames(HOLA_presence_thinned_sub_42) <- c("longitude","latitude","HOLA_pa")

### Remove species column
HOLA_presence_above_42 <- HOLA_presence_above_42[,1:3]

### Combine presences and absences for HOLA
HOLA_presence_thinned_ss <- bind_rows(HOLA_presence_thinned_sub_42,HOLA_presence_above_42)

### Combine both occurrence datasets and the thinned absences
HOLA_p_a_data_ss <- bind_rows(HOLA_presence_thinned_ss, HOLA_absences_thin_df)


################
### Combine thinned presence and absence data for each species
################

BOBO_p_a_data <- bind_rows(BOBO_presence_thinned, BOBO_absences_thin_df)
EAME_p_a_data <- bind_rows(EAME_presence_thinned, EAME_absences_thin_df)
HOLA_p_a_data <- HOLA_p_a_data_ss # from above
SAVS_p_a_data <- bind_rows(SAVS_presence_thinned, SAVS_absences_thin_df)
KILL_p_a_data <- bind_rows(KILL_presence_thinned, KILL_absences_thin_df)
WOTH_p_a_data <- bind_rows(WOTH_presence_thinned, WOTH_absences_thin_df)
EWPE_p_a_data <- bind_rows(EWPE_presence_thinned, EWPE_absences_thin_df)
LEFL_p_a_data <- bind_rows(LEFL_presence_thinned, LEFL_absences_thin_df)


### Save thinned data for future use

write.csv(BOBO_p_a_data, "BOBO_p_a_data.csv")
write.csv(EAME_p_a_data, "EAME_p_a_data.csv")
write.csv(HOLA_p_a_data, "HOLA_p_a_data.csv")
write.csv(SAVS_p_a_data, "SAVS_p_a_data.csv")
write.csv(KILL_p_a_data, "KILL_p_a_data.csv")
write.csv(WOTH_p_a_data, "WOTH_p_a_data.csv")
write.csv(EWPE_p_a_data, "EWPE_p_a_data.csv")
write.csv(LEFL_p_a_data, "LEFL_p_a_data.csv")


############################################################################
### Convert observed presences and absences (2018 & 5yr) into spatial points 
### to eventually test accuracy of model predictions.
############################################################################

### Set new projection to align with raster outputs
geo_proj = "+proj=longlat +ellps=GRS80 +no_defs" 

### Subset thinned data into presence and absence again

BOBO_presence_thinned <- subset(BOBO_p_a_data, BOBO_pa == 1)
EAME_presence_thinned <- subset(EAME_p_a_data, EAME_pa == 1)
HOLA_presence_thinned <- subset(HOLA_p_a_data, HOLA_pa == 1)
SAVS_presence_thinned <- subset(SAVS_p_a_data, SAVS_pa == 1)
KILL_presence_thinned <- subset(KILL_p_a_data, KILL_pa == 1)
WOTH_presence_thinned <- subset(WOTH_p_a_data, WOTH_pa == 1)
EWPE_presence_thinned <- subset(EWPE_p_a_data, EWPE_pa == 1)
LEFL_presence_thinned <- subset(LEFL_p_a_data, LEFL_pa == 1)

BOBO_absences_thinned <- subset(BOBO_p_a_data, BOBO_pa == 0)
EAME_absences_thinned <- subset(EAME_p_a_data, EAME_pa == 0)
HOLA_absences_thinned <- subset(HOLA_p_a_data, HOLA_pa == 0)
SAVS_absences_thinned <- subset(SAVS_p_a_data, SAVS_pa == 0)
KILL_absences_thinned <- subset(KILL_p_a_data, KILL_pa == 0)
WOTH_absences_thinned <- subset(WOTH_p_a_data, WOTH_pa == 0)
EWPE_absences_thinned <- subset(EWPE_p_a_data, EWPE_pa == 0)
LEFL_absences_thinned <- subset(LEFL_p_a_data, LEFL_pa == 0)

##############################
### A) Occurrence observations
##############################

### Extract coordinates

BOBO_TP_coord <- BOBO_presence_thinned[,1:2] 
EAME_TP_coord <- EAME_presence_thinned[,1:2] 
HOLA_TP_coord <- HOLA_presence_thinned[,1:2] 
SAVS_TP_coord <- SAVS_presence_thinned[,1:2] 
KILL_TP_coord <- KILL_presence_thinned[,1:2] 
WOTH_TP_coord <- WOTH_presence_thinned[,1:2] 
EWPE_TP_coord <- EWPE_presence_thinned[,1:2] 
LEFL_TP_coord <- LEFL_presence_thinned[,1:2] 

### Convert to spatial points data frame.

BOBO_TP <- SpatialPointsDataFrame(coords = BOBO_TP_coord, data = BOBO_presence_thinned,
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

EAME_TP <- SpatialPointsDataFrame(coords = EAME_TP_coord, data = EAME_presence_thinned,
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

HOLA_TP <- SpatialPointsDataFrame(coords = HOLA_TP_coord, data = HOLA_presence_thinned,
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

SAVS_TP <- SpatialPointsDataFrame(coords = SAVS_TP_coord, data = SAVS_presence_thinned,
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

KILL_TP <- SpatialPointsDataFrame(coords = KILL_TP_coord, data = KILL_presence_thinned,
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

WOTH_TP <- SpatialPointsDataFrame(coords = WOTH_TP_coord, data = WOTH_presence_thinned,
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

EWPE_TP <- SpatialPointsDataFrame(coords = EWPE_TP_coord, data = EWPE_presence_thinned,
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

LEFL_TP <- SpatialPointsDataFrame(coords = LEFL_TP_coord, data = LEFL_presence_thinned,
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))


### Transform projections

BOBO_TP_use <- spTransform(BOBO_TP, geo_proj)
EAME_TP_use <- spTransform(EAME_TP, geo_proj)
HOLA_TP_use <- spTransform(HOLA_TP, geo_proj)
SAVS_TP_use <- spTransform(SAVS_TP, geo_proj)
KILL_TP_use <- spTransform(KILL_TP, geo_proj)
WOTH_TP_use <- spTransform(WOTH_TP, geo_proj)
EWPE_TP_use <- spTransform(EWPE_TP, geo_proj)
LEFL_TP_use <- spTransform(LEFL_TP, geo_proj)



##################################
### B) Absences for 2018
##################################

### Extract coordinates

BOBO_TA_coord_2018 <- BOBO_absences_thinned[,1:2] 
EAME_TA_coord_2018 <- EAME_absences_thinned[,1:2] 
HOLA_TA_coord_2018 <- HOLA_absences_thinned[,1:2] 
SAVS_TA_coord_2018 <- SAVS_absences_thinned[,1:2] 
KILL_TA_coord_2018 <- KILL_absences_thinned[,1:2] 
WOTH_TA_coord_2018 <- WOTH_absences_thinned[,1:2] 
EWPE_TA_coord_2018 <- EWPE_absences_thinned[,1:2] 
LEFL_TA_coord_2018 <- LEFL_absences_thinned[,1:2] 


### Convert to spatial points data frame

BOBO_TA_2018 <- SpatialPointsDataFrame(coords = BOBO_TA_coord_2018, data = BOBO_absences_thinned,
                                       proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

EAME_TA_2018 <- SpatialPointsDataFrame(coords = EAME_TA_coord_2018, data = EAME_absences_thinned,
                                       proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

HOLA_TA_2018 <- SpatialPointsDataFrame(coords = HOLA_TA_coord_2018, data = HOLA_absences_thinned,
                                       proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

SAVS_TA_2018 <- SpatialPointsDataFrame(coords = SAVS_TA_coord_2018, data = SAVS_absences_thinned,
                                       proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

KILL_TA_2018 <- SpatialPointsDataFrame(coords = KILL_TA_coord_2018, data = KILL_absences_thinned,
                                       proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

WOTH_TA_2018 <- SpatialPointsDataFrame(coords = WOTH_TA_coord_2018, data = WOTH_absences_thinned,
                                       proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

EWPE_TA_2018 <- SpatialPointsDataFrame(coords = EWPE_TA_coord_2018, data = EWPE_absences_thinned,
                                       proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

LEFL_TA_2018 <- SpatialPointsDataFrame(coords = LEFL_TA_coord_2018, data = LEFL_absences_thinned,
                                       proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))


### Transform projections

BOBO_TA_2018_use <- spTransform(BOBO_TA_2018, geo_proj)
EAME_TA_2018_use <- spTransform(EAME_TA_2018, geo_proj)
HOLA_TA_2018_use <- spTransform(HOLA_TA_2018, geo_proj)
SAVS_TA_2018_use <- spTransform(SAVS_TA_2018, geo_proj)
KILL_TA_2018_use <- spTransform(KILL_TA_2018, geo_proj)
WOTH_TA_2018_use <- spTransform(WOTH_TA_2018, geo_proj)
EWPE_TA_2018_use <- spTransform(EWPE_TA_2018, geo_proj)
LEFL_TA_2018_use <- spTransform(LEFL_TA_2018, geo_proj)


############################
### C) Absences over 5 years
############################

### Sum observations across years for each point to get an idea of presence
### over time (i.e., not just within 3 min)

bird_data_agg <- ddply(bird_data, .(route_ID, stop), summarise, 
                       longitude = mean(longitude),
                       latitude = mean(latitude),
                       BOBO = sum(BOBO),
                       EAME = sum(EAME),
                       HOLA = sum(HOLA),
                       SAVS = sum(SAVS),
                       KILL = sum(KILL),
                       WOTH = sum(WOTH),
                       EWPE = sum(EWPE),
                       LEFL = sum(LEFL))

### Extract absences
BOBO_absence_data <- subset(bird_data_agg, BOBO == 0)
EAME_absence_data <- subset(bird_data_agg, EAME == 0)
HOLA_absence_data <- subset(bird_data_agg, HOLA == 0)
SAVS_absence_data <- subset(bird_data_agg, SAVS == 0)
KILL_absence_data <- subset(bird_data_agg, KILL == 0)
WOTH_absence_data <- subset(bird_data_agg, WOTH == 0)
EWPE_absence_data <- subset(bird_data_agg, EWPE == 0)
LEFL_absence_data <- subset(bird_data_agg, LEFL == 0)


### Extract lat and long for each absence
BOBO_coord <- BOBO_absence_data[,3:4] 
EAME_coord <- EAME_absence_data[,3:4]
HOLA_coord <- HOLA_absence_data[,3:4]
SAVS_coord <- SAVS_absence_data[,3:4]
KILL_coord <- KILL_absence_data[,3:4]
WOTH_coord <- WOTH_absence_data[,3:4]
EWPE_coord <- EWPE_absence_data[,3:4]
LEFL_coord <- LEFL_absence_data[,3:4]


### Convert to spatialpointsdataframe
BOBO_TA <- SpatialPointsDataFrame(coords = BOBO_coord, data = BOBO_absence_data,
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

EAME_TA <- SpatialPointsDataFrame(coords = EAME_coord, data = EAME_absence_data,
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

HOLA_TA <- SpatialPointsDataFrame(coords = HOLA_coord, data = HOLA_absence_data,
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

SAVS_TA <- SpatialPointsDataFrame(coords = SAVS_coord, data = SAVS_absence_data,
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

KILL_TA <- SpatialPointsDataFrame(coords = KILL_coord, data = KILL_absence_data,
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

WOTH_TA <- SpatialPointsDataFrame(coords = WOTH_coord, data = WOTH_absence_data,
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

EWPE_TA <- SpatialPointsDataFrame(coords = EWPE_coord, data = EWPE_absence_data,
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

LEFL_TA <- SpatialPointsDataFrame(coords = LEFL_coord, data = LEFL_absence_data,
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))


### Transform
BOBO_TA_use <- spTransform(BOBO_TA, geo_proj)
EAME_TA_use <- spTransform(EAME_TA, geo_proj)
HOLA_TA_use <- spTransform(HOLA_TA, geo_proj)
SAVS_TA_use <- spTransform(SAVS_TA, geo_proj)
KILL_TA_use <- spTransform(KILL_TA, geo_proj)
WOTH_TA_use <- spTransform(WOTH_TA, geo_proj)
EWPE_TA_use <- spTransform(EWPE_TA, geo_proj)
LEFL_TA_use <- spTransform(LEFL_TA, geo_proj)


################################################################################
################# Pre-processing land cover data  ##############################
################################################################################

### Read in raster stacks (already pre-processed)

landcover_2018_ss <- stack("landcover_2018_study_site_final.tif")
landcover_2018_sf_ss <- stack("landcover_2018_sf_study_site_final.tif")

### Rename levels

names(landcover_2018_ss) <- c("water_200m","urban_200m","shrub_200m","wetland_200m",
                              "corn_200m","soy_200m","forest_200m","grassland_pasture_200m",
                              "cereals_200m","lp_veg_200m","fruit_total_200m",
                              "water_1km","urban_1km","shrub_1km","wetland_1km",
                              "corn_1km","soy_1km","forest_1km","grassland_pasture_1km",
                              "cereals_1km","lp_veg_1km","fruit_total_1km","elevation",
                              "slope","aspect","bios_10","bios_18")

names(landcover_2018_sf_ss) <- c("water_200m","urban_200m","shrub_200m",
                                 "wetland_200m","corn_200m","soy_200m",
                                 "coniferous_200m","deciduous_200m",
                                 "mixed_forest_200m","grassland_pasture_200m",
                                 "cereals_200m","lp_veg_200m","fruit_total_200m",
                                 "water_1km","urban_1km","shrub_1km","wetland_1km",
                                 "corn_1km","soy_1km","coniferous_1km",
                                 "deciduous_1km","mixed_forest_1km",
                                 "grassland_pasture_1km","cereals_1km","lp_veg_1km",
                                 "fruit_total_1km","elevation","slope","aspect",
                                 "bios_10","bios_18")

################################################
### Finalize land cover data for grassland species
################################################

### Remove urban 200m as it focuses the predictions too much on the roads.
### Remove fruits and veg as only minor influence.
### Keep forest as total cover and separate out main crop types.

landcover_multiscale_ss <- subset(landcover_2018_ss, 
                                  c("water_200m","shrub_200m",
                                    "wetland_200m","corn_200m","soy_200m",
                                    "forest_200m","grassland_pasture_200m",
                                    "cereals_200m","water_1km","urban_1km",
                                    "shrub_1km","wetland_1km","corn_1km",
                                    "soy_1km","forest_1km","grassland_pasture_1km",
                                    "cereals_1km","elevation","slope","aspect",
                                    "bios_10","bios_18"))


#################################################
### Create land cover data for forest specialists
#################################################

### Collapse 3 main crop types into a single 'crops' layer

landcover_2018_sf_ss$crops_200m <- sum(landcover_2018_sf_ss$soy_200m,
                                       landcover_2018_sf_ss$corn_200m,
                                       landcover_2018_sf_ss$cereals_200m)

landcover_2018_sf_ss$crops_1km <- sum(landcover_2018_sf_ss$soy_1km,
                                      landcover_2018_sf_ss$corn_1km,
                                      landcover_2018_sf_ss$cereals_1km)


### Remove urban 200m as it focuses the predictions too much on the roads.
### Also remove lp veg, fruits and individual crops

landcover_2018_forest_sf_final <- subset(landcover_2018_sf_ss, 
                                         c("water_200m","shrub_200m",
                                           "crops_200m","wetland_200m","deciduous_200m",
                                           "coniferous_200m", "mixed_forest_200m",
                                           "grassland_pasture_200m","water_1km","shrub_1km",
                                           "wetland_1km","deciduous_1km","coniferous_1km",
                                           "mixed_forest_1km","grassland_pasture_1km",
                                           "crops_1km","elevation","slope","aspect","bios_10",
                                           "bios_18","urban_1km"))


################################################################################
############### Fit presence-absence ensemble models ###########################
################################################################################


###########
### A) BOBO
###########

BOBO_ta_final_ses <- ensemble_modelling(c("GLM","GAM","RF"), 
                            cta.args=list(ntree=5000),
                            BOBO_p_a_data, landcover_multiscale_ss,
                            Pcol = "BOBO_pa", Xcol = "longitude", Ycol = "latitude",
                            rep = 10, uncertainty = TRUE, bin.thresh="SES",
                            cv = "holdout", cv.param = c(0.70,5),
                            final.fit.data = "holdout",
                            ensemble.metric = "AUC", ensemble.thresh = 0.70,
                            weight = TRUE, verbose = TRUE)

### Save output
saveRDS(BOBO_ta_final_ses, "BOBO_ta_final_ses.rds")

###########
### B) EAME
###########

EAME_ta_final_ses <- ensemble_modelling(c("GLM","GAM","RF"), 
                         cta.args=list(ntree=5000),
                         EAME_p_a_data, landcover_multiscale_ss,
                         Pcol = "EAME_pa", Xcol = "longitude", Ycol = "latitude",
                         rep = 10, uncertainty = TRUE, bin.thresh="SES",
                         cv = "holdout", cv.param = c(0.7,5),
                         final.fit.data = "holdout",
                         ensemble.metric = "AUC", ensemble.thresh = 0.50,
                         weight = TRUE, verbose = TRUE)

### Save output
saveRDS(EAME_ta_final_ses, "EAME_ta_final_ses.rds")


###########
### C) HOLA
###########

HOLA_ta_final_ses <- ensemble_modelling(c("GLM","GAM","RF"), 
                         cta.args=list(ntree=5000),
                         HOLA_p_a_data, landcover_multiscale_ss,
                         Pcol = "HOLA_pa", Xcol = "longitude", Ycol = "latitude",
                         rep = 10, uncertainty = TRUE, bin.thresh="SES",
                         cv = "holdout", cv.param = c(0.7,5),
                         final.fit.data = "holdout",
                         ensemble.metric = "AUC", ensemble.thresh = 0.50,
                         weight = TRUE, verbose = TRUE)

### Save output
saveRDS(HOLA_ta_final_ses, "HOLA_ta_final_ses.rds")


###########
### D) SAVS
###########

SAVS_ta_final_ses <- ensemble_modelling(c("GLM","GAM","RF"), 
                         cta.args=list(ntree=5000),
                         SAVS_p_a_data, landcover_multiscale_ss,
                         Pcol = "SAVS_pa", Xcol = "longitude", Ycol = "latitude",
                         rep = 10, uncertainty = TRUE, bin.thresh="SES",
                         cv = "holdout", cv.param = c(0.7,5),
                         final.fit.data = "holdout",
                         ensemble.metric = "AUC", ensemble.thresh = 0.40,
                         weight = TRUE, verbose = TRUE)

### Save output
saveRDS(SAVS_ta_final_ses, "SAVS_ta_final_ses.rds")



###########
### E) KILL
###########

KILL_ta_final_ses <- ensemble_modelling(c("GLM","GAM","RF"), 
                        cta.args=list(ntree=5000),
                        KILL_p_a_data, landcover_multiscale_ss,
                        Pcol = "KILL_pa", Xcol = "longitude", Ycol = "latitude",
                        rep = 10, uncertainty = TRUE, bin.thresh="SES",
                        cv = "holdout", cv.param = c(0.7,5),
                        final.fit.data = "holdout",
                        ensemble.metric = "AUC", ensemble.thresh = 0.50,
                        weight = TRUE, verbose = TRUE)

### Save output
saveRDS(KILL_ta_final_ses, "KILL_ta_final_ses.rds")



###########
### F) WOTH
###########

WOTH_ta_final_ses <- ensemble_modelling(c("GLM","GAM","RF"), 
                        cta.args=list(ntree=5000),
                        WOTH_p_a_data, landcover_2018_forest_sf_final,
                        Pcol = "WOTH_pa", Xcol = "longitude", Ycol = "latitude",
                        rep = 10, uncertainty = TRUE, bin.thresh="SES",
                        cv = "holdout", cv.param = c(0.7,5),
                        final.fit.data = "holdout",
                        ensemble.metric = "AUC", ensemble.thresh = 0.5,
                        weight = TRUE, verbose = TRUE)

### Save output
saveRDS(WOTH_ta_final_ses, "WOTH_ta_final_ses.rds")



###########
### G) EWPE
###########

EWPE_ta_final_ses <- ensemble_modelling(c("GLM","GAM","RF"), 
                        cta.args=list(ntree=5000),
                        EWPE_p_a_data, landcover_2018_forest_sf_final,
                        Pcol = "EWPE_pa", Xcol = "longitude", Ycol = "latitude",
                        rep = 10, uncertainty = TRUE, bin.thresh="SES",
                        cv = "holdout", cv.param = c(0.7,5),
                        final.fit.data = "holdout",
                        ensemble.metric = "AUC", ensemble.thresh = 0.5,
                        weight = TRUE, verbose = TRUE)

### Save output
saveRDS(EWPE_ta_final_ses, "EWPE_ta_final_ses.rds")


###########
### H) LEFL
###########

LEFL_ta_final_ses <- ensemble_modelling(c("GLM","GAM","RF"), 
                        cta.args=list(ntree=5000),
                        LEFL_p_a_data, landcover_2018_forest_sf_final,
                        Pcol = "LEFL_pa", Xcol = "longitude", Ycol = "latitude",
                        rep = 10, uncertainty = TRUE, bin.thresh="SES",
                        cv = "holdout", cv.param = c(0.7,5),
                        final.fit.data = "holdout",
                        ensemble.metric = "AUC", ensemble.thresh = 0.5,
                        weight = TRUE, verbose = TRUE)

### Save output
saveRDS(LEFL_ta_final_ses, "LEFL_ta_final_ses.rds")




#########################################
### Convert projections to binary rasters
#########################################

### BOBO

BOBO_ta_final_ses@algorithm.evaluation

BOBO_projection <- BOBO_ta_final_ses@projection
BOBO_reclassified <- reclassify(BOBO_projection, 
                                cbind(0.187, Inf, 1))
BOBO_reclassified_final <- reclassify(BOBO_reclassified, 
                                      cbind(-Inf, 0.187, 0))
writeRaster(BOBO_reclassified_final, "BOBO_ta_reclassified_final.tif", overwrite=T)


### EAME
EAME_ta_final_ses@algorithm.evaluation

EAME_projection <- EAME_ta_final_ses@projection
EAME_reclassified <- reclassify(EAME_projection, 
                                cbind(0.165, Inf, 1))
EAME_reclassified_final <- reclassify(EAME_reclassified, 
                                      cbind(-Inf, 0.165, 0))
writeRaster(EAME_reclassified_final, "EAME_ta_reclassified_final.tif", overwrite=T)


### HOLA

HOLA_ta_final_ses@algorithm.evaluation

HOLA_projection <- HOLA_ta_final_ses@projection
HOLA_reclassified <- reclassify(HOLA_projection, 
                                    cbind(0.13, Inf, 1))
HOLA_reclassified_final <- reclassify(HOLA_reclassified, 
                                          cbind(-Inf, 0.13, 0))
writeRaster(HOLA_reclassified_final, "HOLA_ta_reclassified_final.tif", overwrite=T)


### SAVS

SAVS_ta_final_ses@algorithm.evaluation

SAVS_projection <- SAVS_ta_final_ses@projection
SAVS_reclassified <- reclassify(SAVS_projection, 
                                cbind(0.29, Inf, 1))
SAVS_reclassified_final <- reclassify(SAVS_reclassified, 
                                      cbind(-Inf, 0.29, 0))
writeRaster(SAVS_reclassified_final, "SAVS_ta_reclassified_final.tif", overwrite=T)


### KILL

KILL_ta_final_ses@algorithm.evaluation

KILL_projection <- KILL_ta_final_ses@projection
KILL_reclassified <- reclassify(KILL_projection, 
                                cbind(0.18, Inf, 1))
KILL_reclassified_final <- reclassify(KILL_reclassified, 
                                      cbind(-Inf, 0.18, 0))
writeRaster(KILL_reclassified_final, "KILL_ta_reclassified_final.tif", overwrite=T)


### WOTH

WOTH_ta_final_ses@algorithm.evaluation

WOTH_projection <- WOTH_ta_final_ses@projection
WOTH_reclassified_new <- reclassify(WOTH_projection, 
                                    cbind(0.11, Inf, 1))
WOTH_reclassified_final_new <- reclassify(WOTH_reclassified_new, 
                                          cbind(-Inf, 0.11, 0))
writeRaster(WOTH_reclassified_final_new, "WOTH_reclassified_final.tif", overwrite=T)


### EWPE

EWPE_ta_final_ses@algorithm.evaluation

EWPE_projection <- EWPE_ta_final_ses@projection
EWPE_reclassified_new <- reclassify(EWPE_projection, 
                                    cbind(0.105, Inf, 1))
EWPE_reclassified_final_new <- reclassify(EWPE_reclassified_new, 
                                          cbind(-Inf, 0.105, 0))
writeRaster(EWPE_reclassified_final_new, "EWPE_ta_reclassified_final.tif", overwrite=T)


### LEFL

LEFL_ta_final_ses@algorithm.evaluation

LEFL_projection <- LEFL_ta_final_ses@projection
LEFL_reclassified_new <- reclassify(LEFL_projection, 
                                    cbind(0.095, Inf, 1))
LEFL_reclassified_final_new <- reclassify(LEFL_reclassified_new, 
                                          cbind(-Inf, 0.095, 0))

### Extra step is to mask water

water <- landcover_multiscale_ss$water_200m

LEFL_reclassified_final_wow <- mask(LEFL_reclassified_final_new, water, 
                                    maskvalue=1)

writeRaster(LEFL_reclassified_final_wow, "LEFL_reclassified_final.tif", overwrite=T)

plot(LEFL_reclassified_final_wow)




##################################
### Final individual species plots
##################################

### Read in study site shapefile for border
study_area <- readOGR(dsn = "EasternBoundary", "Eastern")

### Reproject study area shapefile
study_proj <- "+proj=longlat +ellps=GRS80 +no_defs"
study_area_new <- spTransform(study_area, study_proj)


### BOBO map

png(file="BOBO_distribution.png",width=3500,height=2500, res=600)
pal <- colorRampPalette(c("transparent","#993404"))
plot(BOBO_reclassified_final, col = pal(2),
     xlab = "Longitude", ylab = "Latitude", cex.lab=1.5, xaxt='n', yaxt='n')
axis(side = 2, at=c(42,44,46,48), las=1, cex.axis=1.2)
axis(side = 1, at=c(-82,-78,-74,-70), cex.axis=1.2)
plot(study_area_new, add=T, border = "black", lwd=0.3)
dev.off()

### EAME map

png(file="EAME_distribution.png",width=3500,height=2500, res=600)
pal <- colorRampPalette(c("transparent","#993404"))
plot(EAME_reclassified_final, col = pal(2),
     xlab = "Longitude", ylab = "Latitude", cex.lab=1.5, xaxt='n', yaxt='n')
axis(side = 2, at=c(42,44,46,48), las=1, cex.axis=1.2)
axis(side = 1, at=c(-82,-78,-74,-70), cex.axis=1.2)
plot(study_area_new, add=T, border = "black", lwd=0.3)
dev.off()

### HOLA map

png(file="HOLA_distribution.png",width=3500,height=2500, res=600)
pal <- colorRampPalette(c("transparent","#993404"))
plot(HOLA_reclassified_final, col = pal(2),
     xlab = "Longitude", ylab = "Latitude", cex.lab=1.5, xaxt='n', yaxt='n')
axis(side = 2, at=c(42,44,46,48), las=1, cex.axis=1.2)
axis(side = 1, at=c(-82,-78,-74,-70), cex.axis=1.2)
plot(study_area_new, add=T, border = "black", lwd=0.3)
dev.off()


### KILL map

png(file="KILL_distribution.png",width=3500,height=2500, res=600)
pal <- colorRampPalette(c("transparent","#993404"))
plot(KILL_reclassified_final, col = pal(2),
     xlab = "Longitude", ylab = "Latitude", cex.lab=1.5, xaxt='n', yaxt='n')
axis(side = 2, at=c(42,44,46,48), las=1, cex.axis=1.2)
axis(side = 1, at=c(-82,-78,-74,-70), cex.axis=1.2)
plot(study_area_new, add=T, border = "black", lwd=0.3)
dev.off()


### SAVS map

png(file="SAVS_distribution.png",width=3500,height=2500, res=600)
pal <- colorRampPalette(c("transparent","#993404"))
plot(SAVS_reclassified_final, col = pal(2),
     xlab = "Longitude", ylab = "Latitude", cex.lab=1.5, xaxt='n', yaxt='n')
axis(side = 2, at=c(42,44,46,48), las=1, cex.axis=1.2)
axis(side = 1, at=c(-82,-78,-74,-70), cex.axis=1.2)
plot(study_area_new, add=T, border = "black", lwd=0.3)
dev.off()



### WOTH map

png(file="WOTH_distribution.png",width=3500,height=2500, res=600)
pal <- colorRampPalette(c("transparent","#35978f"))
plot(WOTH_reclassified_final, col = pal(2),
     xlab = "Longitude", ylab = "Latitude", cex.lab=1.5, xaxt='n', yaxt='n')
axis(side = 2, at=c(42,44,46,48), las=1, cex.axis=1.2)
axis(side = 1, at=c(-82,-78,-74,-70), cex.axis=1.2)
plot(study_area_new, add=T, border = "black", lwd=0.3)
dev.off()


### EWPE map

png(file="EWPE_distribution.png",width=3500,height=2500, res=600)
pal <- colorRampPalette(c("transparent","#35978f"))
plot(EWPE_reclassified_final, col = pal(2),
     xlab = "Longitude", ylab = "Latitude", cex.lab=1.5, xaxt='n', yaxt='n')
axis(side = 2, at=c(42,44,46,48), las=1, cex.axis=1.2)
axis(side = 1, at=c(-82,-78,-74,-70), cex.axis=1.2)
plot(study_area_new, add=T, border = "black", lwd=0.3)
dev.off()


### LEFL map

png(file="LEFL_distribution.png",width=3500,height=2500, res=600)
pal <- colorRampPalette(c("transparent","#35978f"))
plot(LEFL_reclassified_final_wow, col = pal(2),
     xlab = "Longitude", ylab = "Latitude", cex.lab=1.5, xaxt='n', yaxt='n')
axis(side = 2, at=c(42,44,46,48), las=1, cex.axis=1.2)
axis(side = 1, at=c(-82,-78,-74,-70), cex.axis=1.2)
plot(study_area_new, add=T, border = "black", lwd=0.3)
dev.off()



################################
### Ensemble model correlations
################################

BOBO_ta_final_ses@algorithm.correlation
EAME_ta_final_ses@algorithm.correlation
HOLA_ta_final_ses@algorithm.correlation
SAVS_ta_final_ses@algorithm.correlation
KILL_ta_final_ses@algorithm.correlation
WOTH_ta_final_ses@algorithm.correlation
EWPE_ta_final_ses@algorithm.correlation
LEFL_ta_final_ses@algorithm.correlation


#########################
### Variable importance
#########################

BOBO_ta_final_ses@variable.importance
EAME_ta_final_ses@variable.importance
HOLA_ta_final_ses@variable.importance
SAVS_ta_final_ses@variable.importance
KILL_ta_final_ses@variable.importance
WOTH_ta_final_ses@variable.importance
EWPE_ta_final_ses@variable.importance
LEFL_ta_final_ses@variable.importance



################################################################################
####### Species richness - stacked species distribution models (SSDMs) #########
################################################################################


### Read in final binary outputs for each species

BOBO_reclassified_final <- raster("BOBO_ta_reclassified_final.tif")
EAME_reclassified_final <- raster("EAME_ta_reclassified_final.tif")
HOLA_reclassified_final <- raster("HOLA_ta_reclassified_final.tif")
SAVS_reclassified_final <- raster("SAVS_ta_reclassified_final.tif")
KILL_reclassified_final <- raster("KILL_ta_reclassified_final.tif")
WOTH_reclassified_final <- raster("WOTH_ta_reclassified_final.tif")
EWPE_reclassified_final <- raster("EWPE_ta_reclassified_final.tif")
LEFL_reclassified_final <- raster("LEFL_ta_reclassified_final.tif")


#########################
### A) Grassland species
#########################

### Rename so do not overwrite
ESDM1 <- BOBO_ta_final_ses
ESDM2 <- EAME_ta_final_ses
ESDM3 <- HOLA_ta_final_ses
ESDM4 <- SAVS_ta_final_ses
ESDM5 <- KILL_ta_final_ses

### Name slot so code sees each model as different
ESDM1@name <- "BOBO_ESDM" 
ESDM2@name <- "EAME_ESDM"
ESDM3@name <- "HOLA_ESDM"
ESDM4@name <- "SAVS_ESDM"
ESDM5@name <- "KILL_ESDM"

### Generate uncertainty and variable importance
oc_stack_bin <- stacking(ESDM1,ESDM2,ESDM3,ESDM4,ESDM5,
                           method = "bSSDM", Env = landcover_multiscale_ss, 
                           eval = TRUE, verbose = TRUE)

### Stack binary predictions
oc_stack_opt <- stack(BOBO_reclassified_final,EAME_reclassified_final,
                      HOLA_reclassified_final,SAVS_reclassified_final,
                      KILL_reclassified_final)

### Sum across raster stack
oc_stack_sum <- calc(oc_stack_opt, sum)


### Remove water
water <- landcover_multiscale_ss$water_200m
oc_stack_final <- mask(oc_stack_sum, water, maskvalue= 1)
plot(oc_stack_final)


### Produce figure of grassland specialist richness
png(file="oc_richness_final.png",width=3500,height=2500, res=600)
cuts <- c(-1,0,1,2,3,4,5)
pal <- colorRampPalette(viridis(6, option="magma"))
plot(oc_stack_final, breaks=cuts,  col = pal(6),
     xlab = "Longitude", ylab = "Latitude", cex.lab=1.5, xaxt='n', yaxt='n')
axis(side = 2, at=c(42,44,46,48), las=1, cex.axis=1.2)
axis(side = 1, at=c(-82,-78,-74,-70), cex.axis=1.2)
plot(study_area_new, add=T, border = "black", lwd=0.3)
dev.off()


####################################
### Bar plot for variable importance
####################################

oc_stack_bin@variable.importance

### Create vectors for variables, influence, and sd
variables <- c("water 200m","shrub 200m","wetland 200m","corn 200m","soy 200m",
               "forest 200m","grassland 200m","cereals 200m","water 1km",
               "urban 1km","shrub 1km","wetland 1km","corn 1km","soy 1km",
               "forest 1km","grassland 1km","cereals 1km","elevation",
               "slope","aspect","bios 10","bios 18","crops 200m","crops 1km",
               "deciduous 200m","deciduous 1km","coniferous 200m","coniferous 1km",
               "mixed forest 200m","mixed forest 1km")

influence <- c(1.06,1.38,1.52,5.36,9.82,5.12,25.3,5.90,1.34,1.40,2.60,2.70,
               4.11,5.23,2.25,3.25,5.70,5.46,1.25,1.38,4.93,2.91,0,0,0,0,0,0,0,0)

st_dev <- c(0.70,0.90,0.68,3.63,5.67,3.27,21.15,3.82,0.37,1.09,1.99,2.32,
            2.66,3.70,1.27,2.07,3.69,4.94,0.20,0.39,3.07,1.75,0,0,0,0,0,0,0,0)

### Combine into dataframe
oc_influence <- data.frame(variables,influence,st_dev)


### Make variables into a factor with the correct order
oc_influence$variables_f <- factor(oc_influence$variables, 
                                   levels=c("water 200m","shrub 200m","wetland 200m","grassland 200m",
                                            "crops 200m","cereals 200m","corn 200m","soy 200m",
                                            "forest 200m","deciduous 200m","coniferous 200m",
                                            "mixed forest 200m","water 1km","urban 1km","shrub 1km",
                                            "wetland 1km","grassland 1km","crops 1km","cereals 1km",
                                            "corn 1km","soy 1km","forest 1km","deciduous 1km",
                                            "coniferous 1km","mixed forest 1km","elevation","slope",
                                            "aspect","bios 10","bios 18"))


head(oc_influence)


### Plot variable importance
png(file="oc_richness_barplot.png",width=3500,height=2500, res=600)
oc_barplot <- ggplot(oc_influence) +
              geom_bar(aes(x=variables_f, y=influence), stat="identity", fill="#b2182b", alpha=0.9) +
              labs(x = "Predictors", y = "Variable contribution (%)", cex = 5) +
              geom_errorbar(aes(x=variables_f, ymin=influence-st_dev, ymax=influence+st_dev), 
                            width=0.4, colour="black", alpha=1, size=1.2) +
              theme_classic() +
                  theme(axis.text=element_text(size=12, colour="black"),
                        axis.title=element_text(size=16),
                        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
oc_barplot
dev.off()


##############################
### B) Forest specialists only
##############################

ESDM6_f <- WOTH_ta_final_ses
ESDM7_f <- EWPE_ta_final_ses
ESDM8_f <- LEFL_ta_final_ses

ESDM6_f@name <- "WOTH_ESDM" 
ESDM7_f@name <- "EWPE_ESDM"
ESDM8_f@name <- "LEFL_ESDM"

### Generate uncertainty and variable importance
forest_stack_bin <- stacking(ESDM6_f,ESDM7_f,ESDM8_f,
                               method = "bSSDM", Env = landcover_2018_crop_final_ss, 
                               eval = TRUE, verbose = TRUE)

#saveRDS(forest_stack_ml, "forest_stack.rds")


### Stack binary predictions

forest_stack_opt <- stack(WOTH_reclassified_final,EWPE_reclassified_final,
                          LEFL_reclassified_final)

### Sum across raster stack
forest_stack_sum <- calc(forest_stack_opt, sum)

### Remove water
water <- landcover_multiscale_ss$water_200m
forest_stack_final <- mask(forest_stack_sum, water, maskvalue= 1)

### Create forest specialist richness figure
png(file="forest_richness_final_magma.png",width=3500,height=2500, res=600)
cuts <- c(-1,0,1,2,3)
pal <- colorRampPalette(viridis(4, option="magma"))
plot(forest_stack_final, breaks=cuts, col = pal(4),
     xlab = "Longitude", ylab = "Latitude", cex.lab=1.5, xaxt='n', yaxt='n')
axis(side = 2, at=c(42,44,46,48), las=1, cex.axis=1.2)
axis(side = 1, at=c(-82,-78,-74,-70), cex.axis=1.2)
dev.off()


####################################
### Bar plot for variable importance
####################################

### Create vectors for variables, influence, and sd
variables <- c("water 200m","shrub 200m","crops 200m","wetland 200m","deciduous 200m",
               "coniferous 200m","mixed forest 200m","grassland 200m","water 1km",
               "shrub 1km","wetland 1km","deciduous 1km","coniferous 1km","mixed forest 1km",
               "grassland 1km","crops 1km","elevation","slope","aspect","bios 10","bios 18",
               "urban 1km","cereals 200m","cereals 1km","soy 200m","soy 1km",
               "corn 200m","corn 1km","forest 200m","forest 1km")

influence <- c(2.80,3.67,2.15,2.60,21.00,3.09,9.50,2.38,2.05,2.43,4.25,5.87,3.64,1.62,3.52,
               2.52,3.47,7.32,4.43,7.12,3.10,2.47,0,0,0,0,0,0,0,0)

st_dev <- c(1.47,2.78,1.39,1.75,12.61,2.45,6.78,1.62,1.33,0.71,2.52,6.10,1.75,
            0.51,2.47,1.63,0.94,10.16,2.74,3.69,0.37,1.18,0,0,0,0,0,0,0,0)

for_influence <- data.frame(variables,influence,st_dev)

### Make variables into a factor with the correct order
for_influence$variables_f <- factor(for_influence$variables, 
                                   levels=c("water 200m","shrub 200m","wetland 200m","grassland 200m",
                                            "crops 200m","cereals 200m","corn 200m","soy 200m",
                                            "forest 200m","deciduous 200m","coniferous 200m",
                                            "mixed forest 200m","water 1km","urban 1km","shrub 1km",
                                            "wetland 1km","grassland 1km","crops 1km","cereals 1km",
                                            "corn 1km","soy 1km","forest 1km","deciduous 1km",
                                            "coniferous 1km","mixed forest 1km","elevation","slope",
                                            "aspect","bios 10","bios 18"))


head(for_influence)

### Barplot of variable importance
png(file="for_richness_barplot.png",width=3500,height=2500, res=600)
for_barplot <- ggplot(for_influence) +
  geom_bar(aes(x=variables_f, y=influence), stat="identity", fill="#35978f", alpha=0.9) +
  labs(x = "Predictors", y = "Variable contribution (%)", cex = 5) +
  geom_errorbar(aes(x=variables_f, ymin=influence-st_dev, ymax=influence+st_dev), 
                width=0.4, colour="black", alpha=1, size=1.2) +
  theme_classic() +
  theme(axis.text=element_text(size=12, colour="black"),
        axis.title=element_text(size=16),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
for_barplot
dev.off()



################################
### C) Combined SAR/SOSC species
################################

### Aggregate previous layers to 2kmx2km (increase by factor of 3.7)

oc_agg <- aggregate(oc_stack_final, fact=3.7, fun=max)
forest_agg <- aggregate(forest_stack_final, fact=3.7, fun=max)

### Combine
sar_stack_agg <- stack(oc_agg,forest_agg)

### Sum across layers
sar_stack_agg_sum <- calc(sar_stack_agg, sum)

### Plot SAR/SOC species richness at 2km res 
png(file="sar_richness_2km.png",width=3500,height=2500, res=600)
cuts <- c(-1,0,1,2,3,4,5,6,7,8)
pal <- colorRampPalette(viridis(9, option="viridis"))
plot(sar_stack_agg_sum, breaks=cuts, col = pal(9),
     xlab = "Longitude", ylab = "Latitude", cex.lab=1.5, xaxt='n', yaxt='n')
axis(side = 2, at=c(42,44,46,48), las=1, cex.axis=1.2)
axis(side = 1, at=c(-82,-78,-74,-70), cex.axis=1.2)
dev.off()



####################################################
### Habitat associations for SAR/SOSC richness
####################################################

### Create dataset

landcover_sar_sosc <- subset(landcover_2018_sf_ss, 
                             c("water_200m","shrub_200m",
                               "wetland_200m","corn_200m","soy_200m",
                               "coniferous_200m","deciduous_200m",
                               "mixed_forest_200m","grassland_pasture_200m",
                               "cereals_200m",
                               "water_1km","urban_1km","shrub_1km","wetland_1km",
                               "corn_1km","soy_1km","coniferous_1km",
                               "deciduous_1km","mixed_forest_1km",
                               "grassland_pasture_1km","cereals_1km","elevation",
                               "slope","aspect",
                               "bios_10","bios_18"))

### Extract names
landcover_levels <- names(landcover_sar_sosc)


### Write a loop that creates a raster stack with habitat variables
landcover_sar_sosc_2km <- list()

for(i in 1:length(landcover_levels)) {
  
  # Extract proportion of landcover values and combine in list
  landcover_sar_sosc_2km[[i]] <- aggregate(landcover_sar_sosc[[i]], fact=3.7, fun=mean)
  
  # Create raster stacks from lists
  landcover_sar_sosc_2km_stack <- stack(landcover_sar_sosc_2km, RAT = TRUE)
  
   # Show progress
  print(i)
  
}

### Save each as a Geotiff

#writeRaster(landcover_sar_sosc_2km_stack, filename= "landcover_sar_sosc_2km_stack.tif", 
           #format="GTiff", options = "INTERLEAVE=BAND", overwrite=TRUE)


### Extract pixel values and make a data set
sar_sosc_data <- data.frame(values(sar_stack_agg_sum), values(landcover_sar_sosc_2km_stack))

### Change column name
colnames(sar_sosc_data)[1] <- c("sar_sosc_richness")

### Grab coordinates
sar_sosc_data$latitude <- coordinates(sar_stack_agg_sum)[,2]
sar_sosc_data$longitude <- coordinates(sar_stack_agg_sum)[,1]

### Remove Nas
sar_sosc_data_use <- sar_sosc_data[complete.cases(sar_sosc_data), ]

### Thin datapoints randomly to reduce computer time and avoid spatial
### autocorrelation

sar_sosc_data_reduced <- sample_n(sar_sosc_data_use, 20000)


### Fit GAM model

system.time(sar_richness_gam <- gam(sar_sosc_richness ~ 
                                    + s(wetland_200m,  k=3)
                                    + s(water_200m, k=3)
                                    + s(shrub_200m, k=3)
                                    + s(grassland_pasture_200m, k=3)
                                    + s(deciduous_200m, k=3)
                                    + s(coniferous_200m, k=3)
                                    + s(mixed_forest_200m, k=3)
                                    + s(soy_200m, k=3)
                                    + s(corn_200m, k=3)
                                    + s(cereals_200m,  k=3)
                                    + s(urban_1km, k=3)
                                    + s(wetland_1km, k=3)
                                    + s(water_1km, k=3)
                                    + s(shrub_1km,  k=3)
                                    + s(grassland_pasture_1km, k=3)
                                    + s(deciduous_1km, k=3)
                                    + s(coniferous_1km, k=3)
                                    + s(mixed_forest_1km, k=3)
                                    + s(soy_1km, k=3)
                                    + s(corn_1km, k=3)
                                    + s(cereals_1km, k=3)
                                    + s(elevation, k=3)
                                    + s(aspect, k=3)
                                    + s(slope, k=3)
                                    + s(bios_10, k=3)
                                    + s(bios_18, k=3)
                                    , family=gaussian, data= sar_sosc_data_reduced, 
                                    select=TRUE, method = 'REML',
                                    control= list(nthreads=4)))

summary(sar_richness_gam)


### Plot and save the top 9 habitat variables
theme_set(theme_classic())

png(file="richness_gam_9.png",width=3500,height=2500, res=600)
draw(sar_richness_gam, select = c(3,4,5,12,15,20,21,25,26))
dev.off()



#############################################
### D) Species richness including all species
#############################################

### Read in bird data
bird_data <- read.csv("bird_landcover_weather_terrain_data_2014_2018.csv")

head(bird_data)

### Subset to 2018 only
bird_data_2018 <- subset(bird_data, year==2018)


######################################################
### Stacked species distribution model (SSDM) approach
######################################################

# Read in candidate bird list (pre-filtered)
bird_list <- read.csv("Eastern farmland bird community.csv")

### Select data of interest only
bird_data_2018_use <- bird_data_2018 %>% 
                      select(longitude,latitude,bird_list$species_code)

### Convert abundance into presence absence
bird_data_2018_use[3:156] <- lapply(bird_data_2018_use[3:156], 
                                    function(x) ifelse(x > 0,1,0))

### Separate into forest and open country specialists

oc_list <- subset(bird_list, habitat=="Open")
forest_list <- subset(bird_list, habitat=="Forest")

oc_bird_data <- bird_data_2018_use %>% select(longitude,latitude,oc_list$species_code)
forest_bird_data <- bird_data_2018_use %>% select(longitude,latitude,forest_list$species_code)


### Subset only species with greater than 50 individuals
oc_birds_keep <- which(colSums(oc_bird_data[3:48]) >= 50)
forest_birds_keep <- which(colSums(forest_bird_data[3:110]) >= 50)

oc_birds_use <- oc_bird_data %>% select(longitude,latitude,names(oc_birds_keep))
forest_birds_use <- forest_bird_data %>% select(longitude,latitude,names(forest_birds_keep))

head(oc_birds_use)
head(forest_birds_use)

### Sample size for each species
colSums(oc_birds_use[3:23])
colSums(forest_birds_use[3:57])

### Remove focal SAR/SOC that have already been analyzed

oc_birds_use <- oc_birds_use[,-c(7,13,15,16,19)]
forest_birds_use <- forest_birds_use[,-c(22,29,54)]


#####################
### Grassland species
#####################

### Create a loop that separates each species into presence/absence and then
### thins each.

### Vector of species
oc_species <- colnames(oc_birds_use[3:length(oc_birds_use)])

### Empty lists
oc_species_data <- list()
oc_species_presence <- list()
oc_species_absence <- list()

### For loop
for(i in 1:length(oc_species)) {
  
  # Select species
  oc_species_data[[i]] <- oc_birds_use %>% select(longitude,latitude,oc_species[i])
  
  # Split into presence and absence
  oc_species_presence[[i]] <- subset(oc_species_data[[i]], oc_species_data[[i]][,3] == 1)
  oc_species_absence[[i]] <- subset(oc_species_data[[i]], oc_species_data[[i]][,3] == 0)
  
  # Create a species column for each
  oc_species_presence[[i]]$species <- rep(oc_species[i], length(oc_species_presence[[i]]$longitude))
  oc_species_absence[[i]]$species <- rep(oc_species[i], length(oc_species_absence[[i]]$longitude))
  
}

head(oc_species_presence[[5]])


### Thin loop

oc_presence_thinned <- list()
oc_absence_thinned <- list()

oc_presence_thinned_df <- list()
oc_absence_thinned_df <- list()

oc_species_combined <- list()

set.seed(212)

for(i in 1:length(oc_species_presence)) {
  
  oc_presence_thinned[[i]] <- thin(oc_species_presence[[i]], 
                                   long.col = "longitude", lat.col = "latitude",
                                   spec.col = "species", thin.par = 1, reps = 10, #1km and 10 replication
                                   locs.thinned.list.return = TRUE,
                                   write.files = FALSE, write.log.file = FALSE, 
                                   verbose = TRUE)
  
  oc_presence_thinned_df[[i]] <- data.frame(oc_presence_thinned[[i]][[1]])
  oc_presence_thinned_df[[i]]$pa <- rep(1,length(oc_presence_thinned_df[[i]]$Longitude))
  colnames(oc_presence_thinned_df[[i]]) <- c("longitude","latitude","pa")
  
  oc_absence_thinned[[i]] <- thin(oc_species_absence[[i]], 
                                  long.col = "longitude", lat.col = "latitude",
                                  spec.col = "species", thin.par = 1, reps = 10, #1km and 10 replication
                                  locs.thinned.list.return = TRUE,
                                  write.files = FALSE, write.log.file = FALSE, 
                                  verbose = TRUE)
  
  oc_absence_thinned_df[[i]] <- data.frame(oc_absence_thinned[[i]][[1]])
  oc_absence_thinned_df[[i]]$pa <- rep(0,length(oc_absence_thinned_df[[i]]$Longitude))
  colnames(oc_absence_thinned_df[[i]]) <- c("longitude","latitude","pa")
  
  oc_species_combined[[i]] <- bind_rows(oc_presence_thinned_df[[i]],oc_absence_thinned_df[[i]])
}

### Save thinned datasets
#saveRDS(oc_species_combined, "oc_species_combined_thinned.rds")


### Run ensemble models on each species

oc_ensemble_models <- list()

for(i in 1:length(oc_species_combined)) {
  
    oc_ensemble_models[[i]] <- ensemble_modelling(c("GLM","GAM","RF"), 
                                 cta.args=list(ntree=5000),
                                 oc_species_combined[[i]], landcover_multiscale_ss,
                                 Pcol = "pa", Xcol = "longitude", Ycol = "latitude",
                                 rep = 1, uncertainty = TRUE, bin.thresh="SES",
                                 cv = "holdout", cv.param = c(0.70,5),
                                 final.fit.data = "holdout",
                                 ensemble.metric = "AUC", ensemble.thresh = 0.50,
                                 weight = TRUE, verbose = TRUE, set.seed=67)
  
  saveRDS(oc_ensemble_models[[i]], paste(oc_species[i],"ta_final_ses.rds", sep="_", collapse=NULL))
  
}


### Reclassify and save binary maps for oc species

oc_thresholds <- c()

oc_projections <- list()
oc_reclassified <- list()
oc_reclassified_final <- list()


for(i in 1:length(oc_ensemble_models)) {
  
  oc_thresholds[i] <- max(oc_ensemble_models[[i]]@algorithm.evaluation$threshold)
  
  oc_projections[[i]] <- oc_ensemble_models[[i]]@projection
  
  oc_reclassified[[i]] <- reclassify(oc_projections[[i]], cbind(oc_thresholds[i], Inf, 1))
  
  oc_reclassified_final[[i]] <- reclassify(oc_reclassified[[i]], cbind(-Inf, oc_thresholds[i], 0))
  
  writeRaster(oc_reclassified_final[[i]],paste(oc_species[i],"reclassified_final.tif", sep="_", collapse=NULL))
  
  print(i)
}

### Read in reclassified binary maps to stack

oc_species <- colnames(oc_birds_use[3:length(oc_birds_use)])

oc_reclassified_list <- list()

for(i in 1:length(oc_species)) {
  
  oc_reclassified_list[[i]] <- raster(paste(oc_species[i],"reclassified_final.tif", sep="_", collapse=NULL))
  
}


#################
### Stack and sum
#################


oc_stack_total <- stack(oc_reclassified_list[[1]],oc_reclassified_list[[2]],
                        oc_reclassified_list[[3]],oc_reclassified_list[[4]],
                        oc_reclassified_list[[5]],oc_reclassified_list[[6]],
                        oc_reclassified_list[[7]],oc_reclassified_list[[8]],
                        oc_reclassified_list[[9]],oc_reclassified_list[[10]],
                        oc_reclassified_list[[11]],oc_reclassified_list[[12]],
                        oc_reclassified_list[[13]],oc_reclassified_list[[14]],
                        oc_reclassified_list[[15]],oc_reclassified_list[[16]],
                        BOBO_reclassified_final,
                        EAME_reclassified_final,HOLA_reclassified_final,
                        SAVS_reclassified_final,KILL_reclassified_final)


names(oc_stack_total)

### Rename levels

names(oc_stack_total) <- c("ALFL","AMCR","BARS","BHCO","BRTH","COGR","COYE","EABL",
                           "EAKI","FISP","MODO","RWBL","SWSP","TRSW","TUVU","VESP",
                           "BOBO","EAME","HOLA","SAVS","KILL")

### Save as raster stack

#writeRaster(oc_stack_total, filename= "oc_stack_total.tif", 
#            format="GTiff", options = "INTERLEAVE=BAND", overwrite=TRUE)


### Read in raster stack
oc_stack_total <- stack("oc_stack_total.tif")

### Calculate richness by summing across species
oc_stack_total_sum <- calc(oc_stack_total, sum)

hist(values(oc_stack_total_sum))

### Mask water
water <- landcover_multiscale_ss$water_200m
oc_total_richness_ssdm <- mask(oc_stack_total_sum, water, maskvalue= 1)
plot(oc_total_richness_ssdm)


### Save
#writeRaster(oc_total_richness_ssdm, "oc_total_richness_ssdm.tif",overwrite=T)



##################
### Forest species
##################


### Create a loop that separates each species into presence/absence and then
### thins each.

### Vector of species
forest_species <- colnames(forest_birds_use[3:length(forest_birds_use)])

### Empty vectors
forest_species_data <- list()
forest_species_presence <- list()
forest_species_absence <- list()

### Loop
for(i in 1:length(forest_species)) {
  
  # Select species
  forest_species_data[[i]] <- forest_birds_use %>% select(longitude,latitude,forest_species[i])
  
  # Split into presence and absence
  forest_species_presence[[i]] <- subset(forest_species_data[[i]], forest_species_data[[i]][,3] == 1)
  forest_species_absence[[i]] <- subset(forest_species_data[[i]], forest_species_data[[i]][,3] == 0)
  
  # Create a species column for each
  forest_species_presence[[i]]$species <- rep(forest_species[i], length(forest_species_presence[[i]]$longitude))
  forest_species_absence[[i]]$species <- rep(forest_species[i], length(forest_species_absence[[i]]$longitude))
  
}

head(forest_species_presence[[5]])


### Thin loop
forest_presence_thinned <- list()
forest_absence_thinned <- list()

forest_presence_thinned_df <- list()
forest_absence_thinned_df <- list()

forest_species_combined <- list()

set.seed(353)

for(i in 1:length(forest_species_presence)) {
  
  forest_presence_thinned[[i]] <- thin(forest_species_presence[[i]], 
                                   long.col = "longitude", lat.col = "latitude",
                                   spec.col = "species", thin.par = 1, reps = 10, #1km and 10 replication
                                   locs.thinned.list.return = TRUE,
                                   write.files = FALSE, write.log.file = FALSE, 
                                   verbose = TRUE)
  
  forest_presence_thinned_df[[i]] <- data.frame(forest_presence_thinned[[i]][[1]])
  forest_presence_thinned_df[[i]]$pa <- rep(1,length(forest_presence_thinned_df[[i]]$Longitude))
  colnames(forest_presence_thinned_df[[i]]) <- c("longitude","latitude","pa")
  
  forest_absence_thinned[[i]] <- thin(forest_species_absence[[i]], 
                                  long.col = "longitude", lat.col = "latitude",
                                  spec.col = "species", thin.par = 1, reps = 10, #1km and 10 replication
                                  locs.thinned.list.return = TRUE,
                                  write.files = FALSE, write.log.file = FALSE, 
                                  verbose = TRUE)
  
  forest_absence_thinned_df[[i]] <- data.frame(forest_absence_thinned[[i]][[1]])
  forest_absence_thinned_df[[i]]$pa <- rep(0,length(forest_absence_thinned_df[[i]]$Longitude))
  colnames(forest_absence_thinned_df[[i]]) <- c("longitude","latitude","pa")
  
  forest_species_combined[[i]] <- bind_rows(forest_presence_thinned_df[[i]],forest_absence_thinned_df[[i]])
}

### Save thinned datasets
#saveRDS(forest_species_combined, "forest_species_combined_thinned.rds")



### Run ensemble models on each species

forest_ensemble_models <- list()

for(i in 1:length(forest_species_combined)) {
  
  forest_ensemble_models[[i]] <- ensemble_modelling(c("GLM","GAM","RF"), 
                                    cta.args=list(ntree=5000),
                                    forest_species_combined[[i]], landcover_2018_forest_sf_final,
                                    Pcol = "pa", Xcol = "longitude", Ycol = "latitude",
                                    rep = 1, uncertainty = TRUE, bin.thresh="SES",
                                    cv = "holdout", cv.param = c(0.70,5),
                                    final.fit.data = "holdout",
                                    ensemble.metric = "AUC", ensemble.thresh = 0.50,
                                    weight = TRUE, verbose = TRUE, set.seed=77)
  
    saveRDS(forest_ensemble_models[[i]], paste(forest_species[i],"ta_final_ses.rds", sep="_", collapse=NULL))
  
}


### Reclassify and save binary maps for oc species

forest_thresholds <- c()

forest_projections <- list()
forest_reclassified <- list()
forest_reclassified_final <- list()


for(i in 1:length(forest_ensemble_models)) {
  
  forest_thresholds[i] <- max(forest_ensemble_models[[i]]@algorithm.evaluation$threshold)
  
  forest_projections[[i]] <- forest_ensemble_models[[i]]@projection
  
  forest_reclassified[[i]] <- reclassify(forest_projections[[i]], cbind(forest_thresholds[i], Inf, 1))
  
  forest_reclassified_final[[i]] <- reclassify(forest_reclassified[[i]], cbind(-Inf, forest_thresholds[i], 0))
  
  writeRaster(forest_reclassified_final[[i]],paste(forest_species[i],"reclassified_final.tif", sep="_", collapse=NULL))
  
  print(i)
}


### Read in reclassified binary maps to stack

forest_species <- colnames(forest_birds_use[3:length(forest_birds_use)])

forest_reclassified_list <- list()

for(i in 1:length(forest_species)) {
  
  forest_reclassified_list[[i]] <- raster(paste(forest_species[i],"reclassified_final.tif", sep="_", collapse=NULL))
  
}


#################
### Stack and sum
#################

forest_stack_total <- stack(forest_reclassified_list[[1]],forest_reclassified_list[[2]],
                            forest_reclassified_list[[3]],forest_reclassified_list[[4]],
                            forest_reclassified_list[[5]],forest_reclassified_list[[6]],
                            forest_reclassified_list[[7]],forest_reclassified_list[[8]],
                            forest_reclassified_list[[9]],forest_reclassified_list[[10]],
                            forest_reclassified_list[[11]],forest_reclassified_list[[12]],
                            forest_reclassified_list[[13]],forest_reclassified_list[[14]],
                            forest_reclassified_list[[15]],forest_reclassified_list[[16]],
                            forest_reclassified_list[[17]],forest_reclassified_list[[18]],
                            forest_reclassified_list[[19]],forest_reclassified_list[[20]],
                            forest_reclassified_list[[21]],forest_reclassified_list[[22]],
                            forest_reclassified_list[[23]],forest_reclassified_list[[24]],
                            forest_reclassified_list[[25]],forest_reclassified_list[[26]],
                            forest_reclassified_list[[27]],forest_reclassified_list[[28]],
                            forest_reclassified_list[[29]],forest_reclassified_list[[30]],
                            forest_reclassified_list[[31]],forest_reclassified_list[[32]],
                            forest_reclassified_list[[33]],forest_reclassified_list[[34]],
                            forest_reclassified_list[[35]],forest_reclassified_list[[36]],
                            forest_reclassified_list[[37]],forest_reclassified_list[[38]],
                            forest_reclassified_list[[39]],forest_reclassified_list[[40]],
                            forest_reclassified_list[[41]],forest_reclassified_list[[42]],
                            forest_reclassified_list[[43]],forest_reclassified_list[[44]],
                            forest_reclassified_list[[45]],forest_reclassified_list[[46]],
                            forest_reclassified_list[[47]],forest_reclassified_list[[48]],
                            forest_reclassified_list[[49]],forest_reclassified_list[[50]],
                            forest_reclassified_list[[51]],forest_reclassified_list[[52]],
                            WOTH_reclassified_final,EWPE_reclassified_final,
                            LEFL_reclassified_final)


names(forest_stack_total)

### Rename levels
names(forest_stack_total) <- c("AMGO","AMRE","AMRO","BAOR","BAWW","BBCU","BCCH",
                               "BHVI","BLBW","BLJA","BTBW","BTNW","CEDW","CHSP",
                               "CORA","CSWA","DOWO","EAPH","EATO","GCFL","GRCA",
                               "HAWO","HETH","HOWR","INBU","MAWA","MOWA","MYWA",
                               "NAWA","NOCA","NOFL","NOPA","NOWA","OVEN","PISI",
                               "PIWA","PIWO","PUFI","RBGR","RBNU","REVI","SCTA",
                               "SOSP","SWTH","VEER","WAVI","WBNU","WITU","WIWR",
                               "WTSP","YBSA","YWAR","WOTH","EWPE","LEFL")

### Save as raster stack for later
#writeRaster(forest_stack_total, filename= "forest_stack_total.tif", 
#            format="GTiff", options = "INTERLEAVE=BAND", overwrite=TRUE)

### Read in raster stack
forest_stack_total <- stack("forest_stack_total.tif")

### Calculate forest species richness by summing across layers
forest_stack_total_sum <- calc(forest_stack_total, sum)


### Mask water
water <- landcover_multiscale_ss$water_200m
forest_total_richness_ssdm <- mask(forest_stack_total_sum, water, maskvalue= 1)
plot(forest_total_richness_ssdm)

### Save
#writeRaster(forest_total_richness_ssdm, "forest_total_richness_ssdm.tif")


####################
### Total species
####################

### 2 km resolution

### Aggregate to 2kmx2km (increase by factor of 3.7)

oc_total_richness_agg <- aggregate(oc_total_richness_ssdm, fact=3.7, fun=max)
forest_total_richness_agg <- aggregate(forest_total_richness_ssdm, fact=3.7, fun=max)

plot(oc_total_richness_agg)
plot(forest_total_richness_agg)

# Stack grassland and forest species richness
total_richness_agg <- stack(oc_total_richness_agg,forest_total_richness_agg)

# Sum across layers
total_richness_agg_sum <- calc(total_richness_agg, sum)

### Save
writeRaster(total_richness_agg_sum, "total_richness_agg_sum_2km_res.tif")


### Plot
png(file="total_richness_2km.png",width=3500,height=2500, res=600)
cuts <- c(0,10,20,30,40,50,60)
pal <- colorRampPalette(viridis(6, option="viridis"))
plot(total_richness_agg_sum, breaks=cuts, col = pal(6),
     xlab = "Longitude", ylab = "Latitude", cex.lab=1.5, xaxt='n', yaxt='n')
axis(side = 2, at=c(42,44,46,48), las=1, cex.axis=1.2)
axis(side = 1, at=c(-82,-78,-74,-70), cex.axis=1.2)
dev.off()



##############################################################################
### Finally, determine where higher than average SAR SR overlaps with total SR
##############################################################################

quantile(sar_stack_agg_sum, na.rm=T) 
### Greater than 75% is greater than 5 for SAR

quantile(total_richness_agg_sum, na.rm=T)
### Greater than 75% is greater than 41 species


### Reclassify SAR and total richness stack

### SAR
SAR_reclassified <- reclassify(sar_stack_agg_sum, cbind(-Inf, 5, 0))
SAR_reclassified_final <- reclassify(SAR_reclassified, cbind(5, Inf, 1))

plot(SAR_reclassified_final)

### Total richness

total_richness_2km_reclassified <- reclassify(total_richness_agg_sum, cbind(-Inf, 41, 0))
total_richness_2km_reclassified_final <- reclassify(total_richness_2km_reclassified, cbind(41, Inf, 1))

plot(total_richness_2km_reclassified_final)


### Save these binary qunatile maps
writeRaster(total_richness_2km_reclassified_final, "total_richness_75_quantile_2km.tif")


### Stack the 2km SR quantile and the SAR to sum together

reclassified_stack <- stack(SAR_reclassified_final,total_richness_2km_reclassified_final)
reclassified_stack_sum <- calc(reclassified_stack, sum)

plot(reclassified_stack_sum)


### Plot
png(file="priority_areas_2km.png",width=3500,height=2500, res=600)
cuts <- c(-1,0,1,2)
pal <- colorRampPalette(viridis(3, option="viridis"))
plot(reclassified_stack_sum, breaks=cuts, col = pal(3),
     xlab = "Longitude", ylab = "Latitude", cex.lab=1.5, xaxt='n', yaxt='n')
axis(side = 2, at=c(42,44,46,48), las=1, cex.axis=1.2)
axis(side = 1, at=c(-82,-78,-74,-70), cex.axis=1.2)
dev.off()



################################################################################
############# Fine and coarse scale presence-absence models ####################
################################################################################

#########################
### Create landcover data
#########################

### A) Open country birds

oc_landcover_fine_ss <- subset(landcover_multiscale_ss, 
                          c("water_200m","shrub_200m",
                            "wetland_200m","corn_200m","soy_200m",
                            "forest_200m","grassland_pasture_200m",
                            "cereals_200m","elevation","slope","aspect",
                            "bios_10","bios_18"))


oc_landcover_coarse_ss <- subset(landcover_multiscale_ss, 
                            c("water_1km","shrub_200m","urban_1km",
                              "wetland_1km","corn_1km","soy_1km",
                              "forest_1km","grassland_pasture_1km",
                              "cereals_1km","elevation","slope","aspect",
                              "bios_10","bios_18"))

### B) Forest birds

fs_landcover_fine_ss <- subset(landcover_2018_forest_sf_final, 
                          c("water_200m","shrub_200m",
                            "crops_200m","wetland_200m","deciduous_200m",
                            "coniferous_200m","mixed_forest_200m",
                            "grassland_pasture_200m","elevation","slope",
                            "aspect","bios_10","bios_18"))

fs_landcover_coarse_ss <- subset(landcover_2018_forest_sf_final, 
                            c("water_1km","shrub_1km","urban_1km",
                              "wetland_1km","deciduous_1km","coniferous_1km",
                              "mixed_forest_1km","grassland_pasture_1km",
                              "crops_1km","elevation","slope","aspect","bios_10",
                              "bios_18"))


############################
### Fit fine scale models
############################


############
### A) BOBO
############

BOBO_ta_final_fine <- ensemble_modelling(c("GLM","GAM","RF"), 
                          cta.args=list(ntree=5000),
                          BOBO_p_a_data, oc_landcover_fine_ss,
                          Pcol = "BOBO_pa", Xcol = "longitude", Ycol = "latitude",
                          rep = 10, uncertainty = TRUE, bin.thresh="SES",
                          cv = "holdout", cv.param = c(0.7,5),
                          final.fit.data = "holdout",
                          ensemble.metric = "AUC", ensemble.thresh = 0.50,
                          weight = TRUE, verbose = TRUE, set.seed=76)

### Save output
saveRDS(BOBO_ta_final_fine, "BOBO_ta_final_fine.rds")


###########
### B) EAME
###########

EAME_ta_final_fine <- ensemble_modelling(c("GLM","GAM","RF"), 
                          cta.args=list(ntree=5000),
                          EAME_p_a_data, oc_landcover_fine_ss,
                          Pcol = "EAME_pa", Xcol = "longitude", Ycol = "latitude",
                          rep = 10, uncertainty = TRUE, bin.thresh="SES",
                          cv = "holdout", cv.param = c(0.7,5),
                          final.fit.data = "holdout",
                          ensemble.metric = "AUC", ensemble.thresh = 0.50,
                          weight = TRUE, verbose = TRUE, set.seed=24)

### Save output
saveRDS(EAME_ta_final_fine, "EAME_ta_final_fine.rds")


###########
### C) HOLA
###########

HOLA_ta_final_fine <- ensemble_modelling(c("GLM","GAM","RF"), 
                          cta.args=list(ntree=5000),
                          HOLA_p_a_data, oc_landcover_fine_ss,
                          Pcol = "HOLA_pa", Xcol = "longitude", Ycol = "latitude",
                          rep = 10, uncertainty = TRUE, bin.thresh="SES",
                          cv = "holdout", cv.param = c(0.7,5),
                          final.fit.data = "holdout",
                          ensemble.metric = "AUC", ensemble.thresh = 0.50,
                          weight = TRUE, verbose = TRUE, set.seed=52)

### Save output
saveRDS(HOLA_ta_final_fine, "HOLA_ta_final_fine.rds")


###########
### D) SAVS
###########

SAVS_ta_final_fine <- ensemble_modelling(c("GLM","GAM","RF"), 
                         cta.args=list(ntree=5000),
                         SAVS_p_a_data, oc_landcover_fine_ss,
                         Pcol = "SAVS_pa", Xcol = "longitude", Ycol = "latitude",
                         rep = 10, uncertainty = TRUE, bin.thresh="SES",
                         cv = "holdout", cv.param = c(0.7,5),
                         final.fit.data = "holdout",
                         ensemble.metric = "AUC", ensemble.thresh = 0.50,
                         weight = TRUE, verbose = TRUE, set.seed=323)

### Save output
saveRDS(SAVS_ta_final_fine, "SAVS_ta_final_fine.rds")


###########
### E) KILL
###########

KILL_ta_final_fine <- ensemble_modelling(c("GLM","GAM","RF"), 
                         cta.args=list(ntree=5000),
                         KILL_p_a_data, oc_landcover_fine_ss,
                         Pcol = "KILL_pa", Xcol = "longitude", Ycol = "latitude",
                         rep = 10, uncertainty = TRUE, bin.thresh="SES",
                         cv = "holdout", cv.param = c(0.7,5),
                         final.fit.data = "holdout",
                         ensemble.metric = "AUC", ensemble.thresh = 0.50,
                         weight = TRUE, verbose = TRUE, set.seed=107)

### Save output
saveRDS(KILL_ta_final_fine, "KILL_ta_final_fine.rds")



###########
### F) WOTH
###########

WOTH_ta_final_fine <- ensemble_modelling(c("GLM","GAM","RF"), 
                         cta.args=list(ntree=5000),
                         WOTH_p_a_data, fs_landcover_fine_ss,
                         Pcol = "WOTH_pa", Xcol = "longitude", Ycol = "latitude",
                         rep = 10, uncertainty = TRUE, bin.thresh="SES",
                         cv = "holdout", cv.param = c(0.7,5),
                         final.fit.data = "holdout",
                         ensemble.metric = "AUC", ensemble.thresh = 0.5,
                         weight = TRUE, verbose = TRUE, set.seed=64)

### Save output
saveRDS(WOTH_ta_final_fine, "WOTH_ta_final_fine.rds")



###########
### G) EWPE
###########

EWPE_ta_final_fine <- ensemble_modelling(c("GLM","GAM","RF"), 
                          cta.args=list(ntree=5000),
                          EWPE_p_a_data, fs_landcover_fine_ss,
                          Pcol = "EWPE_pa", Xcol = "longitude", Ycol = "latitude",
                          rep = 10, uncertainty = TRUE, bin.thresh="SES",
                          cv = "holdout", cv.param = c(0.7,5),
                          final.fit.data = "holdout",
                          ensemble.metric = "AUC", ensemble.thresh = 0.5,
                          weight = TRUE, verbose = TRUE, set.seed=37)

### Save output
saveRDS(EWPE_ta_final_fine, "EWPE_ta_final_fine.rds")


###########
### H) LEFL
###########

LEFL_ta_final_fine <- ensemble_modelling(c("GLM","GAM","RF"), 
                          cta.args=list(ntree=5000),
                          LEFL_p_a_data, fs_landcover_fine_ss,
                          Pcol = "LEFL_pa", Xcol = "longitude", Ycol = "latitude",
                          rep = 10, uncertainty = TRUE, bin.thresh="SES",
                          cv = "holdout", cv.param = c(0.7,5),
                          final.fit.data = "holdout",
                          ensemble.metric = "AUC", ensemble.thresh = 0.5,
                          weight = TRUE, verbose = TRUE, set.seed=919)

### Save output
saveRDS(LEFL_ta_final_fine, "LEFL_ta_final_fine.rds")



############################
### Fit coarse scale models
############################


############
### A) BOBO
############

BOBO_ta_final_coarse <- ensemble_modelling(c("GLM","GAM","RF"), 
                            cta.args=list(ntree=5000),
                            BOBO_p_a_data, oc_landcover_coarse_ss,
                            Pcol = "BOBO_pa", Xcol = "longitude", Ycol = "latitude",
                            rep = 10, uncertainty = TRUE, bin.thresh="SES",
                            cv = "holdout", cv.param = c(0.7,5),
                            final.fit.data = "holdout",
                            ensemble.metric = "AUC", ensemble.thresh = 0.50,
                            weight = TRUE, verbose = TRUE, set.seed=245)

### Save output
saveRDS(BOBO_ta_final_coarse, "BOBO_ta_final_coarse.rds")


###########
### B) EAME
###########

EAME_ta_final_coarse <- ensemble_modelling(c("GLM","GAM","RF"), 
                            cta.args=list(ntree=5000),
                            EAME_p_a_data, oc_landcover_coarse_ss,
                            Pcol = "EAME_pa", Xcol = "longitude", Ycol = "latitude",
                            rep = 10, uncertainty = TRUE, bin.thresh="SES",
                            cv = "holdout", cv.param = c(0.7,5),
                            final.fit.data = "holdout",
                            ensemble.metric = "AUC", ensemble.thresh = 0.50,
                            weight = TRUE, verbose = TRUE, set.seed=53)

### Save output
saveRDS(EAME_ta_final_coarse, "EAME_ta_final_coarse.rds")


###########
### C) HOLA
###########

HOLA_ta_final_coarse <- ensemble_modelling(c("GLM","GAM","RF"), 
                            cta.args=list(ntree=5000),
                            HOLA_p_a_data, oc_landcover_coarse_ss,
                            Pcol = "HOLA_pa", Xcol = "longitude", Ycol = "latitude",
                            rep = 10, uncertainty = TRUE, bin.thresh="SES",
                            cv = "holdout", cv.param = c(0.7,5),
                            final.fit.data = "holdout",
                            ensemble.metric = "AUC", ensemble.thresh = 0.50,
                            weight = TRUE, verbose = TRUE, set.seed=232)

### Save output
saveRDS(HOLA_ta_final_coarse, "HOLA_ta_final_coarse.rds")


###########
### D) SAVS
###########

SAVS_ta_final_coarse <- ensemble_modelling(c("GLM","GAM","RF"), 
                            cta.args=list(ntree=5000),
                            SAVS_p_a_data, oc_landcover_coarse_ss,
                            Pcol = "SAVS_pa", Xcol = "longitude", Ycol = "latitude",
                            rep = 10, uncertainty = TRUE, bin.thresh="SES",
                            cv = "holdout", cv.param = c(0.7,5),
                            final.fit.data = "holdout",
                            ensemble.metric = "AUC", ensemble.thresh = 0.50,
                            weight = TRUE, verbose = TRUE, set.seed=541)

### Save output
saveRDS(SAVS_ta_final_coarse, "SAVS_ta_final_coarse.rds")


###########
### E) KILL
###########

KILL_ta_final_coarse <- ensemble_modelling(c("GLM","GAM","RF"), 
                            cta.args=list(ntree=5000),
                            KILL_p_a_data, oc_landcover_coarse_ss,
                            Pcol = "KILL_pa", Xcol = "longitude", Ycol = "latitude",
                            rep = 10, uncertainty = TRUE, bin.thresh="SES",
                            cv = "holdout", cv.param = c(0.7,5),
                            final.fit.data = "holdout",
                            ensemble.metric = "AUC", ensemble.thresh = 0.50,
                            weight = TRUE, verbose = TRUE, set.seed=111)

### Save output
saveRDS(KILL_ta_final_coarse, "KILL_ta_final_coarse.rds")


###########
### F) WOTH
###########

WOTH_ta_final_coarse <- ensemble_modelling(c("GLM","GAM","RF"), 
                            cta.args=list(ntree=5000),
                            WOTH_p_a_data, fs_landcover_coarse_ss,
                            Pcol = "WOTH_pa", Xcol = "longitude", Ycol = "latitude",
                            rep = 10, uncertainty = TRUE, bin.thresh="SES",
                            cv = "holdout", cv.param = c(0.7,5),
                            final.fit.data = "holdout",
                            ensemble.metric = "AUC", ensemble.thresh = 0.5,
                            weight = TRUE, verbose = TRUE, set.seed=222)

### Save output
saveRDS(WOTH_ta_final_coarse, "WOTH_ta_final_coarse.rds")


###########
### G) EWPE
###########

EWPE_ta_final_coarse <- ensemble_modelling(c("GLM","GAM","RF"), 
                            cta.args=list(ntree=5000),
                            EWPE_p_a_data, fs_landcover_coarse_ss,
                            Pcol = "EWPE_pa", Xcol = "longitude", Ycol = "latitude",
                            rep = 10, uncertainty = TRUE, bin.thresh="SES",
                            cv = "holdout", cv.param = c(0.7,5),
                            final.fit.data = "holdout",
                            ensemble.metric = "AUC", ensemble.thresh = 0.5,
                            weight = TRUE, verbose = TRUE, set.seed=343)

### Save output
saveRDS(EWPE_ta_final_coarse, "EWPE_ta_final_coarse.rds")


###########
### H) LEFL
###########

LEFL_ta_final_coarse <- ensemble_modelling(c("GLM","GAM","RF"), 
                            cta.args=list(ntree=5000),
                            LEFL_p_a_data, fs_landcover_coarse_ss,
                            Pcol = "LEFL_pa", Xcol = "longitude", Ycol = "latitude",
                            rep = 10, uncertainty = TRUE, bin.thresh="SES",
                            cv = "holdout", cv.param = c(0.7,5),
                            final.fit.data = "holdout",
                            ensemble.metric = "AUC", ensemble.thresh = 0.5,
                            weight = TRUE, verbose = TRUE, set.seed=123)

### Save output
saveRDS(LEFL_ta_final_coarse, "LEFL_ta_final_coarse.rds")



################################################################################
### Model evaluation 
################################################################################

### Read in model outputs
BOBO_ta_final_ses <- readRDS("BOBO_ta_final_ses.rds")
EAME_ta_final_ses <- readRDS("EAME_ta_final_ses.rds")
HOLA_ta_final_ses <- readRDS("HOLA_ta_final_ses.rds")
SAVS_ta_final_ses <- readRDS("SAVS_ta_final_ses.rds")
KILL_ta_final_ses <- readRDS("KILL_ta_final_ses.rds")
WOTH_ta_final_ses <- readRDS("WOTH_ta_final_ses.rds")
EWPE_ta_final_ses <- readRDS("EWPE_ta_final_ses.rds")
LEFL_ta_final_ses <- readRDS("LEFL_ta_final_ses.rds")


### Initialize empty lists to extract necessary model data
model_data_list <- list()
test_list <- list()

test_p_list <- list()
test_a_list <- list()

test_p_2018_list <- list()
test_a_2018_list <- list()

species_test_p_list <- list()
species_test_a_list <- list()

predictions_test_p_list <- list()
predictions_test_a_list <- list()


### List of ensemble models
ensemble_model_list <- list(BOBO_ta_final_ses,EAME_ta_final_ses,HOLA_ta_final_ses,
                            SAVS_ta_final_ses,KILL_ta_final_ses,WOTH_ta_final_ses,
                            EWPE_ta_final_ses,LEFL_ta_final_ses)


### Fit loop to extract values

for(j in 1:8) { # For each species
  
  for(i in 1:30) { # For each iteration
    
    #Extract data for each component model
    model_data_list[[i]] <- ensemble_model_list[[j]]@sdms[[i]]@data
    
    #Subset training and testing data
    train_list[[i]] <- subset(model_data_list[[i]], train==TRUE)
    
    #Subset occurrence and absence data from each of the training and testing datasets
    test_p_list[[i]] <- subset(test_list[[i]], Presence==1) 
    test_a_list[[i]] <- subset(test_list[[i]], Presence==0) 
    
    #Convert points from each dataset into spatial data points
    test_p_2018_list[[i]] <- SpatialPointsDataFrame(coords = test_p_list[[i]][,1:2], 
                                                    data = test_p_list[[i]],
                                                    proj4string = CRS("+proj=longlat +ellps=GRS80 +no_defs"))
    test_a_2018_list[[i]] <- SpatialPointsDataFrame(coords = test_a_list[[i]][,1:2], 
                                                    data = test_a_list[[i]],
                                                    proj4string = CRS("+proj=longlat +ellps=GRS80 +no_defs"))
   
    #Extract projection values for each point
    species_test_p_list[[i]] <- extract(ensemble_model_list[[j]]@projection, test_p_2018_list[[i]])
    species_test_a_list[[i]] <- extract(ensemble_model_list[[j]]@projection, test_a_2018_list[[i]])
    
  }
  
  #Put final p/a lists for each species into a larger list with all species
  predictions_test_p_list[[j]] <-  species_test_p_list
  predictions_test_a_list[[j]] <-  species_test_a_list
  
}


### Empty lists to calculate evaluation metrics

predicted_test_presences <- list()
predicted_test_absences <- list()

presence_tables <- list()
absence_tables <- list()

omission_rate <- c()
commission_rate <- c()
sensitivity <- c()
specificity <- c()
TSS <- c()
accuracy <- c()


### Create vector of thresholds
thresholds <- c(0.187,0.165,0.13,0.29,0.18,0.11,0.105,0.095)


########
### BOBO
########

### Fit loop through iterations.

for(i in 1:30) {
  
  predicted_test_presences[[i]] <- ifelse(predictions_test_p_list[[1]][[i]] <= thresholds[1],0,1)
  predicted_test_absences[[i]] <- ifelse(predictions_test_a_list[[1]][[i]] <= thresholds[1],0,1)
  
  presence_tables[[i]] <- as.data.frame(table(predicted_test_presences[[i]]))
  absence_tables[[i]] <- as.data.frame(table(predicted_test_absences[[i]]))
  
  omission_rate[i] <- presence_tables[[i]][1,2]/(presence_tables[[i]][1,2] + presence_tables[[i]][2,2])
  commission_rate[i] <- absence_tables[[i]][2,2]/(absence_tables[[i]][1,2] + absence_tables[[i]][2,2])
  sensitivity[i] <- 1 - omission_rate[i]
  specificity[i] <- 1 - commission_rate[i]
  TSS[i] <-  sensitivity[i] + specificity[i] - 1
  accuracy[i] <- (presence_tables[[i]][2,2] + absence_tables[[i]][1,2])/
    (presence_tables[[i]][2,2] + presence_tables[[i]][1,2] + 
       absence_tables[[i]][1,2] + absence_tables[[i]][2,2])

}

# Assemble data
BOBO_evaluation <- data.frame(omission_rate,commission_rate,sensitivity,
                              specificity,TSS,accuracy)

# Select median metric values
median(BOBO_evaluation$omission_rate, na.rm=T)
median(BOBO_evaluation$commission_rate, na.rm=T)
median(BOBO_evaluation$TSS, na.rm=T)
median(BOBO_evaluation$accuracy, na.rm=T)

# TSS_diff
median(abs(TSS_train - BOBO_evaluation$TSS), na.rm=T)


########
### EAME
########


### Loop through iterations

for(i in 1:30) {
  
  predicted_test_presences[[i]] <- ifelse(predictions_test_p_list[[2]][[i]] <= thresholds[2],0,1)
  predicted_test_absences[[i]] <- ifelse(predictions_test_a_list[[2]][[i]] <= thresholds[2],0,1)
  
  presence_tables[[i]] <- as.data.frame(table(predicted_test_presences[[i]]))
  absence_tables[[i]] <- as.data.frame(table(predicted_test_absences[[i]]))
  
  omission_rate[i] <- presence_tables[[i]][1,2]/(presence_tables[[i]][1,2] + presence_tables[[i]][2,2])
  commission_rate[i] <- absence_tables[[i]][2,2]/(absence_tables[[i]][1,2] + absence_tables[[i]][2,2])
  sensitivity[i] <- 1 - omission_rate[i]
  specificity[i] <- 1 - commission_rate[i]
  TSS[i] <-  sensitivity[i] + specificity[i] - 1
  accuracy[i] <- (presence_tables[[i]][2,2] + absence_tables[[i]][1,2])/
    (presence_tables[[i]][2,2] + presence_tables[[i]][1,2] + 
       absence_tables[[i]][1,2] + absence_tables[[i]][2,2])

}


# Assemble data
EAME_evaluation <- data.frame(omission_rate,commission_rate,sensitivity,
                              specificity,TSS,accuracy)

# Select median metric values
median(EAME_evaluation$omission_rate, na.rm=T)
median(EAME_evaluation$commission_rate, na.rm=T)
median(EAME_evaluation$TSS, na.rm=T)
median(EAME_evaluation$accuracy, na.rm=T)

# TSS_diff
median(abs(TSS_train - EAME_evaluation$TSS), na.rm=T)


########
### HOLA
########

### Loop through iterations

for(i in 1:30) {
  
  predicted_test_presences[[i]] <- ifelse(predictions_test_p_list[[3]][[i]] <= thresholds[3],0,1)
  predicted_test_absences[[i]] <- ifelse(predictions_test_a_list[[3]][[i]] <= thresholds[3],0,1)
  
  presence_tables[[i]] <- as.data.frame(table(predicted_test_presences[[i]]))
  absence_tables[[i]] <- as.data.frame(table(predicted_test_absences[[i]]))
  
  omission_rate[i] <- presence_tables[[i]][1,2]/(presence_tables[[i]][1,2] + presence_tables[[i]][2,2])
  commission_rate[i] <- absence_tables[[i]][2,2]/(absence_tables[[i]][1,2] + absence_tables[[i]][2,2])
  sensitivity[i] <- 1 - omission_rate[i]
  specificity[i] <- 1 - commission_rate[i]
  TSS[i] <-  sensitivity[i] + specificity[i] - 1
  accuracy[i] <- (presence_tables[[i]][2,2] + absence_tables[[i]][1,2])/
    (presence_tables[[i]][2,2] + presence_tables[[i]][1,2] + 
       absence_tables[[i]][1,2] + absence_tables[[i]][2,2])

}


# Assemble data
HOLA_evaluation <- data.frame(omission_rate,commission_rate,sensitivity,
                              specificity,TSS,accuracy)

# Select median metric values
median(HOLA_evaluation$omission_rate, na.rm=T)
median(HOLA_evaluation$commission_rate, na.rm=T)
median(HOLA_evaluation$TSS, na.rm=T)
median(HOLA_evaluation$accuracy, na.rm=T)

# TSS_diff
median(abs(TSS_train - HOLA_evaluation$TSS), na.rm=T)


########
### SAVS
########

### Loop through iterations

for(i in 1:30) {
  
  predicted_test_presences[[i]] <- ifelse(predictions_test_p_list[[4]][[i]] <= thresholds[4],0,1)
  predicted_test_absences[[i]] <- ifelse(predictions_test_a_list[[4]][[i]] <= thresholds[4],0,1)
  
  presence_tables[[i]] <- as.data.frame(table(predicted_test_presences[[i]]))
  absence_tables[[i]] <- as.data.frame(table(predicted_test_absences[[i]]))
  
  omission_rate[i] <- presence_tables[[i]][1,2]/(presence_tables[[i]][1,2] + presence_tables[[i]][2,2])
  commission_rate[i] <- absence_tables[[i]][2,2]/(absence_tables[[i]][1,2] + absence_tables[[i]][2,2])
  sensitivity[i] <- 1 - omission_rate[i]
  specificity[i] <- 1 - commission_rate[i]
  TSS[i] <-  sensitivity[i] + specificity[i] - 1
  accuracy[i] <- (presence_tables[[i]][2,2] + absence_tables[[i]][1,2])/
    (presence_tables[[i]][2,2] + presence_tables[[i]][1,2] + 
       absence_tables[[i]][1,2] + absence_tables[[i]][2,2])

}


# Assemble data
SAVS_evaluation <- data.frame(omission_rate,commission_rate,sensitivity,
                              specificity,TSS,accuracy)

# Select median metric values
median(SAVS_evaluation$omission_rate, na.rm=T)
median(SAVS_evaluation$commission_rate, na.rm=T)
median(SAVS_evaluation$TSS, na.rm=T)
median(SAVS_evaluation$accuracy, na.rm=T)

# TSS_diff
median(abs(TSS_train - SAVS_evaluation$TSS), na.rm=T)


########
### KILL
########

### Loop through iterations

for(i in 1:30) {
  
  predicted_test_presences[[i]] <- ifelse(predictions_test_p_list[[5]][[i]] <= thresholds[5],0,1)
  predicted_test_absences[[i]] <- ifelse(predictions_test_a_list[[5]][[i]] <= thresholds[5],0,1)
  
  presence_tables[[i]] <- as.data.frame(table(predicted_test_presences[[i]]))
  absence_tables[[i]] <- as.data.frame(table(predicted_test_absences[[i]]))
  
  omission_rate[i] <- presence_tables[[i]][1,2]/(presence_tables[[i]][1,2] + presence_tables[[i]][2,2])
  commission_rate[i] <- absence_tables[[i]][2,2]/(absence_tables[[i]][1,2] + absence_tables[[i]][2,2])
  sensitivity[i] <- 1 - omission_rate[i]
  specificity[i] <- 1 - commission_rate[i]
  TSS[i] <-  sensitivity[i] + specificity[i] - 1
  accuracy[i] <- (presence_tables[[i]][2,2] + absence_tables[[i]][1,2])/
    (presence_tables[[i]][2,2] + presence_tables[[i]][1,2] + 
       absence_tables[[i]][1,2] + absence_tables[[i]][2,2])

}

# Assemble data
KILL_evaluation <- data.frame(omission_rate,commission_rate,sensitivity,
                              specificity,TSS,accuracy)

# Select median metric values
median(KILL_evaluation$omission_rate, na.rm=T)
median(KILL_evaluation$commission_rate, na.rm=T)
median(KILL_evaluation$TSS, na.rm=T)
median(KILL_evaluation$accuracy, na.rm=T)

#TSS_diff
median(abs(TSS_train - KILL_evaluation$TSS), na.rm=T)


########
### WOTH
########

### Loop through iterations

for(i in 1:30) {
  
  predicted_test_presences[[i]] <- ifelse(predictions_test_p_list[[6]][[i]] <= thresholds[6],0,1)
  predicted_test_absences[[i]] <- ifelse(predictions_test_a_list[[6]][[i]] <= thresholds[6],0,1)
  
  presence_tables[[i]] <- as.data.frame(table(predicted_test_presences[[i]]))
  absence_tables[[i]] <- as.data.frame(table(predicted_test_absences[[i]]))
  
  omission_rate[i] <- presence_tables[[i]][1,2]/(presence_tables[[i]][1,2] + presence_tables[[i]][2,2])
  commission_rate[i] <- absence_tables[[i]][2,2]/(absence_tables[[i]][1,2] + absence_tables[[i]][2,2])
  sensitivity[i] <- 1 - omission_rate[i]
  specificity[i] <- 1 - commission_rate[i]
  TSS[i] <-  sensitivity[i] + specificity[i] - 1
  accuracy[i] <- (presence_tables[[i]][2,2] + absence_tables[[i]][1,2])/
    (presence_tables[[i]][2,2] + presence_tables[[i]][1,2] + 
       absence_tables[[i]][1,2] + absence_tables[[i]][2,2])
  
}

# Assemble data
WOTH_evaluation <- data.frame(omission_rate,commission_rate,sensitivity,
                              specificity,TSS,accuracy)

# Select median metric values
median(WOTH_evaluation$omission_rate, na.rm=T)
median(WOTH_evaluation$commission_rate, na.rm=T)
median(WOTH_evaluation$TSS, na.rm=T)
median(WOTH_evaluation$accuracy, na.rm=T)

#TSS_diff
median(abs(TSS_train - WOTH_evaluation$TSS), na.rm=T)


########
### EWPE
########

### Loop through iterations

for(i in 1:30) {
  
  predicted_test_presences[[i]] <- ifelse(predictions_test_p_list[[7]][[i]] <= thresholds[7],0,1)
  predicted_test_absences[[i]] <- ifelse(predictions_test_a_list[[7]][[i]] <= thresholds[7],0,1)
  
  presence_tables[[i]] <- as.data.frame(table(predicted_test_presences[[i]]))
  absence_tables[[i]] <- as.data.frame(table(predicted_test_absences[[i]]))
  
  omission_rate[i] <- presence_tables[[i]][1,2]/(presence_tables[[i]][1,2] + presence_tables[[i]][2,2])
  commission_rate[i] <- absence_tables[[i]][2,2]/(absence_tables[[i]][1,2] + absence_tables[[i]][2,2])
  sensitivity[i] <- 1 - omission_rate[i]
  specificity[i] <- 1 - commission_rate[i]
  TSS[i] <-  sensitivity[i] + specificity[i] - 1
  accuracy[i] <- (presence_tables[[i]][2,2] + absence_tables[[i]][1,2])/
    (presence_tables[[i]][2,2] + presence_tables[[i]][1,2] + 
       absence_tables[[i]][1,2] + absence_tables[[i]][2,2])
  
}

# Assemble data
EWPE_evaluation <- data.frame(omission_rate,commission_rate,sensitivity,
                              specificity,TSS,accuracy)

# Select median metric values
median(EWPE_evaluation$omission_rate, na.rm=T)
median(EWPE_evaluation$commission_rate, na.rm=T)
median(EWPE_evaluation$TSS, na.rm=T)
median(EWPE_evaluation$accuracy, na.rm=T)

#TSS_diff
median(abs(TSS_train - EWPE_evaluation$TSS), na.rm=T)


########
### LEFL
########

### Loop through iterations

for(i in 1:30) {
  
  predicted_test_presences[[i]] <- ifelse(predictions_test_p_list[[8]][[i]] <= thresholds[8],0,1)
  predicted_test_absences[[i]] <- ifelse(predictions_test_a_list[[8]][[i]] <= thresholds[8],0,1)
  
  presence_tables[[i]] <- as.data.frame(table(predicted_test_presences[[i]]))
  absence_tables[[i]] <- as.data.frame(table(predicted_test_absences[[i]]))
  
  omission_rate[i] <- presence_tables[[i]][1,2]/(presence_tables[[i]][1,2] + presence_tables[[i]][2,2])
  commission_rate[i] <- absence_tables[[i]][2,2]/(absence_tables[[i]][1,2] + absence_tables[[i]][2,2])
  sensitivity[i] <- 1 - omission_rate[i]
  specificity[i] <- 1 - commission_rate[i]
  TSS[i] <-  sensitivity[i] + specificity[i] - 1
  accuracy[i] <- (presence_tables[[i]][2,2] + absence_tables[[i]][1,2])/
    (presence_tables[[i]][2,2] + presence_tables[[i]][1,2] + 
       absence_tables[[i]][1,2] + absence_tables[[i]][2,2])

}

# Assemble data
LEFL_evaluation <- data.frame(omission_rate,commission_rate,sensitivity,
                              specificity,TSS,accuracy)

# Select median metric values
median(LEFL_evaluation$omission_rate, na.rm=T)
median(LEFL_evaluation$commission_rate, na.rm=T)
median(LEFL_evaluation$TSS, na.rm=T)
median(LEFL_evaluation$accuracy, na.rm=T)

#TSS_diff
median(abs(TSS_train - LEFL_evaluation$TSS), na.rm=T)



################################
### Calculate AUC for each model
################################

### Empty lists
predicted_test_presences <- list()
predicted_test_absences <- list()

predicted_presences <- list()
predicted_absences <- list()

predicted <- list()

model_data_list <- list()
test_list <- list()

pred <- list()
auc_tmp<- list()
auc <- list()

auc_all_species <- list()

### Filled lists
ensemble_model_list <- list(BOBO_ta_final_ses,EAME_ta_final_ses,HOLA_ta_final_ses,
                            SAVS_ta_final_ses,KILL_ta_final_ses,WOTH_ta_final_ses,
                            EWPE_ta_final_ses,LEFL_ta_final_ses)


### Loop through to calculate AUC across species and iterations

for(j in 1:8) {
  
  for(i in 1:30) {
    
    #Convert projection into binary
    predicted_test_presences[[i]] <- ifelse(predictions_test_p_list[[j]][[i]] <= thresholds[j],0,1)
    predicted_test_absences[[i]] <- ifelse(predictions_test_a_list[[j]][[i]] <= thresholds[j],0,1)
    
    #Convert test presences and absences to dataframe
    predicted_presences[[i]] <- as.data.frame(predicted_test_presences[[i]])
    predicted_absences[[i]] <- as.data.frame(predicted_test_absences[[i]])
    
    #Rename columns
    colnames(predicted_presences[[i]]) <- "predicted"
    colnames(predicted_absences[[i]])  <- "predicted"
    
    ### Combine predictions
    predicted[[i]] <- bind_rows(predicted_presences[[i]],predicted_absences[[i]])
    
    #Grab model data
    model_data_list[[i]] <- ensemble_model_list[[j]]@sdms[[i]]@data
    
    #Subset training and testing data
    train_list[[i]] <- subset(model_data_list[[i]], train==TRUE)
    
    ### Combine with occurrence data
    test_list[[i]]$predicted <- predicted[[i]][,1]
    
    #Calculate AUC
    pred[[i]] <- prediction(test_list[[i]]$predicted, test_list[[i]]$Presence)
    auc_tmp[[i]] <- performance(pred[[i]],"auc")
    auc[[i]] <- as.numeric(auc_tmp[[i]]@y.values)
    
}
  
  #Put auc vector into list across all species
  auc_all_species[[j]] <-  auc
  auc_all_species_train[[j]] <-  auc_train
  
}

########
### Print AUC values
########

median(unlist(auc_all_species[[1]])) # BOBO
median(unlist(auc_all_species[[2]])) # EAME
median(unlist(auc_all_species[[3]])) # HOLA
median(unlist(auc_all_species[[4]])) # SAVS
median(unlist(auc_all_species[[5]])) # KILL
median(unlist(auc_all_species[[6]])) # WOTH
median(unlist(auc_all_species[[7]])) # EWPE
median(unlist(auc_all_species[[8]])) # LEFL


###########
### AUC dif
###########

# Create empty vectors and lists

AUC_dif <- c()
AUC_dif_all_species <- list()

# Fit loop to extract difference values

for(j in 1:8) {
  
  for(i in 1:30) {
    
    AUC_dif[i] <- median(c(auc_all_species_train[[j]][[i]] - auc_all_species[[j]][[i]]))
    
  }

  AUC_dif_all_species[[j]] <- AUC_dif
}

### Print AUC dif values

median(unlist(AUC_dif_all_species[[1]]))
median(unlist(AUC_dif_all_species[[2]]))
median(unlist(AUC_dif_all_species[[3]]))
median(unlist(AUC_dif_all_species[[4]]))
median(unlist(AUC_dif_all_species[[5]]))
median(unlist(AUC_dif_all_species[[6]]))
median(unlist(AUC_dif_all_species[[7]]))
median(unlist(AUC_dif_all_species[[8]]))



#################################
### Calculate 5yr commission rate
#################################

### Extract 5yr absence data for each predicted point
BOBO_predictions_a_5yr <- extract(BOBO_ta_final_ses@projection, BOBO_TA_use)
EAME_predictions_a_5yr <- extract(EAME_ta_final_ses@projection, EAME_TA_use)
HOLA_predictions_a_5yr <- extract(HOLA_ta_final_ses@projection, HOLA_TA_use)
SAVS_predictions_a_5yr <- extract(SAVS_ta_final_ses@projection, SAVS_TA_use)
KILL_predictions_a_5yr <- extract(KILL_ta_final_ses@projection, KILL_TA_use)
WOTH_predictions_a_5yr <- extract(WOTH_ta_final_ses@projection, WOTH_TA_use)
EWPE_predictions_a_5yr <- extract(EWPE_ta_final_ses@projection, EWPE_TA_use)
LEFL_predictions_a_5yr <- extract(LEFL_ta_final_ses@projection, LEFL_TA_use)

### Combine into a list
predictions_a_5yr <- list(BOBO_predictions_a_5yr,EAME_predictions_a_5yr,
                          HOLA_predictions_a_5yr,SAVS_predictions_a_5yr,
                          KILL_predictions_a_5yr,WOTH_predictions_a_5yr,
                          EWPE_predictions_a_5yr,LEFL_predictions_a_5yr)

### Create empty lists and vectors
predicted_absences_5yr <- list()
confusion_matrices <- list()

cr_5yr <- c()

### Fit loop to calculate 5yr commission rate

for(i in 1:8) {
  
  # Convert to binary using most conservative threshold
  predicted_absences_5yr[[i]] <- ifelse(predictions_a_5yr[[i]] <= thresholds[i],0,1)
  
  # Confusion matrices
  confusion_matrices[[i]] <- table(predicted_absences_5yr[[i]])
  
  # Calculate commission rate
  cr_5yr[i] <- confusion_matrices[[i]][2]/(confusion_matrices[[i]][2]+confusion_matrices[[i]][1])
}

## Print values

cr_5yr


###########
### Note - to evaluate fine and coarse scale models, run them through the
### above code as well.
###########


################################################################################
############################## End of code #####################################
################################################################################
