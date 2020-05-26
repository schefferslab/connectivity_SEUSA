library(ggmap)
library(readr)
library(maptools)
library(sp)
library(raster)
library(rgdal)
library(maps)
library(ggplot2)
library(ggnewscale)
library(smoothr)
library(tidyverse)
library(devtools)
library(gdistance)
library(ggpubr)
library(spatstat)
library(rgeos)
library(ggthemes)
library(viridis)
library(spatialEco)
library(ResistanceGA)
library(devtools)


rm(list = ls())

#setwd("C:/Users/jbaecher/Dropbox (UFL)/UF/Scheffers lab/Projects/Landscape_prioritization")

HSI <- read_rds("C:/Users/jbaecher/Dropbox (UFL)/UF/Scheffers lab/Projects/Landscape_prioritization/HSI_D_aur.rds")
proj4string <- crs(HSI)
plot(HSI)

HSI_tr <- transition(HSI, transitionFunction = mean, 8) %>%
  geoCorrection(type="c",multpl=F)

################################################################################################## 
#################################### Cleaning Gbif Data ##########################################
################################################################################################## 
D_aur <- readRDS("C:/Users/jbaecher/Dropbox (UFL)/UF/Courses/Spring 2020/SDM_project/Project/D_aur.rds")
D_aur_df <- data.frame(D_aur$gbif$data)
# Begin filtering low quality data:
D_aur_fltr <- D_aur_df %>%                                                  # filter out:
  filter(!is.na(Desmognathus_auriculatus.year)) %>%                         # na's in year
  filter(!Desmognathus_auriculatus.year < 1960) %>%                         # obs before 1960
  filter(!Desmognathus_auriculatus.year > 2019) %>%                         # obs before 1960
  filter(!is.na(Desmognathus_auriculatus.longitude)) %>%                    # na's in long
  filter(!is.na(Desmognathus_auriculatus.latitude)) %>%                     # na's in lat
  filter(!round(Desmognathus_auriculatus.longitude)==0) %>%                 # 0's in long
  filter(!round(Desmognathus_auriculatus.latitude)==0) %>%                  # 0's in lat
  filter(!duplicated(round(Desmognathus_auriculatus.longitude,1),           # duplicate longitudes
                     round(Desmognathus_auriculatus.latitude,1)) == TRUE) %>%        # duplicate latitudes
  dplyr::rename(lon = Desmognathus_auriculatus.longitude,                   # renaming coordinate columns
                lat = Desmognathus_auriculatus.latitude)

# Removing points outside the US
D_aur_spdf   <- SpatialPointsDataFrame(dplyr::select(D_aur_fltr,    # create new spdf with removed points
                                                     lon, 
                                                     lat),
                                       data=D_aur_fltr, proj4string = proj4string)
HSI_sp <- as(extent(HSI),"SpatialPolygons")
crs(HSI_sp) <- crs(HSI)

D_aur_over <- over(D_aur_spdf, HSI_sp)    
D_aur_spdf <- D_aur_spdf[!is.na(D_aur_over),]

# Extract data for plotting
D_aur_xy <- data.frame(lon=    D_aur_spdf$lon, 
                       lat=    D_aur_spdf$lat,
                       decade= 10*floor(D_aur_spdf$Desmognathus_auriculatus.year/10))
D_aur_spdf   <- SpatialPointsDataFrame(dplyr::select(D_aur_xy,    # create new spdf with removed points
                                                     lon, 
                                                     lat),
                                       data=D_aur_xy, proj4string = proj4string)
ggplot(D_aur_xy) + geom_point(aes(x=lon, y=lat)) + facet_wrap(~decade)

################################################################################################## 
########################################### 1960 #################################################
################################################################################################## 
combn_1960 <- combn(nrow(D_aur_xy[D_aur_xy$decade == 1960,1:2]),2)
combn_1960_tr <- as.matrix(t(combn_1960))

commutes_1960 <- list()
system.time(
  for (i in 1:nrow(combn_1960_tr)) {
    sites_1960 <- SpatialPoints(rbind(D_aur_spdf[combn_1960_tr[i,1],1:2],
                                      D_aur_spdf[combn_1960_tr[i,2],1:2]),
                                proj4string)
    commutes_1960[[i]] <- passage(HSI_tr, 
                                  origin=sites_1960[1],
                                  goal=sites_1960[2], theta = 0.00001) 
  }
)
commutes_1960_stack <- stack(commutes_1960)
commutes_1960_overlay <- scale(sum(commutes_1960_stack)/nrow(combn_1960_tr))
plot(commutes_1960_overlay, legend=F)
write_rds(commutes_1960_overlay, "D_aur_1960_cum_cir_HSI.rds")

################################################################################################## 
########################################### 1970 #################################################
################################################################################################## 
combn_1970 <- combn(nrow(D_aur_xy[D_aur_xy$decade == 1970,1:2]),2)
combn_1970_tr <- as.matrix(t(combn_1970))

commutes_1970 <- list()
system.time(
  for (i in 1:nrow(combn_1970_tr)) {
    sites_1970 <- SpatialPoints(rbind(D_aur_spdf[combn_1970_tr[i,1],1:2],
                                      D_aur_spdf[combn_1970_tr[i,2],1:2]),
                                proj4string)
    commutes_1970[[i]] <- passage(HSI_tr, 
                                  origin=sites_1970[1],
                                  goal=sites_1970[2], theta = 0.00001) 
  }
)
commutes_1970_stack <- stack(commutes_1970)
commutes_1970_overlay <- scale(sum(commutes_1970_stack)/nrow(combn_1970_tr))
plot(commutes_1970_overlay, legend=F)
write_rds(commutes_1970_overlay, "D_aur_1970_cum_cir_HSI.rds")

################################################################################################## 
########################################### 1980 #################################################
################################################################################################## 
combn_1980 <- combn(nrow(D_aur_xy[D_aur_xy$decade == 1980,1:2]),2)
combn_1980_tr <- as.matrix(t(combn_1980))

commutes_1980 <- list()
system.time(
  for (i in 1:nrow(combn_1980_tr)) {
    sites_1980 <- SpatialPoints(rbind(D_aur_spdf[combn_1980_tr[i,1],1:2],
                                      D_aur_spdf[combn_1980_tr[i,2],1:2]),
                                proj4string)
    commutes_1980[[i]] <- passage(HSI_tr, 
                                  origin=sites_1980[1],
                                  goal=sites_1980[2], theta = 0.00001) 
  }
)
commutes_1980_stack <- stack(commutes_1980)
commutes_1980_overlay <- scale(sum(commutes_1980_stack)/nrow(combn_1980_tr))
plot(commutes_1980_overlay, legend=F)
write_rds(commutes_1980_overlay, "D_aur_1980_cum_cir_HSI.rds")

################################################################################################## 
########################################### 1990 #################################################
################################################################################################## 
combn_1990 <- combn(nrow(D_aur_xy[D_aur_xy$decade == 1990,1:2]),2)
combn_1990_tr <- as.matrix(t(combn_1990))

commutes_1990 <- list()
system.time(
  for (i in 1:nrow(combn_1990_tr)) {
    sites_1990 <- SpatialPoints(rbind(D_aur_spdf[combn_1990_tr[i,1],1:2],
                                      D_aur_spdf[combn_1990_tr[i,2],1:2]),
                                proj4string)
    commutes_1990[[i]] <- passage(HSI_tr, 
                                  origin=sites_1990[1],
                                  goal=sites_1990[2], theta = 0.00001) 
  }
)
commutes_1990_stack <- stack(commutes_1990)
commutes_1990_overlay <- scale(sum(commutes_1990_stack)/nrow(combn_1990_tr))
plot(commutes_1990_overlay, legend=F)
write_rds(commutes_1990_overlay, "D_aur_1990_cum_cir_HSI.rds")


################################################################################################## 
########################################### 2000 #################################################
################################################################################################## 
combn_2000 <- combn(nrow(D_aur_xy[D_aur_xy$decade == 2000,1:2]),2)
combn_2000_tr <- as.matrix(t(combn_2000))
number <- round(runif(1,1,nrow(combn_2000_tr)),0)
sites_2000 <- SpatialPoints(rbind(D_aur_spdf[combn_2000_tr[number,1],1:2],
                                  D_aur_spdf[combn_2000_tr[number,2],1:2]),
                            proj4string)
plot(HSI); points(sites_2000)
test_2000 <- passage(HSI_tr, origin=sites_2000[1], goal=sites_2000[2], theta = 0.0001) 
ggplot(as.data.frame(test_2000,xy=T)) + geom_raster(aes(x=x,y=y,fill=layer)) + scale_fill_viridis_c(na.value=NA)

commutes_2000 <- list()
system.time(
  for (i in 1:nrow(combn_2000_tr)) {
    sites_2000 <- SpatialPoints(rbind(D_aur_spdf[combn_2000_tr[i,1],1:2],
                                      D_aur_spdf[combn_2000_tr[i,2],1:2]),
                                proj4string)
    commutes_2000[[i]] <- passage(HSI_tr, 
                                  origin=sites_2000[1],
                                  goal=sites_2000[2], theta = 0.00001) 
  }
)
commutes_2000_stack <- stack(commutes_2000)
commutes_2000_overlay <- scale(sum(commutes_2000_stack)/nrow(combn_2000_tr))
plot(commutes_2000_overlay, legend=F)
write_rds(commutes_2000_overlay, "D_aur_2000_cum_cir_HSI.rds")

################################################################################################## 
########################################### 2010 #################################################
################################################################################################## 
combn_2010 <- combn(nrow(D_aur_xy[D_aur_xy$decade == 2010,1:2]),2)
combn_2010_tr <- as.matrix(t(combn_2010))

commutes_2010 <- list()
system.time(
  for (i in 1:nrow(combn_2010_tr)) {
    sites_2010 <- SpatialPoints(rbind(D_aur_spdf[combn_2010_tr[i,1],1:2],
                                      D_aur_spdf[combn_2010_tr[i,2],1:2]),
                                proj4string)
    commutes_2010[[i]] <- passage(HSI_tr, 
                                  origin=sites_2010[1],
                                  goal=sites_2010[2], theta = 0.00001) 
  }
)
commutes_2010_stack <- stack(commutes_2010)
commutes_2010_overlay <- scale(sum(commutes_2010_stack)/nrow(combn_2010_tr))
plot(commutes_2010_overlay, legend=F)
write_rds(commutes_2010_overlay, "D_aur_2010_cum_cir_HSI.rds")


commutes_stack <- stack(commutes_1960_overlay,
                        commutes_1970_overlay,
                        commutes_1980_overlay,
                        commutes_1990_overlay,
                        commutes_2000_overlay,
                        commutes_2010_overlay)
names(commutes_stack) <- c("com_1960","com_1970","com_1980","com_1990","com_2000","com_2010")
plot(commutes_stack,zlim=c(0,50))

rds_1960 <- read_rds("C:/Users/jbaecher/Dropbox (UFL)/UF/Scheffers lab/Projects/Landscape_prioritization/D_aur_1960_cum_cir_HSI.rds")
rds_1970 <- read_rds("C:/Users/jbaecher/Dropbox (UFL)/UF/Scheffers lab/Projects/Landscape_prioritization/D_aur_1970_cum_cir_HSI.rds")
rds_1980 <- read_rds("C:/Users/jbaecher/Dropbox (UFL)/UF/Scheffers lab/Projects/Landscape_prioritization/D_aur_1980_cum_cir_HSI.rds")
rds_1990 <- read_rds("C:/Users/jbaecher/Dropbox (UFL)/UF/Scheffers lab/Projects/Landscape_prioritization/D_aur_1990_cum_cir_HSI.rds")
rds_2000 <- read_rds("C:/Users/jbaecher/Dropbox (UFL)/UF/Scheffers lab/Projects/Landscape_prioritization/D_aur_2000_cum_cir_HSI.rds")
rds_2010 <- read_rds("C:/Users/jbaecher/Dropbox (UFL)/UF/Scheffers lab/Projects/Landscape_prioritization/D_aur_2010_cum_cir_HSI.rds")
