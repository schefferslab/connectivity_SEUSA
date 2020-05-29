library(readr)
library(maptools)
library(sp)
library(raster)
library(rgdal)
library(maps)
library(ggplot2)
library(tidyverse)
library(smoothr)
library(devtools)
library(ENMTools)
library(spatstat)
library(rgeos)
library(lattice)
library(reshape2)
library(ggthemes)
library(viridis)
library(gridExtra)
library(rasterVis)
library(spatialEco)
library(dismo)

rm(list = ls())

setwd("C:/Users/jbaecher/Dropbox (UFL)/UF/Scheffers lab/Projects/Landscape_prioritization")

HSI <- read_rds("C:/Users/jbaecher/Dropbox (UFL)/UF/Scheffers lab/Projects/Landscape_prioritization/HSI_D_aur.rds")
raster <- read_rds("data/circuitscape/se_raster_stack.rds") %>%
  crop(HSI) %>%
  mask(HSI)
proj4string <- crs(raster)
plot(raster)

# write_rds(D_aur_spdf, "D_aur_spdf.rds")
# write_rds(D_aur_xy, "D_aur_xy.rds")
D_aur_spdf <- read_rds("D_aur_spdf.rds")                                        # Load GBIF data (precleaned), spdf version
D_aur_xy   <- read_rds("D_aur_xy.rds") 
# Load df version
D_aur_coords_decade <- split(cbind(D_aur_xy[,c('lon','lat')]), 
                             f = D_aur_xy$decade)                               # Split occurrence data by decade
names(D_aur_coords_decade[2:4]) <- names(D_aur_coords_decade[1])

############# NLCD processing ############# 
# Original values:
# 1:  	Open Water
# 2:  	Urban/Developed
# 3:  	Intentionally Left Blank 
# 4:  	Intentionally Left Blank 
# 5:  	Intentionally Left Blank 
# 6:   	Mining
# 7: 	  Barren
# 8: 	  Deciduous Forest
# 9: 	  Evergreen Forest
# 10:	  Mixed Forest
# 11: 	Grassland
# 12: 	Shrubland
# 13: 	Cultivated Cropland
# 14: 	Hay/Pasture
# 15: 	Herbaceous Wetland
# 16: 	Woody Wetland
# 17: 	Perennial Ice/Snow

# Reclassify rasters based on resistance to salamander movement 
## First, create reclassification matrix
reclassify <- matrix(
  c(    0,   1,   7,      #Open water                    7   
        1,   6,   6,      #Developed                     6
        6,   7,   5,      #Barren                        5
        7,  10,   1,      #Forest                        1
       10,  12,   3,      #Grassland/Shrubland           3
       12,  14,   4,      #Planted/cultivated cropland   4
       14,  16,   2,      #Wetlands                      2
       16, Inf,   NA),    #Ice/Snow (none in SE)         NA
   ncol=3, byrow=T)

# Create labels for attribute table manipulations during raster processing 
nlcdclass <- c("Forest", "Wetlands", "Grass/Shrub/Herbaceous", "Planted/Cultivated", "Barren", "Developed", "Open Water")
classdf <- data.frame(classvalue1 = c(1,2,3,4,5,6,7), classnames1 = nlcdclass)

# Load in  backcasted NLCD rasters to project, reclassify, crop, mask, and stack
## First, identify files
all_rasters = list.files(
  path="C:/Users/jbaecher/Dropbox (UFL)/UF/Scheffers lab/Projects/Landscape_prioritization/data/NLCD/Historic",
  pattern = "\\.tif$", 
  full.names = TRUE)

## Next, remove files outside study time period
sub_rasters <- all_rasters[c(23:55)]
sub_rasters_list <- list()

decades <- c("CONUS_Backcasting_y196",
             "CONUS_Backcasting_y197",
             "CONUS_Backcasting_y198",
             "CONUS_Backcasting_y199")

decade_layers <- list()
decade_modes <- list()

# Big as muthafuckin' for loop to... load, reproject, reclassify, crop, mask, stack,and calculate modes of NLCD rasters from 1960 until 1992

for (i in 1:length(sub_rasters)){
  sub_rasters_list[[i]] <- raster(sub_rasters[i]) %>%                       # Loading rasters in from file names
    projectRaster(HSI) %>%                                                  # Reprojecting 
    crop(HSI) %>%                                                           # Croping to study extent
    mask(HSI) %>%                                                           # Masking to study polygon
    reclassify(reclassify) %>%                                              # Reclassifying habitat types based on resistance value
    ratify()                                                                # Reorganizing attributes table 
  rat <- levels(sub_rasters_list[[i]])                                      # Setting levels of raster values
  rat$landcover <- nlcdclass                                                # Assigning category names to attribute table
  levels(sub_rasters_list[[i]]) <- rat                                      # Saving attribute table in raster object
  if(length(sub_rasters_list) == length(sub_rasters)){                      # Testing if raster processing is complete
    print(names(sub_rasters_list[[i]]))                                     # If raster processing is complete, printing final raster name
    sub_rasters_stack <- stack(sub_rasters_list)                            # If raster processing is complete, stacking list into a raster brick
    decade_layers <- lapply(                                                # Begin function to...
      decades, function(x)                                                  # use a list of decades from study period...
        which(grepl(tolower(x),                                             # to find raster years
                    tolower(names(sub_rasters_stack)))))                    # and return those rasters to decadal groupings
    decade_layers[-c(5,6)]                                                  # removing unnecessary elements from list
    for (j in 1:length(decade_layers)){                                     # Begin for loop to...
      decade_modes[[j]] <- modal(sub_rasters_stack[[                        # calculate the mode of a raster brick...
        decade_layers[[j]]                                                  # across layers...
        ]],                                                                 # representing decade groupings from previous...
        ties="random",freq=F)                                               # setting ties to a random outcome
      if(length(decade_modes) == length(decade_layers)){                    # Testing if mode calculation is complete
        print(names(decade_modes[[j]]))                                     # If mode calculation is complete, printing final name of calculated mode layer
        decade_stack <- stack(decade_modes)                                 # If mode calculation is complete, stacking list of mode layers into a raster brick
        names(decade_stack) <- c("NLCD_1960_1969","NLCD_1970_1979",
                                 "NLCD_1980_1989","NLCD_1990_1992")
      } else(                                                               # If mode calculation is incomplete...
        print(names(decade_modes[[j]])))                                    # print progress...
    }
  } else(                                                                   # If raster processing is incomplete...
    print(names(sub_rasters_list[[i]])))                                    # print progress...
}

barplot(decade_modes[[1]],axes=F,col=plasma(7));axis(1,labels=nlcdclass,at=c(1:7))
################################################# Maxent #################################################
plot(decade_stack)
##
# args to pass to Maxent 
args <- list(
  c("-J", "-P", "-q", "-p", "-h", "replicates=3", "randomtestpoints=27", "betamultiplier=1", 
    "askoverwrite=false", "threads=6"),
  c("-J", "-P", "-q", "-p", "-h", "replicates=3", "randomtestpoints=38", "betamultiplier=1", 
    "askoverwrite=false", "threads=6"),
  c("-J", "-P", "-q", "-p", "-h", "replicates=3", "randomtestpoints=8", "betamultiplier=1", 
    "askoverwrite=false", "threads=6"),
  c("-J", "-P", "-q", "-p", "-h", "replicates=3", "randomtestpoints=5", "betamultiplier=1", 
    "askoverwrite=false", "threads=6"))

D_aur_maxent_list <- list()
D_aur_preds_list <- list()
for(k in 1:nlayers(decade_stack)){
   D_aur_maxent_list[[k]] <- maxent(x=stack(decade_stack[[k]],raster),       # Run MaxEnt on each decade
                                p=coordinates(D_aur_coords_decade[[k]]),     # Partition occurrence by decade
                                args=args[[k]])                              # Pass decade-specific arguments to MaxEnt
  if(length(D_aur_maxent_list) == nlayers(decade_stack)){                   # Test if Maxent calculations are done
    D_aur_preds_list[[k]] <- mean(predict(D_aur_maxent_list[[k]],                # If complete, calculate predictions...
                                     stack(decade_stack[[k]],raster)))       # from decadal raster data...
  } else(                                                                   # If incomplete, 
    print(names(decade_modes[[k]])))                                        # print progress...
}

D_aur_preds_stack <- stack(D_aur_preds_list)
plot(D_aur_preds_stack, zlim=c(0,1))

response(D_aur_maxent_list[[1]], var="NLCD_1960_1969")

         