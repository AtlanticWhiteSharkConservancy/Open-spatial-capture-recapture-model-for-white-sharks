############################################################################
### Construct spatial grid
############################################################################

##Load libraries##
library(readxl)
library(tidyverse)
library(tidyr)
library(sf)
library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(viridis)
library(gridExtra)
library(rasterVis)
library(lwgeom)
library(ggspatial)
# devtools::install_github("ehrscape/R-project/AnalyzeGPS")
library(analyzeGPS)
library(plot3D)
library(rgl)
library(lunar)
library(oce)
library(here)
here( )
# read in support functions
source("Gridding_Support_Functions.R")

##Add Massachusetts shapefile##
ma <- st_read("data/OUTLINE25K_POLY.shp")
mass <- st_as_sf(ma)

##Survey area##
survey_area <- readOGR("data/Survey Area.shp")
proj4string(survey_area)
# # Convert data to MASP 26986 (MAss State Plane, in meters)
survey_area <- spTransform(survey_area, CRS("+init=epsg:26986"))
# # Convert to sf object to calculate area
sa <- st_as_sf(survey_area)
st_area(sa)
# Convert to km because very large in m
st_area(sa)/1000/1000

# need to add a buffer to study area because clipping out sightings close to shore
# create the buffers object
sa <- st_buffer(sa, 25*1000) # 350 m buffer contains all

# clip by Mass to crop out land
# first insert Mass a bit so no one gets cut out
st_erase = function(x, y) st_difference(x, st_union(st_combine(y)))
sa <- st_erase(sa, mass)
plot(sa)
sa <- st_buffer(sa, 350)

# convert back
survey_area <- as(sa,'Spatial')
##Establish desired area##
area_km2 <- 100 #100 km2 gridding to start
area_m2 <- km2_to_m2(area_km2) 

##Hexagonal - with clip##
hex_grid_c <- make_grid(survey_area, type = "hexagonal", cell_area = area_m2, clip = TRUE) 
plot(survey_area, col = "grey50", bg = "light blue", axes = FALSE)
plot(hex_grid_c, border = "orange", add = TRUE)
box()

# make sure area matches with survey area
# In meters squared
sum( area(hex_grid_c) )
# In kilometers squared
sum( area(hex_grid_c) )/1000/1000

# Make clipped hex grid and sf object for other operations
hexsf <- as( hex_grid_c, 'sf' )

# calc distance from the coast

# read in Cape survey area for clipping gps tracks in next step
survey_poly <- readOGR("data/Survey Area.shp")
proj4string(survey_poly)
# # Convert data to MASP 26986 (MAss State Plane, in meters)
survey_poly <- spTransform(survey_poly, CRS("+init=epsg:26986"))
# # Convert to sf object to calculate area
spoly <- st_as_sf(survey_poly)
st_area(spoly)/1000/1000
# need to add a buffer to study area so doesn't clip out sightings/tracks close to shore
spoly <- st_buffer(spoly, 350) # 350 m buffer contains all
st_area(spoly)
st_area(hexsf)/1000/1000

#End#