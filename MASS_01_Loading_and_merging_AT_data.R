## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## 
## INTRODUCTION TO THE ANALYSIS OF ACOUSTIC TELEMETRY DATA IN R
## 1. LOADING AND MERGING AT DATA
##
## Mallorca Science School (MASS) - Palma (Balearic Islands)
##
## Author: Eneko Aspillaga (IMEDEA, CSIC-UIB, Spain)
## Contact: aspillaga@imedea.uib-csic.es
## Date: October 21-25 2024
## 
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Set working directory (directory of the current script)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Install required libraries:
# install.packages(c("data.table", "lubridate", "sp", "sf", "igraph"))

# Load libraries
library(data.table) # Extension of data.frame for faster management
library(lubridate) # Dealing with dates and times
library(sp) # Dealing with spatial data
library(sf) # Also for dealing with spatial data


# 1. LOAD DATA FILES ===========================================================

# We will use the "fread" function of the package "data.table" to load all the
# data files into R. It works in a similar way to the "read.table" functions.
# (note that when dates and times have a recognizable format, this function 
# automatically assigns the corresponding class).

# Load deployment metadata (dates and coordinates of the deployed acoustic
# receivers)
deploy <- fread("./data/deployment_metadata.csv")

head(deploy)
# Check the class and the timezone of deployment times
class(deploy$deployment_time_utc)
tz(deploy$deployment_time_utc)

# Load the metadata of fish tagged with acoustic transmitters
fish_ref <- fread("./data/fish_metadata.csv")
head(fish_ref)

# Load detection data
detect <- fread("./data/detection_data.csv.gz")
head(detect)


# 2. PLOT THE RECEIVER NETWORK =================================================

# Receiver coordinates are in decimal latitude and longitude (geographic 
# coordinates)
deploy[, c("latitude", "longitude")]

# We can define this projection using the proj4 system (library for performing 
# conversions between cartographic projections)
proj_longlat <- CRS("+proj=longlat +datum=WGS84")
# Reference for the projection: https://epsg.io/4326

# For plotting and analysis, we will use a UTM system (projected coordinates)
proj_utm <- CRS("+proj=utm +datum=WGS84 +zone=31")
# Reference: https://epsg.io/32631

# To transform from one coordinate system to another, we have first to create
# a "SpatialPoints object
sp <- SpatialPoints(coords = deploy[, c("longitude", "latitude")], 
                    proj4string = proj_longlat)
coordinates(sp)
plot(sp, axes = T)

# We can transform them using the "spTransform" function and providing the
# desired coordinate system
sp_t <- spTransform(sp, CRSobj = proj_utm)
coordinates(sp_t)
plot(sp_t, axes = T)

# Add the new coordinates to the "deploy" data.table
deploy$x <- coordinates(sp_t)[, 1]
deploy$y <- coordinates(sp_t)[, 2]

# Another advantage of projected systems is that we can calculate distances
# between coordinates directly (Euclidean distance)
dist_mat <- dist(deploy[, c("x", "y")])
dist_mat

# Load coastline shapefile using the "sf" package
coast <- sf::read_sf("./data/coastline_shp/mallorca.shp")

# Transfrom to a sp object (easier to plot)
coast <- as_Spatial(coast)
plot(coast, xlim = range(deploy$x), ylim = range(deploy$y), 
     col = "gray80", border = "gray50", main = "Receiver array")
points(deploy$x, deploy$y, pch = 21, bg = "gold", cex = 1.4)
axis(1)
axis(2)
box()


# 3. ASSIGN STATION AND FISH IDs TO DETECTIONS =================================

# 3.1. Assign station IDs ------------------------------------------------------

# We will assign the station ID to each detection, based on the receiver
# IDs and the time they were in the water (time between "date_in" and 
# "date_out")
head(deploy)

# We will use a "for" loop to iterate over each deployment (each row) and look
# for all the detections that match the same receiver ID and occurred between
# the deployment and retrieval times.

# First, generate an empty variable in data where we will store the station ID
detect$station <- NA
head(detect)

# For loop, iterating over every row in the deployment metadata
for (i in 1:nrow(deploy)) {
  
  # To control the progress, this line will print the "i" value in the console
  cat(i, "\n") 
  
  # Find the detections that match the receiver ID and the deployment and 
  # recovery using logical operators
  indx <- (detect$receiver == deploy$receiver_id[i] & 
             detect$date_time > deploy$deployment_time_utc[i] & 
             detect$date_time < deploy$recovery_time_utc[i])
  
  # Assign the station
  detect$station[indx] <- deploy$station[i]
  
}
head(detect)

# Number of detections at different stations
table(detect$station)


# 3.2. Assign fish IDs ---------------------------------------------------------

# We will assign the fish ID to each detection, based on the tag ID and the
# tagging date
head(fish_ref)

# First, generate a new variable to store the fish id
detect$fish_id <- NA

# For loop, iterating over every row in the fish metadata
for (i in 1:nrow(fish_ref)) {
  
  # To control the progress, this line will print the "i" value in the console
  # and the total number of rows
  cat(i, "/", nrow(fish_ref), "\n", sep = "") 

  # Find the data that match the transmitter ID and the release time
  indx <- (detect$tag_id == fish_ref$tag_id[i] & 
             detect$date_time >= fish_ref$release_time_utc[i])
  
  # Assign ID to the detections
  detect$fish_id[indx] <- fish_ref$fish[i]

}

# Now we have our data ready to be analyzed!
head(detect)

# Number of detections of each tagged individual
table(detect$fish_id)


# 4. EXPORT DATA ===============================================================

# Export data to be used in the next script
save(detect, fish_ref, deploy, coast, file = "./data/data_for_analysis.rda", 
     compress = "xz")

