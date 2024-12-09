## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## 
## INTRODUCTION TO THE ANALYSIS OF ACOUSTIC TELEMETRY DATA IN R
## 2. ANALYZING ACOUSTIC TELEMETRY DATA
##
## Mallorca Science School (MASS) - Palma (Balearic Islands)
##
## Author: Eneko Aspillaga (IMEDEA, CSIC-UIB, Spain)
## Contact: aspillaga@imedea.uib-csic.es
## Date: October 21-25 2024
## 
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Remove all the objects in the environment.
rm(list = ls())

# Set working directory (directory of the current script)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Install required libraries:
# install.packages(c("data.table", "lubridate", "sp", "sf", "igraph"))

# Load libraries
library(data.table) # Extension of data.frame for faster management
library(lubridate) # Dealing with dates and times
library(sp) # Dealing with spatial data
library(sf) # Also for dealing with spatial data
library(igraph) # Network analysis

# Load the data prepared in the previous script
load("./data/data_for_analysis.rda")

# Function to rescale variables
rescale <- function(values, orig_range = NULL, new_range) {
  if (is.null(orig_range)) orig_range <- c(0, max(values, na.rm = TRUE))
  val <- (new_range[1] + (values - orig_range[1]) * 
            (new_range[2] - new_range[1]) / (orig_range[2] - orig_range[1]))
  return(val)
}


# 1. DATA EXPLORATION ==========================================================

# Check detections not assigned to any individual (ghost detections)
false_indx <- is.na(detect$fish_id) | detect$fish_id == ""
false_ids <- table(detect$tag_id[false_indx])
false_ids[order(false_ids, decreasing = TRUE)][1:50]
table(detect$protocol[false_indx])

# We can get rid of those detections
detect <- detect[!false_indx]

# Number of detections of each individual
n_det_ind <- table(detect$fish_id)
n_det_ind


# 2. PLOT NUMBER OF DETECTIONS IN EACH STATION =================================

# Convert station to a factor (so that also stations without detections
# are taken into account)
detect$station <- factor(detect$station, levels = deploy$station)

# Count detections at each station
n_det_st <- table(detect$station)
n_det_st

# Re-scale detection count to specific range (cex size of points)
log_det <- log10(n_det_st + 1)
cex <- rescale(log_det, new_range = c(0.8, 4))

# Define colours for stations. A different colour will be given to stations 
# without detections
col <- ifelse(n_det_st > 0, "gold", "blue")

# Plot map
par(mar = c(3, 3, 1, 1))
plot(coast, col = "gray80", border = "gray50", xlim = range(deploy$x),
     ylim = range(deploy$y))

# Add stations
points(deploy$x, deploy$y, pch = 21, cex = cex, bg = col)

# Add legend (has to be rescaled as well)
legend_values <- c(1, 3, 5)
legend_cex <- rescale(legend_values, orig_range = c(0, max(log_det)),
                      new_range =c(0.8, 4))
legend("bottomleft", legend = 10^legend_values, pch = 21, pt.bg = "gold",
       pt.cex = legend_cex, y.intersp = 2, x.intersp = 2, 
       title = "No. of detections", inset = c(0.02, 0.02))
axis(1)
axis(2)
box()


# 3. RESIDENCY ANALYSIS ========================================================

# 3.1. Calculate the residency-related parameters of each individual -----------

# We will use the special syntax of "data.table" to easily calculate, for each
# individual ("by" argument), the total number of detections ("n_det"), the 
# number of stations with detections ("n_st"), the number of days with
# detections ("days_detected") and the date of the last detection ("last_date")
residency <- detect[, .(n_det = .N,
                        n_st = length(unique(station)),
                        days_detected = length(unique(as.Date(date_time))),
                        last_date = max(as.Date(date_time))),
                    by = fish_id]
head(residency)

# Match this table with the fish reference table
residency <- residency[match(fish_ref$fish_id, fish_id), ]

# Add the tagging date from the fish metadata data frame
residency$tag_date <- as.Date(fish_ref$release_time_utc)
head(residency)

# Set the end of the tracking period at the date at which receivers were 
# retrieved
track_end <- as.Date(max(deploy$recovery_time_utc))
track_end

# Calculate tracking period of each individual as the number of days elapsed
# between the tagging date and the end of the tracking period
residency$tracking_period <- as.numeric(track_end - residency$tag_date + 1)

# Calculate residency index (RI) as the ratio between the days with detections
# and the tracking period
residency$ri <- residency$days_detected / residency$tracking_period
head(residency)

# Compare tracking period of different species
par(mar = c(3.1, 4.8, 1, 1))
boxplot(residency$ri ~ fish_ref$species, col = 2:4, ylab = "Residency index")


# 3.2. Plot days with detections -----------------------------------------------

# Function to identify the station with most detections
maxSt <- function(st) {
  tab <- table(st)
  sample(names(tab[tab == max(tab)]), 1)
}


# Count, for each tracking day and individual ("by" argument), the total number 
# of detections ("n_det"), the station with most detections ("main_st") and the
# number of stations with detections
daily_detect <- detect[, .(n_det = .N,
                           n_st = length(unique(station)),
                           main_st = maxSt(station)),
                       by = list(date = as.Date(date_time), fish_id)]
head(daily_detect)

# Convert fish_id to a factor
daily_detect[, fish_id := factor(fish_id, levels = fish_ref$fish_id)]

# Add species information
daily_detect[, species := fish_ref$species[match(fish_id, fish_ref$fish_id)]]

# Simple plot of detection days (coloured by species)
plot(as.numeric(fish_id) ~ date, data = daily_detect, col = factor(species),
     ylim = c(nrow(fish_ref), 1), pch = "|", ylab = "Fish ID")

# Colour depending on the zone (station codes)
st_zones <- unique(substr(deploy$station, 1, 3))
col <- c("forestgreen", "dodgerblue", "orange", "firebrick")
names(col) <- st_zones

par(mar = c(3.1, 6.2, 1, 1))
plot(as.numeric(fish_id) ~ date, data = daily_detect, type = "n", axes = FALSE,
     ann = FALSE, ylim = c(nrow(fish_ref), 1))
abline(h = 1:nrow(fish_ref), col = "gray60", lty = 2)
points(as.numeric(fish_id) ~ date, data = daily_detect, pch = "|", 
       col = col[substr(main_st, 1, 3)])

# Add tagging date
points(as.Date(fish_ref$release_time_utc), 1:nrow(fish_ref), pch = 16, cex = 0.9)

axis(2, at = 1:nrow(fish_ref), labels = fish_ref$fish_id, las = 1)
axis.Date(1, at = pretty(daily_detect$date, 10))
box()
title(ylab = "Fish ID", line = 4.7)

legend("bottomleft", legend = st_zones, fill = col, inset = c(0.02, 0.02),
       xpd = TRUE)


# 4. NETWORK ANALYSIS ==========================================================

# We will contruct the networks from the daily detection data

# Convert the station variable to factor so that all the stations are
# taken into account
daily_detect$main_st <- factor(daily_detect$main_st, levels = deploy$station)


# 4.1. Transition matrices -----------------------------------------------------

# Generate transition matrices. These square matrices indicate the number of
# movements observed between stations.
tr_mat <- lapply(fish_ref$fish_id, function(ind) {
  
  # Subset data for the individual
  det_sub <- daily_detect[fish_id == ind]
  
  # Count transitions between stations
  trans <- data.frame(from = det_sub$main_st[-nrow(det_sub)],
                      to = det_sub$main_st[-1])
  trans <- as(table(trans$from, trans$to, dnn = c("from", "to")), "matrix")
  return(trans)
})
names(tr_mat) <- fish_ref$fish_id

head(tr_mat[[1]][1:10, 1:10])

# 4.2. Generate networks using igraph ------------------------------------------

net_list <- lapply(seq_along(tr_mat), function(i) {
  
  net <- graph_from_adjacency_matrix(tr_mat[[i]], mode = "directed", 
                                     weighted = TRUE)
  
  # Add coordinate data to nodes
  V(net)$x <- deploy$x[match(V(net)$name, deploy$station)]
  V(net)$y <- deploy$y[match(V(net)$name, deploy$station)]
  
  # Add total number of detections per station
  tot_det <- table(detect$station[detect$fish_id == names(tr_mat)[i]])
  V(net)$n_det <- tot_det[V(net)$name]
  
  return(net)
  
})
names(net_list) <- names(tr_mat)


# 4.3. Plot networks -----------------------------------------------------------

for (i in seq_along(net_list)) {
  
  net <- net_list[[i]]
  
  # Remove nodes without edges and loops
  net <- delete_vertices(net, degree(net) == 0)
  net <- simplify(net, remove.multiple = TRUE, remove.loops = TRUE)
  
  # Set limits for the plot
  xlim <- range(deploy$x[deploy$station %in% V(net)$name])
  ylim <- range(deploy$y[deploy$station %in% V(net)$name])
  xlim <- xlim + c(-1, 1) * diff(xlim) / 10
  ylim <- ylim + c(-1, 1) * diff(ylim) / 10
  
  par(mar = c(1, 1, 3, 1))
  
  plot(deploy$x, deploy$y, asp = 1, xlim = xlim, ylim = ylim, 
       type = "n", ann = FALSE, axes = FALSE)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = "#AAE3F9", border = "transparent")
  
  # Plot coastline
  plot(coast, add = TRUE, col = "gray80", border = "gray30")
  
  # Plot all the deployments
  points(y ~ x, data = deploy, pch = 16, cex = 0.8)
  
  # Add receivers with detections
  points(V(net)$x, V(net)$y, cex = rescale(V(net)$n_det, new_range = c(1.2, 4)),
         pch = 21, bg = "#5ADB5A")
  
  plot(net, rescale = FALSE, add = TRUE,
       edge.curved = 0.5, edge.arrow.size = 0.8,
       edge.color = adjustcolor("#2C1F89", 0.9),
       edge.width = rescale(E(net)$weight, range(E(net)$weight),
                            new_range = c(1.8, 3)),
       vertex.label = "", vertex.color = "transparent",
       vertex.frame.color = "transparent")
  
  title(main = names(net_list)[i])
}


# 4.4. Extract network metrics -------------------------------------------------

# NETWORK-LEVEL METRICS:

# Edge density (observed edges / possible edges)
edge_dens <- sapply(net_list, edge_density)
# Compare the edge density  of different species
par(mar = c(3.1, 4.8, 1, 1))
boxplot(edge_dens ~ fish_ref$species, col = 2:4, ylab = "Residency index")


# NODE-LEVEL METRICS
net <- net_list[["DP0003"]]

# Degree (number of edges of each node
deg <- degree(net)
deg

# Strength (sum of all the edge weights of a node, represents tge association
# rate per node)
str <- strength(net)
str

# Betweenenss centrality (count of the number of shortest paths going through
# each node, indicates how important a node is in connecting different parts of 
# the network)
btw <- betweenness(net)
btw

# Eigenvector centrality (sum of the centralities of the neighbours of a node, 
# indicates the influence of a node in a nerwork, high centralities are reached 
# by having a high degree or by being connected to associates with a high 
# degree)
eig <- eigen_centrality(net)
eig$vector


