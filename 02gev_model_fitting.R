############################################################
#
# This script contains the main real-data analysis used in
# the accompanying report on extreme sea surface temperatures
# (SSTs) in the Red Sea.
#
# The aim of this part of the analysis is to model annual
# extreme SSTs using Generalized Extreme Value (GEV) theory,
# while allowing the model parameters to vary with
# geographical and environmental covariates.
#
# The script includes:
#   1. Import of the Red Sea SST dataset and bathymetry data
#   2. Exploratory data analysis of spatial and temporal SST
#      variation
#   3. Single-location GEV fitting across all 108 locations
#   4. Construction of annual block maxima and alignment with
#      location-based covariates
#   5. Exploration of location, scale, and shape parameters
#   6. Comparison of candidate non-stationary GEV models
#      using AIC
#   7. Selection of a final model for the Red Sea SST data
#
# Notes:
#   - Required external files include:
#       * Red_Sea_data.RData
#       * bathymetry.xlsx
############################################################




############################################################
# Import Real Data
#
############################################################

library(readxl)

# Load Red Sea SST data
load("Red_Sea_data.RData")

# Load bathymetry data
bathymetry <- read_excel("bathymetry.xlsx")

# Store SST data in a separate object
RedSST <- data






############################################################
# EDA
#
# This section carries out exploratory data analysis for the
# Red Sea SST dataset. The aim is to examine the main spatial
# and temporal features of the data before fitting the GEV
# models.
############################################################

library(ggplot2)

# mean SST at each location
mean_by_loc <- colMeans(RedSST)

# organize all variables into one dataframe
mean_SST_loc_df <- data.frame(
  MeanTemp = mean_by_loc,
  Longitude = loc[,1],
  Latitude = loc[,2],
  Depth = bathymetry[,3]
)


############################################################
# Spatial variability
#
# The first plot shows the spatial pattern of mean SST across
# locations in the Red Sea.
############################################################

# spatical variability
ggplot(mean_SST_loc_df, aes(x = Longitude, y = Latitude, color = MeanTemp)) +
  geom_point(size = 1.5) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Spatial Variability in Sea Surface Temperature (SST)",
       x = "Longitude",
       y = "Latitude",
       color = "Mean SST (deg C)")


############################################################
# Temporal variability
#
# The next plots show the daily mean SST across all locations,
# first over the full sample and then over a shorter period
# for closer inspection.
############################################################

# mean SST for each day
mean_by_day <- rowMeans(RedSST)

mean_SST_day_df <- data.frame(
  MeanTemp = mean_by_day,
  DayIndex = 1:5840)

ggplot(mean_SST_day_df, aes(x = DayIndex, y = MeanTemp)) +
  geom_line(color = "blue", size = 1) +
  labs(title = "Mean Daily Sea Surface Temperature (SST)",
       x = "Day Index",
       y = "Mean SST (deg C)") +
  scale_x_continuous(breaks = seq(0, 5840, by = 180)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(mean_SST_day_df[1:730,], aes(x = DayIndex, y = MeanTemp)) +
  geom_line(color = "blue", size = 1) +
  labs(title = "Mean Daily Sea Surface Temperature (SST)",
       x = "Day Index",
       y = "Mean SST (deg C)") +
  scale_x_continuous(breaks = seq(0, 5840, by = 180)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1))


############################################################
# Distribution of SST values
#
# A histogram is used to inspect the overall distribution of
# SST values in the dataset.
############################################################

hist(RedSST, 
     main = "Histogram of Sea Surface Temperature (SST)",
     xlab = "SST (deg C)",
     ylab = "Frequency")


############################################################
# Seasonal separation
#
# This part separates the daily mean SST into summer and
# winter periods in order to compare seasonal behaviour more
# directly.
############################################################
library(dplyr)
library(ggplot2)

# mean daily SST in summer and winter
mean_SST_day_df <- rowMeans(RedSST)
summer <- mean_SST_day_df[month %in% c(6,7,8)]
winter <- mean_SST_day_df[month %in% c(12,1,2)]

# line plot of summer mean SST
y_limits <- range(c(max(summer), min(winter)), na.rm = TRUE)
plot(summer, type = "l", col = "red", 
     xlab = "Day", ylab = "Mean SST (deg C)",
     main = "Mean Daily Sea Surface Temperature (SST) in Summer and Winter", ylim = y_limits)
# add a line for winter mean SST
lines(winter, type = "l", col = "blue")
legend("topright", legend = c("Summer", "Winter"), col = c("red", "blue"), lty = 1)


# zoom in to the first six months of data
first_six_summer <- summer[1:184]
first_six_winter <- winter[1:180]

plot(first_six_summer, type = "l", col = "red", 
     xlab = "Day", ylab = "Mean SST (deg C)",
     main = "Mean Daily Sea Surface Temperature (SST) in Summer and Winter", ylim = y_limits)
# add a line for winter mean SST
lines(first_six_winter, type = "l", col = "blue")
legend("topright", legend = c("Summer", "Winter"), col = c("red", "blue"), lty = 1)






############################################################
# Single Location
#
# The next part fits GEV models separately at each location
# using annual block maxima, in order to obtain location-wise
# parameter estimates that can later be compared with spatial
# covariates.
############################################################


## Function Combining BM and GEV Fit

####################################################################################################
# Description:
# Fits a Generalized Extreme Value (GEV) distribution to annual block maxima from daily SST data.
#
# Parameters:
# SST_data: Numeric vector of daily SST values for a single location
#
# Returns:
# GEV fit results
#
# Process:
# 1. Calculate the number of years from SST_data
# 2. Extract annual maxima and store in `block_maxima`.
# 3. Fit GEV distribution to `block_maxima` using `gev.fit`.
# 4. Return GEV fit results.
####################################################################################################

library(ismev)

fit_gev <- function(SST_data) {
  num_years <- 16
  block_maxima <- numeric(num_years)
  
  # find the maxima for each year
  for (i in 1:num_years) {
    start_index <- 365*i - 364
    end_index <- 365 * i
    block_maxima[i] <- max(SST_data[start_index:end_index])
  }
  
  # fit the GEV distribution to the block maxima (use identity link for all parameters)
  gev_fit <- gev.fit(
    xdat = block_maxima,
    method = "Nelder-Mead",
    maxit = 10000
  )
  
  return(gev_fit)
}



############################################################
# GEV fits for all 108 locations
#
# A separate GEV model is fitted at each location, and the
# maximum likelihood estimates are stored for later spatial
# analysis.
############################################################

num_locations <- 108

# initialize a list to store MLEs for each location
mles <- matrix(NA, nrow = 108, ncol = 3)

# loop through each location to apply the fit_gev function
for (i in 1:108) {
  RedSST_loc <- RedSST[, i]
  capture.output(gev_fit_single <- fit_gev(RedSST_loc))
  
  # store the mle results
  mles[i, ] <- gev_fit_single$mle
}

# create a dataframe to store the mles for 108 locations
mles_data <- as.data.frame(mles)
colnames(mles_data) <- c("Location", "Scale", "Shape")
print(mles_data)




############################################################
# Incorporate Locations
#
# This section prepares the block maxima and the associated
# spatial/environmental covariates so that non-stationary GEV
# models can be fitted across locations.
############################################################

## Data Alignment 

### Block Maxima

# calculate block maxima for each location
block_maxima <- sapply(1:108, function(loc) {
  sapply(1:16, function(year) {
    start_index <- 365*year - 364
    end_index <- year * 365
    max(data[start_index:end_index, loc])
  })
})

block_maxima_df <- as.data.frame(block_maxima)
colnames(block_maxima_df) <- paste("Location", 1:108)
rownames(block_maxima_df) <- paste("Maxima for Year", 1:16)
print(block_maxima_df)




############################################################
# Location data
#
# The corresponding longitude, latitude, depth, and log-depth
# values are aligned with the block maxima for each location.
############################################################

location <- data.frame(loc)

# combine the corresponding location data for each block maxima
location_data <- data.frame(
  longitude = rep(location$lon, each = 16),
  latitude = rep(location$lat, each = 16),
  depth = rep(bathymetry$dep, each = 16),
  log_depth = rep(log(bathymetry$dep), each = 16)
)
print(location_data)

# combine loc and mles
loc_mle <- cbind(loc, mles_data)

# add water depth into loc_mle
loc_mle$dep <- bathymetry$dep

# add log (scale) 
loc_mle$log_dep <- log(loc_mle$dep)



############################################################
# Alignment
#
# The block maxima and covariates are reshaped into formats
# suitable for model fitting with covariates.
############################################################

# sequence: Max for year 1 in location 1; max for year 2 in location 1; ... max for year 16 in location 1; max for year 1 in location 2; ... 
block_maxima_vector <- as.vector(t(t(block_maxima_df)))
print(block_maxima_vector)

# location coordinate for location 1 x 16 times; for location 2 x 16 times
location_data <- as.matrix(location_data)
print(location_data)




############################################################
# Location
#
# This section explores how the estimated location parameter
# may vary with longitude, latitude, water depth, and
# log-depth.
############################################################

## ggplots

library(ggplot2)

par(mfrow = c(3, 1))

# scatter plot for location parameter with both lon and lat
ggplot(loc_mle, aes(x = lon, y = lat, color = Location)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Location Parameter vs Longitude and Latitude", 
       x = "Longitude", 
       y = "Latitude", 
       color = "Location MLEs") 

# against longitude
ggplot(loc_mle, aes(x = lon, y = Location)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") + 
  labs(title = "Location Parameter vs Longitude", x = "Longitude", y = "Location")

# against latitude
ggplot(loc_mle, aes(x = lat, y = Location)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") + 
  labs(title = "Location Parameter vs Latitude", x = "Latitude", y = "Location")

# against depth
ggplot(loc_mle, aes(x = dep, y = Location)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") + 
  labs(title = "Location Parameter vs Water Depth", x = "Water Depth", y = "Location")

# against log_depth
ggplot(loc_mle, aes(x = log_dep, y = Location)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") + 
  labs(title = "Location Parameter vs Log Water Depth", x = "Log Water Depth", y = "Location")






############################################################
# Candidate models for the location parameter
#
# Several non-stationary GEV models are fitted using
# different covariate combinations in the location parameter,
# and AIC is used for comparison.
############################################################


## Ismev   
library(ismev)

# location with both lon and lat
gev_fit1 <- gev.fit(
  xdat = block_maxima_vector,
  ydat = location_data,
  mul = c(1,2),  
  method = "Nelder-Mead",
  maxit = 10000
)


# location with only lon
gev_fit2 <- gev.fit(
  xdat = block_maxima_vector,
  ydat = location_data,
  mul = c(1),  
  method = "Nelder-Mead",
  maxit = 10000
)


# location with only latitude
gev_fit3 <- gev.fit(
  xdat = block_maxima_vector,
  ydat = location_data,
  mul = c(2),  
  method = "Nelder-Mead",
  maxit = 10000
)

gev_fit4 <- gev.fit(
  xdat = block_maxima_vector,
  ydat = location_data,
  mul = c(1,2,4),  
  method = "Nelder-Mead",
  maxit = 10000
)


gev_fit5 <- gev.fit(
  xdat = block_maxima_vector,
  ydat = location_data,
  mul = c(1,4),  
  method = "Nelder-Mead",
  maxit = 10000
)

gev_fit6 <- gev.fit(
  xdat = block_maxima_vector,
  ydat = location_data,
  mul = c(2,4),  
  method = "Nelder-Mead",
  maxit = 10000
)

gev_fit7 <- gev.fit(
  xdat = block_maxima_vector,
  ydat = location_data,
  mul = c(4),  
  method = "Nelder-Mead",
  maxit = 10000
)


AIC_ismev1 <- 2 * length(gev_fit1$mle) + 2 * gev_fit1$nllh
AIC_ismev2 <- 2 * length(gev_fit2$mle) + 2 * gev_fit2$nllh
AIC_ismev3 <- 2 * length(gev_fit3$mle) + 2 * gev_fit3$nllh
AIC_ismev4 <- 2 * length(gev_fit4$mle) + 2 * gev_fit4$nllh
AIC_ismev5 <- 2 * length(gev_fit5$mle) + 2 * gev_fit5$nllh
AIC_ismev6 <- 2 * length(gev_fit6$mle) + 2 * gev_fit6$nllh
AIC_ismev7 <- 2 * length(gev_fit7$mle) + 2 * gev_fit7$nllh

AIC_ismev <- data.frame(
  Model = c("ismev_fit1", "ismev_fit2", "ismev_fit3", "ismev_fit4", "ismev_fit5", "ismev_fit6", "ismev_fit7"),
  AIC = c(AIC_ismev1, AIC_ismev2, AIC_ismev3, AIC_ismev4, AIC_ismev5, AIC_ismev6, AIC_ismev7)
)
print(AIC_ismev)







############################################################
# Scale
#
# This section explores how the estimated scale parameter may
# vary spatially and whether water depth may help explain its
# behaviour.
############################################################

## ggplots

library(ggplot2)
par(mfrow = c(3, 1))
# add log (scale) 
loc_mle$log_Scale <- log(loc_mle$Scale)

ggplot(loc_mle, aes(x = lon, y = lat, color = Scale)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Scale Parameter vs Longitude and Latitude", 
       x = "Longitude", 
       y = "Latitude", 
       color = "Scale") 

ggplot(loc_mle, aes(x = lon, y = log_Scale)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Log Scale Parameter vs Longitude", x = "Longitude", y = "Log Scale")


ggplot(loc_mle, aes(x = lat, y = log_Scale)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") + 
  labs(title = "Log Scale Parameter vs Latitude", x = "Latitude", y = "Log Scale")


############################################################
# Bathymetry
#
# These plots investigate the relationship between the scale
# parameter and water depth.
############################################################

# plot of water depth at each location
ggplot(bathymetry, aes(x = lon, y = lat, color = dep)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Water Depth at Each Location", 
       x = "Longitude", 
       y = "Latitude", 
       color = "Depth (m)")

# plot of the relationship of scale parameter and the water depth
ggplot(loc_mle, aes(x = dep, y = log_Scale)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Log Scale Parameter vs Water Depth", x = "Water Depth", y = "Log Scale")

# plot of the relationship of scale parameter and the log water depth
ggplot(loc_mle, aes(x = log_dep, y = log_Scale)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Log Scale Parameter vs Log Water Depth", x = "Log Water Depth", y = "Log Scale")



############################################################
# Candidate models for the scale parameter
#
# Models are fitted using the texmex package to allow
# different covariate combinations in the scale parameter,
# and AIC is again used for comparison.
############################################################

library(texmex)

# combine block maxima with the location data and water depth
location_data_df <- as.data.frame(cbind(BlockMaxima = block_maxima_vector, location_data))

# create gev models using texmex

# scale with both lon and lat
evm_fit1 <- evm(BlockMaxima, data=location_data_df, family=gev, 
                mu = ~ longitude + latitude + log_depth,
                phi = ~ longitude + latitude)
print(evm_fit1)

# scale with only lon
evm_fit2 <- evm(BlockMaxima, data=location_data_df, family=gev, 
                mu = ~ longitude + latitude + log_depth,
                phi = ~ longitude)
print(evm_fit2)

# scale with only lat
evm_fit3 <- evm(BlockMaxima, data=location_data_df, family=gev, 
                mu = ~ longitude + latitude + log_depth,
                phi = ~ latitude)
print(evm_fit3)

# scale with lon, lat and log_depth
evm_fit4 <- evm(BlockMaxima, data=location_data_df, family=gev, 
                mu = ~ longitude + latitude + log_depth,
                phi = ~ longitude + latitude + log_depth)
print(evm_fit4)


evm_fit5 <- evm(BlockMaxima, data=location_data_df, family=gev, 
                mu = ~ longitude + latitude + log_depth,
                phi = ~ longitude + log_depth)
print(evm_fit5)

evm_fit6 <- evm(BlockMaxima, data=location_data_df, family=gev, 
                mu = ~ longitude + latitude + log_depth,
                phi = ~ latitude + log_depth)
print(evm_fit6)


evm_fit7 <- evm(BlockMaxima, data=location_data_df, family=gev, 
                mu = ~ longitude + latitude + log_depth,
                phi = ~ log_depth)
print(evm_fit7)


AIC_evm1 <- AIC(evm_fit1)
AIC_evm2 <- AIC(evm_fit2)
AIC_evm3 <- AIC(evm_fit3)
AIC_evm4 <- AIC(evm_fit4)
AIC_evm5 <- AIC(evm_fit5)
AIC_evm6 <- AIC(evm_fit6)
AIC_evm7 <- AIC(evm_fit7)

# create a data frame to store the AIC values
AIC_evm <- data.frame(
  Model = c("evm_fit1", "evm_fit2", "evm_fit3", "evm_fit4", "evm_fit5", "evm_fit6", "evm_fit7"),
  AIC = c(AIC_evm1, AIC_evm2, AIC_evm3, AIC_evm4, AIC_evm5, AIC_evm6, AIC_evm7)
)

print(AIC_evm)





############################################################
# Combined AIC comparison
#
# The AIC values from the ismev and texmex models are brought
# together to compare all candidate specifications.
############################################################

## All AIC

AIC_results_combined <- rbind(AIC_ismev, AIC_evm)

# Order the combined results from smallest to largest AIC value
AIC_results_combined <- AIC_results_combined[order(AIC_results_combined$AIC), ]

# Print the combined and ordered AIC results
print(AIC_results_combined)


evm_fit5 <- evm(BlockMaxima, data=location_data_df, family=gev, 
                mu = ~ longitude + latitude + log_depth,
                phi = ~ longitude + log_depth)
print(evm_fit3)





############################################################
# Shape
#
# The shape parameter is also inspected visually, although in
# environmental applications it is often assumed to be common
# across locations.
############################################################

## ggplots

par(mfrow = c(3, 1))
ggplot(loc_mle, aes(x = lon, y = lat, color = Shape)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Shape Parameter vs Longitude and Latitude", 
       x = "Longitude", 
       y = "Latitude", 
       color = "Shape") 

ggplot(loc_mle, aes(x = lon, y = Shape)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") + # Linear trend line
  labs(title = "Shape Parameter vs Longitude", x = "Longitude", y = "Shape")


ggplot(loc_mle, aes(x = lat, y = Shape)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") + # Linear trend line
  labs(title = "Shape Parameter vs Latitude", x = "Latitude", y = "Shape")





############################################################
# Final Model
#
# A preferred model is then fitted using the selected
# covariates for the location and scale parameters.
############################################################

## final model1

gev_final1 <- (evm(BlockMaxima, data=location_data_df, family=gev, 
                   mu = ~ longitude + latitude + log_depth,
                   phi = ~ longitude + log_depth))
print(gev_final1)




############################################################
# Multicollinearity check
#
# This part introduces distance from a reference point as an
# alternative covariate and explores whether it may reduce
# issues of collinearity among the explanatory variables.
############################################################

# define the coordinates of the reference location
location1 <- c(34.75,28.75)
location_matrix <- as.matrix(location)

# compute distances from the reference to all other locations
distances <- dist(rbind(location1, location_matrix))

distances_matrix <- as.matrix(distances)

# extract distances from reference
distances_to_L1 <- distances_matrix[1, -1] 

# add the column of distances to the dataframe
location_data <- as.data.frame(location_data)
location_data$distance <- rep(distances_to_L1, each = 16)
location_data$log_distance <- log(location_data$distance)
print(location_data)

ggplot(location_data, aes(x = longitude, y = latitude, color = distance)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Distance from each Coordinate to the Reference Point", 
       x = "Longitude", 
       y = "Latitude", 
       color = "Distance") 

ggplot(location_data, aes(x = longitude, y = latitude, color = log_distance)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Log Distance from each Coordinate to the Reference Point", 
       x = "Longitude", 
       y = "Latitude", 
       color = "Log Distance") 

loc_mle$log_distance <- log(distances_to_L1)
loc_mle$distance <- distances_to_L1

ggplot(loc_mle, aes(x = distance, y = log_Scale)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") + 
  labs(title = "Log Scale Parameter vs Latitude", x = "Distance", y = "Log Scale")

ggplot(loc_mle, aes(x = log_distance, y = log_Scale)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") + 
  labs(title = "Log Scale Parameter vs Latitude", x = "Log_Distance", y = "Log Scale")

ggplot(loc_mle, aes(x = distance, y = Location)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") + 
  labs(title = "Log Scale Parameter vs Latitude", x = "Distance", y = "Location")

ggplot(loc_mle, aes(x = log_distance, y = Location)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") + 
  labs(title = "Log Scale Parameter vs Latitude", x = "Log_Distance", y = "Location")

ggplot(loc_mle, aes(x = dep, y = Location)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") + 
  labs(title = "Log Scale Parameter vs Latitude", x = "Depth", y = "Location")

ggplot(loc_mle, aes(x = log_dep, y = Location)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") + 
  labs(title = "Log Scale Parameter vs Latitude", x = "Log_depth", y = "Location")



############################################################
# Alternative models with distance covariates
#
# Additional models are fitted using distance and log-depth
# covariates in order to check whether these lead to a better
# final specification.
############################################################

## ismev
library(ismev)
location_data <- as.matrix(location_data)
# location: with both lon and lat
gev_fit8 <- gev.fit(
  xdat = block_maxima_vector,
  ydat = location_data,
  mul = c(5),  
  method = "Nelder-Mead",
  maxit = 10000
)

gev_fit9 <- gev.fit(
  xdat = block_maxima_vector,
  ydat = location_data,
  mul = c(3,5),  
  method = "Nelder-Mead",
  maxit = 10000
)


gev_fit10 <- gev.fit(
  xdat = block_maxima_vector,
  ydat = location_data,
  mul = c(3),  
  method = "Nelder-Mead",
  maxit = 10000
)

gev_fit11 <- gev.fit(
  xdat = block_maxima_vector,
  ydat = location_data,
  mul = c(4),  
  method = "Nelder-Mead",
  maxit = 10000
)

gev_fit12 <- gev.fit(
  xdat = block_maxima_vector,
  ydat = location_data,
  mul = c(4,6),  
  method = "Nelder-Mead",
  maxit = 10000
)

AIC_ismev8 <- 2 * length(gev_fit8$mle) + 2 * gev_fit8$nllh
AIC_ismev9 <- 2 * length(gev_fit9$mle) + 2 * gev_fit9$nllh
AIC_ismev10 <- 2 * length(gev_fit10$mle) + 2 * gev_fit10$nllh
AIC_ismev11 <- 2 * length(gev_fit11$mle) + 2 * gev_fit11$nllh
AIC_ismev12 <- 2 * length(gev_fit12$mle) + 2 * gev_fit12$nllh


AIC_ismev_new <- data.frame(
  Model = c("ismev_fit8", "ismev_fit9", "ismev_fit10", "ismev_fit11", "ismev_fit12"),
  AIC = c(AIC_ismev8, AIC_ismev9, AIC_ismev10, AIC_ismev11, AIC_ismev12))
print(AIC_ismev_new)



## texmex

library(texmex)

# combine block maxima with the location data and water depth
location_data_df <- as.data.frame(cbind(BlockMaxima = block_maxima_vector, location_data))


# create gev models using texmex
evm_fit8 <- evm(BlockMaxima, data=location_data_df, family=gev, 
                mu = ~ log_distance + log_depth,
                phi = ~ log_distance + log_depth)
print(evm_fit8)


evm_fit9 <- evm(BlockMaxima, data=location_data_df, family=gev, 
                mu = ~ log_distance + log_depth,
                phi = ~ log_distance)
print(evm_fit9)


evm_fit10 <- evm(BlockMaxima, data=location_data_df, family=gev, 
                 mu = ~ log_distance + log_depth,
                 phi = ~ log_depth)
print(evm_fit10)


evm_fit11 <- evm(BlockMaxima, data=location_data_df, family=gev, 
                 mu = ~ log_distance,
                 phi = ~ log_distance + log_depth)
print(evm_fit11)


evm_fit12 <- evm(BlockMaxima, data=location_data_df, family=gev, 
                 mu = ~ log_distance,
                 phi = ~ log_distance)
print(evm_fit12)


evm_fit13 <- evm(BlockMaxima, data=location_data_df, family=gev, 
                 mu = ~ log_distance,
                 phi = ~ log_depth)
print(evm_fit13)


evm_fit14 <- evm(BlockMaxima, data=location_data_df, family=gev, 
                 mu = ~ log_depth,
                 phi = ~ log_distance + log_depth)
print(evm_fit14)


evm_fit15 <- evm(BlockMaxima, data=location_data_df, family=gev, 
                 mu = ~ log_depth,
                 phi = ~ log_distance)
print(evm_fit15)


evm_fit16 <- evm(BlockMaxima, data=location_data_df, family=gev, 
                 mu = ~ log_depth,
                 phi = ~ log_depth)
print(evm_fit16)



AIC_evm8 <- AIC(evm_fit8)
AIC_evm9 <- AIC(evm_fit9)
AIC_evm10 <- AIC(evm_fit10)
AIC_evm11 <- AIC(evm_fit11)
AIC_evm12 <- AIC(evm_fit12)
AIC_evm13 <- AIC(evm_fit13)
AIC_evm14 <- AIC(evm_fit14)
AIC_evm15 <- AIC(evm_fit15)
AIC_evm16 <- AIC(evm_fit16)

# create a data frame to store the AIC values
AIC_evm_new <- data.frame(
  Model = c("evm_fit8", "evm_fit9", "evm_fit10","evm_fit11", "evm_fit12", "evm_fit13","evm_fit14", "evm_fit15", "evm_fit16"),
  AIC = c(AIC_evm8, AIC_evm9, AIC_evm10, AIC_evm11, AIC_evm12, AIC_evm13, AIC_evm14, AIC_evm15, AIC_evm16)
)

print(AIC_evm_new)




############################################################
# Final model with distance covariate
#
# This final specification uses distance and log-depth
# covariates and is kept as an alternative preferred model.
############################################################

## final model 2

library(texmex)
location_data_df <- as.data.frame(cbind(BlockMaxima = block_maxima_vector, location_data))

final_model2 <- evm(BlockMaxima, data=location_data_df, family=gev, 
                    mu = ~ log_distance,
                    phi = ~ log_distance + log_depth)
print(final_model2)


