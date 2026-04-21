############################################################
# Diagnostics
#
# This section evaluates the fitted GEV model using several
# graphical diagnostics. Selected locations are highlighted
# across the Red Sea in order to compare model fit in
# different parts of the study region.
#
# The diagnostics include:
#   - Q-Q plots
#   - density plots
#   - return level plots
#   - exploratory multivariate scatter plots
############################################################

############################################################
# Diagnostics Script Setup
#
# The objects below are loaded or reconstructed so that the 
# diagnostics can be run independently from the earlier scripts.
############################################################

# Required packages
library(readxl)
library(ggplot2)
library(ismev)
library(texmex)
library(evd)
library(dplyr)

############################################################
# Load data
############################################################

# Load Red Sea SST data
load("Red_Sea_data.RData")
# Load bathymetry data
bathymetry <- read_excel("bathymetry.xlsx")
# Store SST data in a separate object
RedSST <- data

############################################################
# Recreate location-level GEV fits
#
# These objects are needed for the later diagnostic and
# return-level calculations.
############################################################

# Function to fit a GEV model to annual block maxima
fit_gev <- function(SST_data) {
  num_years <- 16
  block_maxima <- numeric(num_years)
  
  for (i in 1:num_years) {
    start_index <- 365 * i - 364
    end_index <- 365 * i
    block_maxima[i] <- max(SST_data[start_index:end_index])
  }
  
  gev_fit <- gev.fit(
    xdat = block_maxima,
    method = "Nelder-Mead",
    maxit = 10000
  )
  
  return(gev_fit)
}

# Fit separate GEV models at each location
num_locations <- 108
mles <- matrix(NA, nrow = 108, ncol = 3)

for (i in 1:108) {
  RedSST_loc <- RedSST[, i]
  capture.output(gev_fit_single <- fit_gev(RedSST_loc))
  mles[i, ] <- gev_fit_single$mle
}

mles_data <- as.data.frame(mles)
colnames(mles_data) <- c("Location", "Scale", "Shape")


############################################################
# Recreate block maxima and aligned covariates
############################################################

# Annual block maxima for each location
block_maxima <- sapply(1:108, function(loc) {
  sapply(1:16, function(year) {
    start_index <- 365 * year - 364
    end_index <- year * 365
    max(data[start_index:end_index, loc])
  })
})

block_maxima_df <- as.data.frame(block_maxima)
colnames(block_maxima_df) <- paste("Location", 1:108)
rownames(block_maxima_df) <- paste("Maxima for Year", 1:16)

# Location information
location <- data.frame(loc)

location_data <- data.frame(
  longitude = rep(location$lon, each = 16),
  latitude = rep(location$lat, each = 16),
  depth = rep(bathymetry$dep, each = 16),
  log_depth = rep(log(bathymetry$dep), each = 16)
)

# Combine location information with MLEs
loc_mle <- cbind(loc, mles_data)
loc_mle$dep <- bathymetry$dep
loc_mle$log_dep <- log(loc_mle$dep)

# Vectorised block maxima for model fitting
block_maxima_vector <- as.vector(t(t(block_maxima_df)))


############################################################
# Recreate distance covariate used in the final model
############################################################

# Define the reference location
location1 <- c(34.75, 28.75)
location_matrix <- as.matrix(location)

# Distances from the reference point
distances <- dist(rbind(location1, location_matrix))
distances_matrix <- as.matrix(distances)
distances_to_L1 <- distances_matrix[1, -1]

# Add distance covariates
location_data$distance <- rep(distances_to_L1, each = 16)
location_data$log_distance <- log(location_data$distance)

loc_mle$distance <- distances_to_L1
loc_mle$log_distance <- log(distances_to_L1)


############################################################
# Recreate final fitted model
############################################################

# Prepare data frame for texmex model fitting
location_data_df <- as.data.frame(cbind(BlockMaxima = block_maxima_vector, location_data))

# Final model used for diagnostics
final_model2 <- evm(BlockMaxima, data = location_data_df, family = gev, 
                    mu = ~ log_distance,
                    phi = ~ log_distance + log_depth)

############################################################
# Recreate fitted parameter values used later in diagnostics
############################################################

sigma <- exp(0.6909 + (-0.2856 * loc_mle$log_distance) + (-0.07912 * loc_mle$log_dep))
xi <- -0.366
mu <- 28.11 + (1.711 * loc_mle$log_distance)








############################################################
# Quantile Plot
#
# The first step highlights a small set of representative
# locations. These locations are then used for diagnostic
# comparisons between empirical block maxima and the fitted
# GEV distributions.
############################################################

## Quantile Plot

library(ggplot2)

# mean SST at each location
mean_by_loc <- colMeans(RedSST)

# specify the row numbers to highlight
highlight_rows <- c(1, 22, 44, 70, 90, 108)

# Create a data frame for plotting
mean_SST_loc_df <- data.frame(
  MeanTemp = mean_by_loc,
  Longitude = location[,1],
  Latitude = location[,2],
  Depth = bathymetry[,3],
  Highlight = ifelse(seq_len(nrow(location)) %in% highlight_rows, "Highlight", "Normal")
)

# Plotting spatial variability
ggplot(mean_SST_loc_df, aes(x = Longitude, y = Latitude, color = MeanTemp)) +
  geom_point(size = 1.5) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Spatial Variability in Sea Surface Temperature (SST)",
       x = "Longitude",
       y = "Latitude",
       color = "Mean SST (¡ãC)") +
  geom_point(data = mean_SST_loc_df[highlight_rows, ], 
             aes(x = Longitude, y = Latitude), 
             color = "black", 
             size = 3, 
             shape = 21, 
             fill = "yellow") +  # Customize the highlight color and shape
  theme_minimal()



############################################################
# Example: Q-Q plot for a single location
#
# A Q-Q plot is first produced for Location 1 in order to
# compare empirical quantiles of the observed block maxima
# with the corresponding theoretical quantiles from the
# fitted GEV model.
############################################################

### loc1

sigma <- exp(0.6909 + (-0.2856 * loc_mle$log_distance) + (-0.07912 * loc_mle$log_dep))
xi <- -0.366
mu <- 28.11 + (1.711 * loc_mle$log_distance)

# Coefficient estimates for location 1
sigma_loc1 <- sigma[1]  
xi_loc1 <- -0.366
mu_loc1 <- mu[1]        

library(evd)  

# the probability for each quantile
probability_values <- 1:16 / 17
# sort block maxima for location 1
sorted_bm_loc1 <- sort(block_maxima_df[,1]) 

# empirical quantile
empirical_quantiles_loc1 <- quantile(sorted_bm_loc1, probs = probability_values)

# theoretical quantiles using GEV parameters
theoretical_quantiles_loc1 <- qgev(probability_values, mu_loc1, sigma_loc1, xi_loc1)

# define the axis range
axis_limits <- range(c(theoretical_quantiles_loc1, empirical_quantiles_loc1))

# Q-Q plot
plot(theoretical_quantiles_loc1, empirical_quantiles_loc1,
     main = "Q-Q Plot for Location 1",
     xlab = "Theoretical Quantiles (Model)",
     ylab = "Empirical Quantiles",
     pch = 1,  
     xlim = axis_limits, 
     ylim = axis_limits) 
abline(0, 1, col = "blue")  # diagonal line for reference



############################################################
# Q-Q plots for all selected locations
#
# The same diagnostic is then repeated for the selected
# locations in order to examine whether the fitted GEV model
# provides a reasonable description across space.
############################################################

### All

par(mfrow = c(2, 3))

# Coefficient estimates directly from your provided data
sigma_loc1 <- sigma[1]  # assuming 'sigma' vector is defined in your workspace
xi_loc1 <- -0.366
mu_loc1 <- mu[1]        # assuming 'mu' vector is defined in your workspace

library(evd)  # Assuming 'evd' package is used for qgev

# Assume sigma_loc1, xi_loc1, and mu_loc1 are already calculated or provided

# Calculate the probabilities
probability_values <- 1:16 / 17

sorted_bm_loc1 <- sort(block_maxima_df[,1]) 
# Assuming 'sorted_bm_loc1' is your sorted block maxima data for location 1
# You need to have this data sorted and loaded into your environment
empirical_quantiles_loc1 <- quantile(sorted_bm_loc1, probs = probability_values)

# Calculate the theoretical quantiles using GEV parameters
theoretical_quantiles_loc1 <- qgev(probability_values, mu_loc1, sigma_loc1, xi_loc1)

# Define the axis limits based on the quantiles
axis_limits <- range(c(theoretical_quantiles_loc1, empirical_quantiles_loc1))

# Plotting the Q-Q plot
plot(theoretical_quantiles_loc1, empirical_quantiles_loc1,
     main = "Q-Q Plot for Location 1",
     xlab = "Theoretical Quantiles (Model)",
     ylab = "Empirical Quantiles",
     pch = 1,  # Using a solid circle for better visibility
     xlim = axis_limits, # Set x-axis limits
     ylim = axis_limits) # Set y-axis limits
abline(0, 1, col = "blue")  # diagonal line for reference




# Coefficient estimates directly from your provided data for Location 22
sigma_loc22 <- sigma[22]  # assuming 'sigma' vector is defined in your workspace
xi_loc22 <- -0.366        # same xi as before, assuming it remains constant across locations
mu_loc22 <- mu[22]        # assuming 'mu' vector is defined in your workspace

# Calculate the probabilities
probability_values <- 1:16 / 17

sorted_bm_loc22 <- sort(block_maxima_df[,22]) 
# Assuming 'sorted_bm_loc22' is your sorted block maxima data for location 22
empirical_quantiles_loc22 <- quantile(sorted_bm_loc22, probs = probability_values)

# Calculate the theoretical quantiles using GEV parameters
theoretical_quantiles_loc22 <- qgev(probability_values, mu_loc22, sigma_loc22, xi_loc22)

# Define the axis limits based on the quantiles
axis_limits <- range(c(theoretical_quantiles_loc22, empirical_quantiles_loc22))

# Plotting the Q-Q plot for Location 22
plot(theoretical_quantiles_loc22, empirical_quantiles_loc22,
     main = "Q-Q Plot for Location 22",
     xlab = "Theoretical Quantiles (Model)",
     ylab = "Empirical Quantiles",
     pch = 1,  # Using a solid circle for better visibility
     xlim = axis_limits, # Set x-axis limits
     ylim = axis_limits) # Set y-axis limits
abline(0, 1, col = "blue")  # diagonal line for reference



# Coefficient estimates for Location 44
sigma_loc44 <- sigma[44]  # 'sigma' vector must be defined in your workspace
xi_loc44 <- -0.366        # Assuming xi is consistent across locations
mu_loc44 <- mu[44]        # 'mu' vector must be defined in your workspace

# Calculate probabilities for the quantiles
probability_values <- 1:16 / 17

# Sorted block maxima data for location 44
sorted_bm_loc44 <- sort(block_maxima_df[,44]) 
# Calculate empirical quantiles
empirical_quantiles_loc44 <- quantile(sorted_bm_loc44, probs = probability_values)

# Calculate theoretical quantiles using the GEV model
theoretical_quantiles_loc44 <- qgev(probability_values, mu_loc44, sigma_loc44, xi_loc44)

# Define the axis limits based on the quantiles
axis_limits <- range(c(theoretical_quantiles_loc44, empirical_quantiles_loc44))

# Plotting the Q-Q plot for Location 44
plot(theoretical_quantiles_loc44, empirical_quantiles_loc44,
     main = "Q-Q Plot for Location 44",
     xlab = "Theoretical Quantiles (Model)",
     ylab = "Empirical Quantiles",
     pch = 1,  # Using a solid circle for better visibility
     xlim = axis_limits, # Set x-axis limits
     ylim = axis_limits) # Set y-axis limits
abline(0, 1, col = "blue")  # diagonal line for reference



sigma_loc70 <- sigma[70]  # 'sigma' vector must be defined in your workspace
xi_loc70 <- -0.366        # Assuming xi is consistent across locations
mu_loc70 <- mu[70]        # 'mu' vector must be defined in your workspace

# Calculate probabilities for the quantiles
probability_values <- 1:16 / 17

# Sorted block maxima data for location 70
sorted_bm_loc70 <- sort(block_maxima_df[,70]) 
# Calculate empirical quantiles
empirical_quantiles_loc70 <- quantile(sorted_bm_loc70, probs = probability_values)

# Calculate theoretical quantiles using the GEV model
theoretical_quantiles_loc70 <- qgev(probability_values, mu_loc70, sigma_loc70, xi_loc70)

# Define the axis limits based on the quantiles
axis_limits <- range(c(theoretical_quantiles_loc70, empirical_quantiles_loc70))

# Plotting the Q-Q plot for Location 70
plot(theoretical_quantiles_loc70, empirical_quantiles_loc70,
     main = "Q-Q Plot for Location 70",
     xlab = "Theoretical Quantiles (Model)",
     ylab = "Empirical Quantiles",
     pch = 1,
     xlim = axis_limits,
     ylim = axis_limits)
abline(0, 1, col = "blue")  # diagonal line for reference




sigma_loc90 <- sigma[90]  # 'sigma' vector must be defined in your workspace
xi_loc90 <- -0.366        # Assuming xi is consistent across locations
mu_loc90 <- mu[90]        # 'mu' vector must be defined in your workspace

# Sorted block maxima data for location 90
sorted_bm_loc90 <- sort(block_maxima_df[,90]) 
# Calculate empirical quantiles
empirical_quantiles_loc90 <- quantile(sorted_bm_loc90, probs = probability_values)

# Calculate theoretical quantiles using the GEV model
theoretical_quantiles_loc90 <- qgev(probability_values, mu_loc90, sigma_loc90, xi_loc90)

# Plotting the Q-Q plot for Location 90
plot(theoretical_quantiles_loc90, empirical_quantiles_loc90,
     main = "Q-Q Plot for Location 90",
     xlab = "Theoretical Quantiles (Model)",
     ylab = "Empirical Quantiles",
     pch = 1,
     xlim = axis_limits,
     ylim = axis_limits)
abline(0, 1, col = "blue")




sigma_loc108 <- sigma[108]  # 'sigma' vector must be defined in your workspace
xi_loc108 <- -0.366        # Assuming xi is consistent across locations
mu_loc108 <- mu[108]

# Calculate the probabilities
probability_values <- 1:16 / 17

sorted_bm_loc108 <- sort(block_maxima_df[,108]) 

# Empirical quantiles for location 108
empirical_quantiles_loc108 <- quantile(sorted_bm_loc108, probs = probability_values)

# Theoretical quantiles for location 108
# Note: Make sure mu_loc108, phi_loc108, and xi_loc108 are calculated correctly beforehand
theoretical_quantiles_loc108 <- qgev(probability_values, mu_loc108, sigma_loc108, xi_loc108)

# Define the axis limits based on the quantiles
axis_limits <- range(c(theoretical_quantiles_loc108, empirical_quantiles_loc108))

# Plotting the Q-Q plot
plot(theoretical_quantiles_loc108, empirical_quantiles_loc108,
     main = "Q-Q Plot for Location 108",
     xlab = "Theoretical Quantiles (Model)",
     ylab = "Empirical Quantiles",
     pch = 1,  # Using a solid circle for better visibility
     xlim = axis_limits, # Set x-axis limits
     ylim = axis_limits) # Set y-axis limits
abline(0, 1, col = "blue")  # diagonal line for reference, changed to red for better visibility






############################################################
# Density plot
#
# Density plots are used as a second diagnostic tool. For
# each selected location, the empirical distribution of the
# block maxima is compared with the fitted GEV density.
############################################################

### All

par(mfrow = c(2, 3))

# generate a sequence of 100 equally spaced values 
x_vals_loc1 <- seq(min(sorted_bm_loc1), max(sorted_bm_loc1), length.out = 100)
# GEV density for location 1
gev_density_loc1 <- dgev(x_vals_loc1, mu_loc1, sigma_loc1, xi_loc1)

# empirical density histogram
hist(sorted_bm_loc1, 
     probability = TRUE,  
     breaks = 10,         
     main = "Density Plot for Location 1", 
     xlab = "Block Maxima", 
     ylab = "Density",
     col = "lightgray",   
     border = "black")    

# theoretical density curve
lines(x_vals_loc1, gev_density_loc1, col = "blue", lwd = 2)

# empirical density curve
lines(density(sorted_bm_loc1), col = "black", lwd = 2)

# legend to the plot
legend("topright", legend = c("Empirical", "Model"), col = c("black", "blue"), lwd = 2)




# GEV density for location 22
x_vals_loc22 <- seq(min(sorted_bm_loc22), max(sorted_bm_loc22), length.out = 100)
gev_density_loc22 <- dgev(x_vals_loc22, mu_loc22, sigma_loc22, xi_loc22)

# empirical density histogram
hist(sorted_bm_loc22, 
     probability = TRUE,  
     breaks = 10,         
     main = "Density Plot for Location 22", 
     xlab = "Block Maxima", 
     ylab = "Density",
     col = "lightgray",   
     border = "black")    

# theoretical density curve
lines(x_vals_loc22, gev_density_loc22, col = "blue", lwd = 2)

# empirical density curve
lines(density(sorted_bm_loc22), col = "black", lwd = 2)

# legend to the plot
legend("topright", legend = c("Empirical", "Model"), col = c("black", "blue"), lwd = 2)




# GEV density for location 44
x_vals_loc44 <- seq(min(sorted_bm_loc44), max(sorted_bm_loc44), length.out = 100)
gev_density_loc44 <- dgev(x_vals_loc44, mu_loc44, sigma_loc44, xi_loc44)

# empirical density histogram
hist(sorted_bm_loc44, 
     probability = TRUE,  
     breaks = 10,         
     main = "Density Plot for Location 44", 
     xlab = "Block Maxima", 
     ylab = "Density",
     col = "lightgray",   
     border = "black")    

# theoretical density curve
lines(x_vals_loc44, gev_density_loc44, col = "blue", lwd = 2)

# empirical density curve
lines(density(sorted_bm_loc44), col = "black", lwd = 2)

# legend to the plot
legend("topright", legend = c("Empirical", "Model"), col = c("black", "blue"), lwd = 2)




# GEV density for location 70
x_vals_loc70 <- seq(min(sorted_bm_loc70), max(sorted_bm_loc70), length.out = 100)
gev_density_loc70 <- dgev(x_vals_loc70, mu_loc70, sigma_loc70, xi_loc70)

# empirical density histogram
hist(sorted_bm_loc70, 
     probability = TRUE,  
     breaks = 10,         
     main = "Density Plot for Location 70", 
     xlab = "Block Maxima", 
     ylab = "Density",
     col = "lightgray",   
     border = "black")    

# theoretical density curve
lines(x_vals_loc70, gev_density_loc70, col = "blue", lwd = 2)

# empirical density curve
lines(density(sorted_bm_loc70), col = "black", lwd = 2)

# legend to the plot
legend("topright", legend = c("Empirical", "Model"), col = c("black", "blue"), lwd = 2)




# GEV density for location 90
x_vals_loc90 <- seq(min(sorted_bm_loc90), max(sorted_bm_loc90), length.out = 100)
gev_density_loc90 <- dgev(x_vals_loc90, mu_loc90, sigma_loc90, xi_loc90)

# empirical density histogram
hist(sorted_bm_loc90, 
     probability = TRUE,  
     breaks = 10,         
     main = "Density Plot for Location 90", 
     xlab = "Block Maxima", 
     ylab = "Density",
     col = "lightgray",   
     border = "black")    

# theoretical density curve
lines(x_vals_loc90, gev_density_loc90, col = "blue", lwd = 2)

# empirical density curve
lines(density(sorted_bm_loc90), col = "black", lwd = 2)

# legend to the plot
legend("topright", legend = c("Empirical", "Model"), col = c("black", "blue"), lwd = 2)




# GEV density for location 108
x_vals_loc108 <- seq(min(sorted_bm_loc108), max(sorted_bm_loc108), length.out = 100)
gev_density_loc108 <- dgev(x_vals_loc108, mu_loc108, sigma_loc108, xi_loc108)

# empirical density histogram
hist(sorted_bm_loc108, 
     probability = TRUE,  
     breaks = 10,         
     main = "Density Plot for Location 108", 
     xlab = "Block Maxima", 
     ylab = "Density",
     col = "lightgray",   
     border = "black")    

# theoretical density curve
lines(x_vals_loc108, gev_density_loc108, col = "blue", lwd = 2)

# empirical density curve
lines(density(sorted_bm_loc108), col = "black", lwd = 2)

# legend to the plot
legend("topright", legend = c("Empirical", "Model"), col = c("black", "blue"), lwd = 2)






############################################################
# Return Level
#
# The next part calculates estimated return levels from the
# fitted GEV model and compares them across locations.
############################################################


# function for calculating return levels for 1 year block size
gev_rl <- function(return_period) {
  
  # calculate the probability corresponding to the return period
  p <- 1 - 1 / return_period
  
  # the estimated coefficients
  sigma <- exp(0.6909 + (-0.2856 * loc_mle$log_distance) + (-0.07912 * loc_mle$log_dep))
  xi <- -0.366
  mu <- 28.11 + (1.711 * loc_mle$log_distance)
  
  # calculate the corresponding quantile of the gev distribution
  return_level <- evd::qgev(p, mu, sigma, xi)
  
  return(return_level)
}


############################################################
# Return levels across all locations
#
# Return levels are first calculated for several return
# periods and stored for all locations.
############################################################

## rls

return_periods <- c(10, 20, 50, 100, 200, 500)

# initialize a return matrix
return_levels <- matrix(NA, nrow = length(return_periods), ncol = nrow(loc_mle))

# calculate return levels for each return periods
for (i in seq_along(return_periods)) {
  return_levels[i, ] <- round(gev_rl(return_periods[i]),3)
}

# set rownames as the return period and columnnames as location
rownames(return_levels) <- paste(return_periods, "years")
colnames(return_levels) <- paste("Location", 1:nrow(loc_mle))

return_levels <- as.data.frame(return_levels)
print(return_levels)

print(return_levels[,c(1,22,44,70,90,108)])



############################################################
# Return level plots for selected locations
#
# Return level curves are plotted for the same selected
# locations, together with approximate confidence bands and
# empirical return levels based on block maxima.
############################################################

## All

par(mfrow = c(2, 3))

return_periods_continuous <- seq(10, 500, by = 1)

# calculate return levels for the continuous return periods
return_levels_1 <- sapply(return_periods, function(period) {
  p <- 1 - 1 / period
  evd::qgev(p, mu[1], sigma[1], xi)
})

# calculate confidence intervals for return levels
alpha <- 0.05  # 95% confidence interval
z_value <- qnorm(1 - alpha / 2)  # Z-value for the confidence interval

# calculate standard errors (SE) for scale and add confidence intervals
se_scale <- z_value * sigma[1] * 0.1  # avoid non-positive values

lower_bound_1 <- sapply(return_periods_continuous, function(period) {
  p <- 1 - 1 / period
  qgev(p, mu[1], pmax(0.0001, sigma[1] - se_scale), xi)  # avoid negative values
})

upper_bound_1 <- sapply(return_periods_continuous, function(period) {
  p <- 1 - 1 / period
  qgev(p, mu[1], sigma[1] + se_scale, xi)
})

# calculate the empirical return levels using block maxima
empirical_return_levels_1 <- quantile(block_maxima_df[,1], probs = 1 - 1/return_periods_continuous)

# set y-axis limits to ensure all points are visible
y_limits <- range(c(return_levels_1, empirical_return_levels_1, lower_bound_1, upper_bound_1))

# Return level plot
plot(return_periods, return_levels_1, type = "l", col = "black", lwd = 2,
     xlab = "Return Period (Years)", ylab = "Return Level",
     main = "Location 1", ylim = y_limits)

# add confidence interval dashed lines
lines(return_periods_continuous, lower_bound_1, col = "red", lty = 2)
lines(return_periods_continuous, upper_bound_1, col = "red", lty = 2)
points(return_periods_continuous, empirical_return_levels_1, col = "blue", pch = 19, cex = 0.5)

legend("topright", 
       legend = c("Return Level", "95% CI", "Empirical"), 
       col = c("black", "red", "blue"), 
       lty = c(1, 2, NA), 
       lwd = c(2, 1, NA), 
       pch = c(NA, NA, 19), 
       cex = 0.55)




## Location 22

# Calculate return levels for the specified return periods
return_levels_22 <- sapply(return_periods, function(period) {
  p <- 1 - 1 / period
  evd::qgev(p, mu_loc22, sigma_loc22, xi_loc22)
})

# Calculate confidence intervals for return levels
alpha <- 0.05  # 95% confidence interval
z_value <- qnorm(1 - alpha / 2)  # Z-value for the confidence interval

# Calculate standard errors (SE) for scale and add confidence intervals
se_scale_22 <- z_value * sigma_loc22 * 0.1  # Adjust as necessary to avoid non-positive values

lower_bound_22 <- sapply(return_periods_continuous, function(period) {
  p <- 1 - 1 / period
  qgev(p, mu_loc22, pmax(0.0001, sigma_loc22 - se_scale_22), xi_loc22)
})

upper_bound_22 <- sapply(return_periods_continuous, function(period) {
  p <- 1 - 1 / period
  qgev(p, mu_loc22, sigma_loc22 + se_scale_22, xi_loc22)
})

# Calculate the empirical return levels using block maxima
empirical_return_levels_22 <- quantile(block_maxima_df[, 22], probs = 1 - 1/return_periods_continuous)

# Set y-axis limits to ensure all points are visible
y_limits_22 <- range(c(return_levels_22, empirical_return_levels_22, lower_bound_22, upper_bound_22))

# Return level plot for Location 22
plot(return_periods, return_levels_22, type = "l", col = "black", lwd = 2,
     xlab = "Return Period (Years)", ylab = "Return Level",
     main = "Location 22", ylim = y_limits_22)

# Add confidence interval dashed lines
lines(return_periods_continuous, lower_bound_22, col = "red", lty = 2)
lines(return_periods_continuous, upper_bound_22, col = "red", lty = 2)
points(return_periods_continuous, empirical_return_levels_22, col = "blue", pch = 19, cex = 0.5)

# Add a legend to the plot
legend("topright", 
       legend = c("Return Level", "95% CI", "Empirical"), 
       col = c("black", "red", "blue"), 
       lty = c(1, 2, NA), 
       lwd = c(2, 1, NA), 
       pch = c(NA, NA, 19), 
       cex = 0.55)



# location 44
# Calculate return levels for the specified return periods
return_levels_44 <- sapply(return_periods, function(period) {
  p <- 1 - 1 / period
  evd::qgev(p, mu_loc44, sigma_loc44, xi_loc44)
})

# Calculate confidence intervals for return levels
alpha <- 0.05  # 95% confidence interval
z_value <- qnorm(1 - alpha / 2)  # Z-value for the confidence interval

# Calculate standard errors (SE) for scale and add confidence intervals
se_scale_44 <- z_value * sigma_loc44 * 0.1  # Adjust as necessary to avoid non-positive values

lower_bound_44 <- sapply(return_periods_continuous, function(period) {
  p <- 1 - 1 / period
  qgev(p, mu_loc44, pmax(0.0001, sigma_loc44 - se_scale_44), xi_loc44)
})

upper_bound_44 <- sapply(return_periods_continuous, function(period) {
  p <- 1 - 1 / period
  qgev(p, mu_loc44, sigma_loc44 + se_scale_44, xi_loc44)
})

# Calculate the empirical return levels using block maxima
empirical_return_levels_44 <- quantile(block_maxima_df[, 44], probs = 1 - 1/return_periods_continuous)

# Set y-axis limits to ensure all points are visible
y_limits_44 <- range(c(return_levels_44, empirical_return_levels_44, lower_bound_44, upper_bound_44))

# Return level plot for Location 44
plot(return_periods, return_levels_44, type = "l", col = "black", lwd = 2,
     xlab = "Return Period (Years)", ylab = "Return Level",
     main = "Location 44", ylim = y_limits_44)

# Add confidence interval dashed lines
lines(return_periods_continuous, lower_bound_44, col = "red", lty = 2)
lines(return_periods_continuous, upper_bound_44, col = "red", lty = 2)
points(return_periods_continuous, empirical_return_levels_44, col = "blue", pch = 19, cex = 0.5)

# Add a legend to the plot
legend("topright", 
       legend = c("Return Level", "95% CI", "Empirical"), 
       col = c("black", "red", "blue"), 
       lty = c(1, 2, NA), 
       lwd = c(2, 1, NA), 
       pch = c(NA, NA, 19), 
       cex = 0.55)



# location 70 

# Calculate return levels for the specified return periods
return_levels_70 <- sapply(return_periods, function(period) {
  p <- 1 - 1 / period
  evd::qgev(p, mu_loc70, sigma_loc70, xi_loc70)
})

# Calculate confidence intervals for return levels
alpha <- 0.05  # 95% confidence interval
z_value <- qnorm(1 - alpha / 2)  # Z-value for the confidence interval

# Calculate standard errors (SE) for scale and add confidence intervals
se_scale_70 <- z_value * sigma_loc70 * 0.1  # Adjust as necessary to avoid non-positive values

lower_bound_70 <- sapply(return_periods_continuous, function(period) {
  p <- 1 - 1 / period
  qgev(p, mu_loc70, pmax(0.0001, sigma_loc70 - se_scale_70), xi_loc70)
})

upper_bound_70 <- sapply(return_periods_continuous, function(period) {
  p <- 1 - 1 / period
  qgev(p, mu_loc70, sigma_loc70 + se_scale_70, xi_loc70)
})

# Calculate the empirical return levels using block maxima
empirical_return_levels_70 <- quantile(block_maxima_df[, 70], probs = 1 - 1/return_periods_continuous)

# Set y-axis limits to ensure all points are visible
y_limits_70 <- range(c(return_levels_70, empirical_return_levels_70, lower_bound_70, upper_bound_70))

# Return level plot for Location 70
plot(return_periods, return_levels_70, type = "l", col = "black", lwd = 2,
     xlab = "Return Period (Years)", ylab = "Return Level",
     main = "Location 70", ylim = y_limits_70)

# Add confidence interval dashed lines
lines(return_periods_continuous, lower_bound_70, col = "red", lty = 2)
lines(return_periods_continuous, upper_bound_70, col = "red", lty = 2)
points(return_periods_continuous, empirical_return_levels_70, col = "blue", pch = 19, cex = 0.5)

# Add a legend to the plot
legend("topright", 
       legend = c("Return Level", "95% CI", "Empirical"), 
       col = c("black", "red", "blue"), 
       lty = c(1, 2, NA), 
       lwd = c(2, 1, NA), 
       pch = c(NA, NA, 19), 
       cex = 0.55)



# location 90
return_levels_90 <- sapply(return_periods, function(period) {
  p <- 1 - 1 / period
  evd::qgev(p, mu_loc90, sigma_loc90, xi_loc90)
})

# Calculate confidence intervals for return levels
alpha <- 0.05  # 95% confidence interval
z_value <- qnorm(1 - alpha / 2)  # Z-value for the confidence interval

# Calculate standard errors (SE) for scale and add confidence intervals
se_scale_90 <- z_value * sigma_loc90 * 0.1  # Adjust as necessary to avoid non-positive values

lower_bound_90 <- sapply(return_periods_continuous, function(period) {
  p <- 1 - 1 / period
  qgev(p, mu_loc90, pmax(0.0001, sigma_loc90 - se_scale_90), xi_loc90)
})

upper_bound_90 <- sapply(return_periods_continuous, function(period) {
  p <- 1 - 1 / period
  qgev(p, mu_loc90, sigma_loc90 + se_scale_90, xi_loc90)
})

# Calculate the empirical return levels using block maxima
empirical_return_levels_90 <- quantile(block_maxima_df[, 90], probs = 1 - 1/return_periods_continuous)

# Set y-axis limits to ensure all points are visible
y_limits_90 <- range(c(return_levels_90, empirical_return_levels_90, lower_bound_90, upper_bound_90))

# Return level plot for Location 90
plot(return_periods, return_levels_90, type = "l", col = "black", lwd = 2,
     xlab = "Return Period (Years)", ylab = "Return Level",
     main = "Location 90", ylim = y_limits_90)

# Add confidence interval dashed lines
lines(return_periods_continuous, lower_bound_90, col = "red", lty = 2)
lines(return_periods_continuous, upper_bound_90, col = "red", lty = 2)
points(return_periods_continuous, empirical_return_levels_90, col = "blue", pch = 19, cex = 0.5)

# Add a legend to the plot
legend("topright", 
       legend = c("Return Level", "95% CI", "Empirical"), 
       col = c("black", "red", "blue"), 
       lty = c(1, 2, NA), 
       lwd = c(2, 1, NA), 
       pch = c(NA, NA, 19), 
       cex = 0.55)




# location 108

return_levels_108 <- sapply(return_periods, function(period) {
  p <- 1 - 1 / period
  evd::qgev(p, mu_loc108, sigma_loc108, xi_loc108)
})

# Calculate confidence intervals for return levels
alpha <- 0.05  # 95% confidence interval
z_value <- qnorm(1 - alpha / 2)  # Z-value for the confidence interval

# Calculate standard errors (SE) for scale and add confidence intervals
se_scale_108 <- z_value * sigma_loc108 * 0.1  # Adjust as necessary to avoid non-positive values

lower_bound_108 <- sapply(return_periods_continuous, function(period) {
  p <- 1 - 1 / period
  qgev(p, mu_loc108, pmax(0.0001, sigma_loc108 - se_scale_108), xi_loc108)
})

upper_bound_108 <- sapply(return_periods_continuous, function(period) {
  p <- 1 - 1 / period
  qgev(p, mu_loc108, sigma_loc108 + se_scale_108, xi_loc108)
})

# Calculate the empirical return levels using block maxima
empirical_return_levels_108 <- quantile(block_maxima_df[, 108], probs = 1 - 1/return_periods_continuous)

# Set y-axis limits to ensure all points are visible
y_limits_108 <- range(c(return_levels_108, empirical_return_levels_108, lower_bound_108, upper_bound_108))

# Return level plot for Location 108
plot(return_periods, return_levels_108, type = "l", col = "black", lwd = 2,
     xlab = "Return Period (Years)", ylab = "Return Level",
     main = "Location 108", ylim = y_limits_108)

# Add confidence interval dashed lines
lines(return_periods_continuous, lower_bound_108, col = "red", lty = 2)
lines(return_periods_continuous, upper_bound_108, col = "red", lty = 2)
points(return_periods_continuous, empirical_return_levels_108, col = "blue", pch = 19, cex = 0.5)

# Add a legend to the plot
legend("topright", 
       legend = c("Return Level", "95% CI", "Empirical"), 
       col = c("black", "red", "blue"), 
       lty = c(1, 2, NA), 
       lwd = c(2, 1, NA), 
       pch = c(NA, NA, 19), 
       cex = 0.55)






############################################################
# Multivariate - Further Research
#
# The final part provides a simple exploratory comparison of
# SST values between selected locations, focusing on summer
# observations and pairwise relationships.
############################################################


# Assuming summer_SST is a matrix where columns represent different locations
summer_SST <- RedSST[month %in% c(6, 7, 8), ]
summer_SST <- as.matrix(summer_SST)

# Extract SST data for location 1
sst_location1 <- summer_SST[, 1]

# List of locations to compare with location 1
locations_to_compare <- c(2, 22, 44, 70, 90, 108)

# Set up plotting area for multiple plots
par(mfrow = c(2, 3))  # Adjust the layout as needed

# Loop through the list of locations and create scatter plots
for (i in locations_to_compare) {
  # Extract SST data for the current location
  sst_current_location <- summer_SST[, i]
  
  # Create the scatter plot
  plot(sst_location1, sst_current_location,
       xlab = "SST at Location 1",
       ylab = paste("SST at Location", i),
       main = paste("SST at Location 1 vs. Location", i),
       col = "blue",
       pch = 19,  # Point type
       cex = 0.7)  # Point size
}




############################################################
# Additional pairwise comparisons
#
# These plots compare locations with very different water
# depths in order to inspect possible differences in their
# SST behaviour.
############################################################

# 23 has highest 1755m depth, 25 has only 4m depth but they are very close
par(mfrow = c(1, 3))

sst_location23 <- summer_SST[, 23]
sst_location25 <- summer_SST[, 25]
plot(sst_location23, sst_location25,
     xlab = "SST at Location 23",
     ylab = paste("SST at Location 25"),
     main = paste("SST at Location 23 vs. Location 25"),
     col = "blue",
     pch = 19,  # Point type
     cex = 0.7)  # Point size


# 66 has depth of 2m
sst_location23 <- summer_SST[, 23]
sst_location66 <- summer_SST[, 66]
plot(sst_location23, sst_location66,
     xlab = "SST at Location 23",
     ylab = paste("SST at Location 66"),
     main = paste("SST at Location 23 vs. Location 66"),
     col = "blue",
     pch = 19,  # Point type
     cex = 0.7)  # Point size

# 104 has depth of 2m
sst_location23 <- summer_SST[, 23]
sst_location104 <- summer_SST[, 104]
plot(sst_location23, sst_location104,
     xlab = "SST at Location 23",
     ylab = paste("SST at Location 104"),
     main = paste("SST at Location 23 vs. Location 104"),
     col = "blue",
     pch = 19,  # Point type
     cex = 0.7)  # Point size











