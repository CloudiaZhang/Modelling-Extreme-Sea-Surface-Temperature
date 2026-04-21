############################################################
# Simulate Data ---- Report Page 12-17
#
# This section performs a simulation study to examine how the
# GEV model behaves under different sample sizes and block
# sizes. Daily SST values are simulated from a normal
# distribution with mean 20 and standard deviation 5.
#
# For each simulated dataset:
#   - block maxima are extracted
#   - GEV models are fitted
#   - diagnostic plots are inspected
#
# The simulation is then extended to compare estimated return
# levels with the corresponding true normal quantiles.
############################################################




############################################################
# Extract Block Maxima Function
#
# This function divides a simulated SST time series into
# blocks of equal length and extracts the maximum value from
# each block. These maxima are then used to fit a GEV model.
############################################################
# function to extract block maxima and plot histogram
Extract_Block_Maxima <- function(num_years, block_size, simulated_SST) {
  num_blocks <- num_years / block_size
  block_size_days <- 365 * block_size
  
  # initialize a numeric vector to store block maxima
  block_maxima <- numeric(num_blocks)
  for (i in 1:num_blocks) {
    # calculate the start and end indices for each block
    start_index <- block_size_days * (i - 1) + 1
    end_index <- block_size_days * i
    block_maxima[i] <- max(simulated_SST[start_index:end_index])
  }
  return(block_maxima)
}




############################################################
# GEV-Fitting Function
#
# This function fits a GEV model to a given set of block
# maxima using the ismev package.
############################################################
# function to fit GEV model
fit_gev_simulate <- function(block_maxima, block_size) {
  library(ismev)
  gev_fit <- gev.fit(xdat = block_maxima, method = "Nelder-Mead", maxit = 10000)
  return(gev_fit)
}





############################################################
# Simulation with 10 years of daily data
#
# A short sample is used first to see how the GEV model
# performs when only a limited amount of data is available.
############################################################
set.seed(5623)
num_years <- 10
num_days <- 365 * num_years

# simulated SST data from normal distribution
simulated_SST_10 <- rnorm(365 * 10, mean = 20, sd = 5)

# block size
block_sizes <- c(1, 2, 5)

# create a list to store the fitted results
results_10 <- list()

par(mfrow = c(1, 3))
# loop through block sizes
for (block_size in block_sizes) {
  # extract block maxima
  maxima <- Extract_Block_Maxima(10, block_size, simulated_SST_10)
  # fit GEV model to the block maxima
  results_10[[paste(block_size, "years")]] <- fit_gev_simulate(maxima)
}

# diagnostic plots
par(mfrow = c(1, 3)) 
for (block_size in block_sizes) {
  gev_fit_10 <- results_10[[paste(block_size, "years")]]
  gev.diag(gev_fit_10)
}





############################################################
# Simulation with 20 years of daily data
#
# The same procedure is repeated with a larger sample size to
# compare the effect of having more data.
############################################################
set.seed(5623)
num_years <- 20
num_days <- 365 * num_years

# simulated SST data from normal distribution
simulated_SST_20 <- rnorm(365 * 20, mean = 20, sd = 5)

# block size
block_sizes <- c(1, 2, 5)
# create a list to store the fitted results
results_20 <- list()

par(mfrow = c(1, 3))
# loop through block sizes
for (block_size in block_sizes) {
  # extract block maxima
  maxima <- Extract_Block_Maxima(20, block_size, simulated_SST_20)
  # fit GEV model to the block maxima
  results_20[[paste(block_size, "years")]] <- fit_gev_simulate(maxima)
}

# diagnostic plots
par(mfrow = c(1, 3)) 
for (block_size in block_sizes) {
  gev_fit_20 <- results_20[[paste(block_size, "years")]]
  gev.diag(gev_fit_20)
}






############################################################
# Simulation with 50 years of daily data
#
# This case allows the model fit to be examined under a much
# larger sample of block maxima.
############################################################
set.seed(5623)
num_years <- 50
num_days <- 365 * num_years

# simulated SST data from normal distribution
simulated_SST_50 <- rnorm(365 * 50, mean = 20, sd = 5)

# block size
block_sizes <- c(1, 2, 5)
# create a list to store the fitted results
results_50 <- list()

par(mfrow = c(1, 3))
# loop through block sizes
for (block_size in block_sizes) {
  # extract block maxima
  maxima <- Extract_Block_Maxima(50, block_size, simulated_SST_50)
  # fit GEV model to the block maxima
  results_50[[paste(block_size, "years")]] <- fit_gev_simulate(maxima)
}

# diagnostic plots
par(mfrow = c(1, 3)) 
for (block_size in block_sizes) {
  gev_fit_50 <- results_50[[paste(block_size, "years")]]
  gev.diag(gev_fit_50)
}






############################################################
# Simulation with 100 years of daily data
#
# This is the largest simulated dataset and is used later for
# the return-level study.
############################################################
set.seed(5623)
num_years <- 100
num_days <- 365 * num_years

# simulated SST data from normal distribution
simulated_SST_100 <- rnorm(365 * 100, mean = 20, sd = 5)

# block size
block_sizes <- c(1, 2, 5, 10)
# create a list to store the fitted results
results_100 <- list()

par(mfrow = c(1, 3))
# loop through block sizes
for (block_size in block_sizes) {
  # extract block maxima
  maxima <- Extract_Block_Maxima(100, block_size, simulated_SST_100)
  # fit GEV model to the block maxima
  results_100[[paste(block_size, "years")]] <- fit_gev_simulate(maxima)
}

# diagnostic plots
par(mfrow = c(1, 3)) 
for (block_size in block_sizes) {
  gev_fit_100 <- results_100[[paste(block_size, "years")]]
  gev.diag(gev_fit_100)
}






############################################################
# Return Level
#
# The next part compares estimated GEV return levels with the
# corresponding true return levels from the underlying normal
# distribution used in the simulation.
############################################################


############################################################
# Return levels for 1-year block size
#
# A first example is considered using annual block maxima.
############################################################
library(evd)
library(ismev)

block_maxima_100_1yr <- Extract_Block_Maxima(100, 1, simulated_SST_100)

gev_fit_100_1yr <- fit_gev_simulate(block_maxima_100_1yr)

return_periods <- c(2, 5, 10, 20, 50, 100)

# return levels
gev_return_levels <- evd::qgev(1 - 1/return_periods, loc = gev_fit_100_1yr$mle[1], scale = gev_fit_100_1yr$mle[2], shape = gev_fit_100_1yr$mle[3])

names(gev_return_levels) <- return_periods

print(gev_return_levels)

true_return_levels <- sapply(return_periods, function(T) {
  qnorm((1 - 1/T)^(1/365), mean = 20, sd = 5)
})
names(true_return_levels) <- return_periods

print(true_return_levels)

plot(return_periods, gev_return_levels, type = "b", pch = 19, xlab = "Return Period (years)", ylab = "Return Level (¡ãC)", main = "Return Level Plot", col = "blue")
lines(return_periods, gev_return_levels, col = "blue", lwd = 2)
points(return_periods, true_return_levels, col = "red", pch = 19)
lines(return_periods, true_return_levels, col = "red", lwd = 2)
legend("bottomright", legend = c("GEV Return Levels", "True Normal Quantiles"), col = c("blue", "red"), pch = 19)

comparison <- data.frame(
  Return_Period = return_periods,
  GEV_Return_Level = gev_return_levels,
  True_Return_Level = true_return_levels
)

print(comparison)





############################################################
# Return levels for 10-year block size
#
# The same comparison is repeated using a much wider block
# size, which changes both the number of block maxima and the
# fitted return levels.
############################################################
library(evd)

block_maxima_100_10yr <- Extract_Block_Maxima(100, 10, simulated_SST_100)

gev_fit_100_10yr <- fit_gev_simulate(block_maxima_100_10yr)

return_periods <- c(50, 100, 200, 500)

# return levels
gev_return_levels <- evd::qgev(1 - 1/return_periods, loc = gev_fit_100_10yr$mle[1], scale = gev_fit_100_10yr$mle[2], shape = gev_fit_100_10yr$mle[3])

names(gev_return_levels) <- return_periods

print(gev_return_levels)

true_return_levels <- sapply(return_periods, function(T) {
  qnorm((1 - 1/T)^(1/365), mean = 20, sd = 5)
})
names(true_return_levels) <- return_periods

print(true_return_levels) 

# plot return levels
plot(return_periods, gev_return_levels, type = "b", pch = 19, xlab = "Return Period (years)", ylab = "Return Level (¡ãC)", main = "Return Level Plot", col = "blue", ylim = range(c(gev_return_levels, true_return_levels)))
lines(return_periods, gev_return_levels, col = "blue", lwd = 2)

# add true return levels
points(return_periods, true_return_levels, col = "red", pch = 19)
lines(return_periods, true_return_levels, col = "red", lwd = 2)

# add legend
legend("bottomright", legend = c("GEV Return Levels", "True Normal Quantiles"), col = c("blue", "red"), pch = 19, lwd = 2)

comparison <- data.frame(
  Return_Period = return_periods,
  GEV_Return_Level = gev_return_levels,
  True_Return_Level = true_return_levels
)

print(comparison)





############################################################
# Function for repeated return-level simulation
#
# To study the variability of estimated return levels, the
# simulation is repeated 100 times for each block size.
############################################################
library(evd)

# create a function that generate return levels for each return period for each simulated dataset
simulate_gev_rl <- function(num_years, block_size, return_periods, data) {
  
  # create a matrix to store the return levels for each return period
  sim_rl_list <- matrix(NA, nrow = 100, ncol = length(return_periods))
  colnames(sim_rl_list) <- return_periods
  
  for (i in 1:100) {
    # choose the corresponding dataset
    simulated_SST <- sim_total[[i]]
    # extract block maxima
    block_maxima <- Extract_Block_Maxima(num_years, block_size, simulated_SST)
    # fit GEV model
    gev_fit <- gev.fit(block_maxima, method = "Nelder-Mead", maxit = 10000)
    
    # calculate return levels
    sim_rl <- round(evd::qgev(1 - 1/((return_periods/block_size)), loc = gev_fit$mle[1], scale = gev_fit$mle[2], shape = gev_fit$mle[3]), 2)
    sim_rl_list[i, ] <- sim_rl
  }
  
  return(as.data.frame(sim_rl_list))
}


set.seed(1111)
# create a list of 100 sets of simulated data
sim_total <- lapply(1:100, function(x) rnorm(num_years * 365, mean = 20, sd = 5))

# parameters
num_years <- 100
block_sizes <- c(1, 2, 5, 10)
return_periods <- c(50, 100, 200, 500)
mean <- 20
sd <- 5
data <- sim_total

# true return levels for the normal distribution
true_return_levels <- sapply(return_periods, function(T) {
  qnorm((1 - 1/T)^(1/365), mean = mean, sd = sd)
})
names(true_return_levels) <- return_periods

# return levels for each block size
rl_1 <- simulate_gev_rl(num_years, 1, return_periods, sim_total)
rl_2 <- simulate_gev_rl(num_years, 2, return_periods, sim_total)
rl_5 <- simulate_gev_rl(num_years, 5, return_periods, sim_total)
rl_10 <- simulate_gev_rl(num_years, 10, return_periods, sim_total)





############################################################
# Plots of simulated return levels
#
# Histograms are used to compare the distribution of
# estimated return levels with the corresponding true normal
# quantiles, for each block size.
############################################################

# for 1 year block size
par(mfrow = c(2, 2))
# estimated mean return level
mean_rl_1 <- colMeans(rl_1)

# histogram of return levels
for (i in 1:length(return_periods)) {
  hist(rl_1[, i], breaks = 30, main = paste("Return Level for", return_periods[i], "years\nBlock Size: 1 year"), xlab = "Return Level", col = "lightblue", xlim = c(min(rl_1), max(rl_1)))
  abline(v = true_return_levels[i], col = "red", lwd = 2)
  abline(v = mean_rl_1[i], col = "blue", lwd = 2, lty = 2)
  legend("topright", legend = c("True Normal Quantile", "Mean Estimated Return Level"), col = c("red","blue"), lwd = 2, lty = c(1, 2), cex = 1, bty = "n", inset = c(-0.35, -0.1))
}

par(mfrow = c(2, 2))
mean_rl_2 <- colMeans(rl_2)

# for 2-year block size
for (i in 1:length(return_periods)) {
  hist(rl_2[, i], breaks = 30, main = paste("Return Level for", return_periods[i], "years\nBlock Size: 2 years"), xlab = "Return Level", col = "lightblue", xlim = c(min(rl_2), max(rl_2)))
  abline(v = true_return_levels[i], col = "red", lwd = 2)
  abline(v = mean_rl_2[i], col = "blue", lwd = 2, lty = 2)
  legend("topright", legend = c("True Normal Quantile", "Mean Estimated Return Level"), col = c("red","blue"), lwd = 2, lty = c(1, 2), cex = 1, bty = "n", inset = c(-0.35, -0.1))
}

par(mfrow = c(2, 2))
mean_rl_5 <- colMeans(rl_5)

# for 5-year block size
for (i in 1:length(return_periods)) {
  hist(rl_5[, i], breaks = 30, main = paste("Return Level for", return_periods[i], "years\nBlock Size: 5 years"), xlab = "Return Level", col = "lightblue")
  abline(v = true_return_levels[i], col = "red", lwd = 2)
  abline(v = mean_rl_5[i], col = "blue", lwd = 2, lty = 2)
  legend("topright", legend = c("True Normal Quantile", "Mean Estimated Return Level"), col = c("red","blue"), lwd = 2, lty = c(1, 2), cex = 1, bty = "n", inset = c(-0.35, -0.1))
}

par(mfrow = c(2, 2))
mean_rl_10 <- colMeans(rl_10)

# for 10-year block size
for (i in 1:length(return_periods)) {
  hist(rl_10[, i], breaks = 30, main = paste("Return Level for", return_periods[i], "years\nBlock Size: 10 years"), xlab = "Return Level", col = "lightblue")
  abline(v = true_return_levels[i], col = "red", lwd = 2)
  abline(v = mean_rl_10[i], col = "blue", lwd = 2, lty = 2)
  legend("topright", legend = c("True Normal Quantile", "Mean Estimated Return Level"), col = c("red","blue"), lwd = 2, lty = c(1, 2), cex = 1, bty = "n", inset = c(-0.4, -0.1))
}
