# --------------------------------------------
# WFCS Probabalistic Risk Assessment
#
# Code Author: David Dayan
# Code QAQC: Brenda Hanley
# Project Supervisor: Krysten Schuler
#
# Project: Waterfowl Contaminants in the NE US. 
# 
# Location: Cornell Wildlife Health Lab
# Date: February 2024
# License: MIT
# 
# Description: This code runs the probabilistic risk assessment.
# 
# Funding: This work was supported under the project 
# name 'Contaminant Loads in Waterfowl of the Northeast Atlantic
# Flyway: New Threats and Outdated Advisories' of the 
# Multistate Conservation Grant Program.
#
# --------------------------------------------

# Load required packages.
library(tidyverse) # Version 2.0.0
library(knitr)	# Version 1.43
library(fitdistrplus) # Version 1.1.11

# To ensure reproducibility.
set.seed(720)

# Load the data.
risk_data <- read_csv("WFCS_group_risk_levels.csv")

#--------------------------------------------------------------------------
# Fitting HQ_TEQ ----------------------------------------------------------

# Load the data.
HQ_TEQ=c(na.omit(risk_data$HQ_TEQ))

# Divide data to go from 2 meals/month to 1 meal/year.
HQ_TEQ_MEAL=HQ_TEQ/24 

# Look at the shape of the data.
summary(HQ_TEQ_MEAL)
hist(HQ_TEQ_MEAL)

# Note: There are only 4 zeros. Extreme J shape. 
# See if we can transform the J shape to normal. 
# Go up or down the roots to see if we can get the best fit RE: diagnostics.
HQ_TEQ_MEAL_S=HQ_TEQ_MEAL^(1/5)
hist(HQ_TEQ_MEAL_S)

# Fit the transformed data to a normal distribution.
fit_TEQ_HQ=fitdist(HQ_TEQ_MEAL_S, "norm")

# Look at the fit metrics.
gofstat(fit_TEQ_HQ)
summary(fit_TEQ_HQ)

# Look at the diagnostics.
plot(fit_TEQ_HQ)

# Save the plot as a PNG with 1000 dpi and 6in x 6in dimensions
png("TEQ_HQ_diagnostics.png", width = 6.5, height = 6, units = "in", res = 1000)
plot(fit_TEQ_HQ)  # Plot the diagnostics
dev.off()  # Close the PNG device

#--------------------------------------------------------------------------
# Fitting CR_TEQ ----------------------------------------------------------

# Load the data.
CR_TEQ=c(na.omit(risk_data$CR_TEQ))

# Divide data to go from 2 meals/month to 1 meal/year.
CR_TEQ_MEAL=CR_TEQ/24 

# Look at the shape of the data.
summary(CR_TEQ_MEAL)
hist(CR_TEQ_MEAL)

# See if we can transform the J shape to normal. 
# Go up or down the roots to see if we can get the best fit RE diagnostics.
CR_TEQ_MEAL_S=CR_TEQ_MEAL^(1/5)
hist(CR_TEQ_MEAL_S)

# Fit the transformed data to a normal distribution.
fit_TEQ_CR=fitdist(CR_TEQ_MEAL_S, "norm")

# Look at the fit metrics.
gofstat(fit_TEQ_CR)
summary(fit_TEQ_CR)

# Look at the diagnostics.
plot(fit_TEQ_CR)

# Save the plot as a PNG with 1000 dpi and 6in x 6in dimensions
png("TEQ_CR_diagnostics.png", width = 6.5, height = 6, units = "in", res = 1000)
plot(fit_TEQ_CR)  # Plot the diagnostics
dev.off()  # Close the PNG device

#--------------------------------------------------------------------------
# Simulating TEQ Risk -----------------------------------------------------

# Number of iterations.
num_iterations <- 10000

# Initialize a dataframe to store results.
TEQ_risks <- data.frame(
  num_samples = numeric(),
  mean_HQ = numeric(),percentile_95_HQ = numeric(),percentile_99_HQ = numeric(),
  mean_CR = numeric(),percentile_95_CR = numeric(),percentile_99_CR = numeric(),
  stringsAsFactors = FALSE)

# Loop over different values of num_samples_per_iteration.
for (num_samples_per_iteration in 1:60) {
  # Initialize vectors to store sums of HQ and CR values
  sums_HQ <- numeric(num_iterations)
  sums_CR <- numeric(num_iterations)
  # Perform simulations.
  for (i in 1:num_iterations) {
    # Simulate HQ values.
    samples_TEQ_HQ <- rnorm(num_samples_per_iteration, mean = fit_TEQ_HQ$estimate[1], sd = fit_TEQ_HQ$estimate[2])
    samples_TEQ_HQ <- samples_TEQ_HQ^5
    sum_TEQ_HQ <- sum(samples_TEQ_HQ)
    # Simulate CR values.
    samples_TEQ_CR <- rnorm(num_samples_per_iteration, mean = fit_TEQ_CR$estimate[1], sd = fit_TEQ_CR$estimate[2])
    samples_TEQ_CR <- samples_TEQ_CR^5
    sum_TEQ_CR <- sum(samples_TEQ_CR)
    # Store sums.
    sums_HQ[i] <- sum_TEQ_HQ
    sums_CR[i] <- sum_TEQ_CR
  } # End num samples per iteration.
  
  # Calculate statistics for HQ sums.
  mean_sum_HQ <- mean(sums_HQ)
  percentile_95_HQ <- quantile(sums_HQ, 0.95)
  percentile_99_HQ <- quantile(sums_HQ, 0.99)
  
  # Calculate statistics for CR sums.
  mean_sum_CR <- mean(sums_CR)
  percentile_95_CR <- quantile(sums_CR, 0.95)
  percentile_99_CR <- quantile(sums_CR, 0.99)
  
  # Add results to data frame.
  TEQ_risks <- rbind(TEQ_risks, list(num_samples_per_iteration,
                                     mean_sum_HQ, percentile_95_HQ, percentile_99_HQ, 
                                     mean_sum_CR, percentile_95_CR, percentile_99_CR))
  
} # End loop over different values of num_samples_per_iteration.

# Name columns.
colnames(TEQ_risks) <- c("num_samples", "mean_HQ", "percentile_95_HQ", "percentile_99_HQ", "mean_CR", "percentile_95_CR", "percentile_99_CR")

# Display results.
view(TEQ_risks)

#--------------------------------------------------------------------------
# Fitting HQ_TPCB ---------------------------------------------------------

# Load the data.
HQ_TPCB=c(na.omit(risk_data$HQ_TPCB))

# Divide the data by the number of meals to get the by-meal vector.
HQ_TPCB_MEAL=HQ_TPCB/24

# Look at the shape of the data.
summary(HQ_TPCB_MEAL)
hist(HQ_TPCB_MEAL)

# See if we can transform the J shape to normal. 
# Go up or down the roots to see if we can get the best fit RE diagnostics.
HQ_TPCB_MEAL_S=HQ_TPCB_MEAL^(1/200)
hist(HQ_TPCB_MEAL_S)

# Fit the transformed data to a normal distribution.
fit_TPCB_HQ=fitdist(HQ_TPCB_MEAL_S, "norm")

# Look at the metrics of the fit.
gofstat(fit_TPCB_HQ)
summary(fit_TPCB_HQ)

# Look at the diagnostics.
plot(fit_TPCB_HQ)

# Save the plot as a PNG with 1000 dpi and 6in x 6in dimensions
png("PCB_HQ_diagnostics.png", width = 6.5, height = 6, units = "in", res = 1000)
plot(fit_TPCB_HQ)  # Plot the diagnostics
dev.off()  # Close the PNG device

#--------------------------------------------------------------------------
# Fitting CR_TPCB ---------------------------------------------------------

# Load the data.
CR_TPCB=c(na.omit(risk_data$CR_TPCB))

# Divide the data by the number of meals to get the by-meal vector.
CR_TPCB_MEAL=CR_TPCB/24

# Look at the shape of the data.
summary(CR_TPCB_MEAL)
hist(CR_TPCB_MEAL)

# See if we can transform the J shape to normal. 
# Go up or down the roots to see if we can get the best fit RE diagnostics.
CR_TPCB_MEAL_S=CR_TPCB_MEAL^(1/100)
hist(CR_TPCB_MEAL_S)

# Fit the transformed data.
fit_TPCB_CR=fitdist(CR_TPCB_MEAL_S, "norm")

# Look at the metrics of the fit. 
gofstat(fit_TPCB_CR)
summary(fit_TPCB_CR)

# Look at the diagnostics.
plot(fit_TPCB_CR)

# Save the plot as a PNG with 1000 dpi and 6in x 6in dimensions
png("PCB_CR_diagnostics.png", width = 6.5, height = 6, units = "in", res = 1000)
plot(fit_TPCB_CR)  # Plot the diagnostics
dev.off()  # Close the PNG device

#--------------------------------------------------------------------------
# Simulating TPCB Risk -----------------------------------------------------

# Number of iterations.
num_iterations <- 10000

# Initialize a dataframe to store results.
TPCB_risks <- data.frame(
  num_samples = numeric(),
  mean_HQ = numeric(),
  percentile_95_HQ = numeric(),
  percentile_99_HQ = numeric(),
  mean_CR = numeric(),
  percentile_95_CR = numeric(),
  percentile_99_CR = numeric(),
  stringsAsFactors = FALSE)

# Loop over different values of num_samples_per_iteration.
for (num_samples_per_iteration in 1:60) {
  # Initialize vectors to store sums of HQ and CR values.
  sums_HQ <- numeric(num_iterations)
  sums_CR <- numeric(num_iterations)
  # Perform simulations.
  for (i in 1:num_iterations) {
    # Simulate HQ values.
    samples_TPCB_HQ <- rnorm(num_samples_per_iteration, mean = fit_TPCB_HQ$estimate[1], sd = fit_TPCB_HQ$estimate[2])
    samples_TPCB_HQ <- samples_TPCB_HQ^200
    sum_TPCB_HQ <- sum(samples_TPCB_HQ)
    # Simulate CR values.
    samples_TPCB_CR <- rnorm(num_samples_per_iteration, mean = fit_TPCB_CR$estimate[1], sd = fit_TPCB_CR$estimate[2])
    samples_TPCB_CR <- samples_TPCB_CR^100
    sum_TPCB_CR <- sum(samples_TPCB_CR)
    # Store sums.
    sums_HQ[i] <- sum_TPCB_HQ
    sums_CR[i] <- sum_TPCB_CR
  } # End iterations.
  
  # Calculate statistics for HQ sums.
  mean_sum_HQ <- mean(sums_HQ)
  percentile_95_HQ <- quantile(sums_HQ, 0.95)
  percentile_99_HQ <- quantile(sums_HQ, 0.99)
  
  # Calculate statistics for CR sums.
  mean_sum_CR <- mean(sums_CR)
  percentile_95_CR <- quantile(sums_CR, 0.95)
  percentile_99_CR <- quantile(sums_CR, 0.99)
  
  # Add results to dataframe.
  TPCB_risks <- rbind(TPCB_risks, list(num_samples_per_iteration,mean_sum_HQ, percentile_95_HQ, percentile_99_HQ, mean_sum_CR, percentile_95_CR, percentile_99_CR))
  
} # End num_samples_per_iteration loop. 

# Name columns.
colnames(TPCB_risks) <- c("num_samples", "mean_HQ", "percentile_95_HQ", "percentile_99_HQ", "mean_CR", "percentile_95_CR", "percentile_99_CR")

# Display results.
view(TPCB_risks)

#--------------------------------------------------------------------------
# Fitting Sum_OCPs_HQ -----------------------------------------------------

# Load the data
Sum_OCPs_HQ=c(na.omit(risk_data$Sum_OCPs_HQ))

# Divide the data by the number of meals to get the by-meal vector.
Sum_OCPs_HQ_MEAL=Sum_OCPs_HQ/24

# Look at the shape of the data.
summary(Sum_OCPs_HQ_MEAL)
hist(Sum_OCPs_HQ_MEAL)

# See if we can transform the J shape to normal. 
# Go up or down the roots to see if we can get the best fit RE diagnostics.
Sum_OCPs_HQ_MEAL_S=Sum_OCPs_HQ_MEAL^(1/200)
hist(Sum_OCPs_HQ_MEAL_S)

# Fit the transformed data to a normal distribution.
fit_OCP_HQ=fitdist(Sum_OCPs_HQ_MEAL_S, "norm")

# Look at the metrics of the fit. 
gofstat(fit_OCP_HQ)
summary(fit_OCP_HQ)

# Look at the diagnostics.
plot(fit_OCP_HQ)

#--------------------------------------------------------------------------
# Fitting Sum_OCPs_CR -----------------------------------------------------

# Load the data.
Sum_OCPs_CR=c(na.omit(risk_data$Sum_OCPs_CR))

# Divide the data by the number of meals to get the by-meal vector.
Sum_OCPs_CR_MEAL=Sum_OCPs_CR/24

# Look at the shape of the data.
summary(Sum_OCPs_CR_MEAL)
hist(Sum_OCPs_CR_MEAL)

# See if we can transform the J shape to normal. 
# Go up or down the roots to see if we can get the best fit RE diagnostics.
Sum_OCPs_CR_MEAL_S=Sum_OCPs_CR_MEAL^(1/300)
hist(Sum_OCPs_CR_MEAL_S)

# Fit the transformed data to the normal distribution.
fit_OCP_CR=fitdist(Sum_OCPs_CR_MEAL_S, "norm")

# Look at the metrics of the fit.
gofstat(fit_OCP_CR)
summary(fit_OCP_CR)

# Look at the diagnostics.
plot(fit_OCP_CR)

# Save the plot as a PNG with 1000 dpi and 6in x 6in dimensions
png("OCP_CR_diagnostics.png", width = 6.5, height = 6, units = "in", res = 1000)
plot(fit_OCP_CR)  # Plot the diagnostics
dev.off()  # Close the PNG device

#--------------------------------------------------------------------------
# Simulating OCP Risk -----------------------------------------------------

# Number of iterations.
num_iterations <- 10000

# Initialize a dataframe to store results.
OCP_risks <- data.frame(
  num_samples = numeric(),
  mean_HQ = numeric(),
  percentile_95_HQ = numeric(),
  percentile_99_HQ = numeric(),
  mean_CR = numeric(),
  percentile_95_CR = numeric(),
  percentile_99_CR = numeric(),
  stringsAsFactors = FALSE)

# Loop over different values of num_samples_per_iteration.
for (num_samples_per_iteration in 1:60) {
  # Initialize vectors to store sums of HQ and CR values.
  sums_HQ <- numeric(num_iterations)
  sums_CR <- numeric(num_iterations)
  # Perform simulations.
  for (i in 1:num_iterations) {
    # Simulate HQ values.
    samples_OCP_HQ <- rnorm(num_samples_per_iteration, mean = fit_OCP_HQ$estimate[1], sd = fit_OCP_HQ$estimate[2])
    samples_OCP_HQ <- samples_OCP_HQ^200
    sum_OCP_HQ <- sum(samples_OCP_HQ)
    # Simulate CR values.
    samples_OCP_CR <- rnorm(num_samples_per_iteration, mean = fit_OCP_CR$estimate[1], sd = fit_OCP_CR$estimate[2])
    samples_OCP_CR <- samples_OCP_CR^300
    sum_OCP_CR <- sum(samples_OCP_CR)
    # Store sums.
    sums_HQ[i] <- sum_OCP_HQ
    sums_CR[i] <- sum_OCP_CR
  } # End num iterations.
  
  # Calculate statistics for HQ sums.
  mean_sum_HQ <- mean(sums_HQ)
  percentile_95_HQ <- quantile(sums_HQ, 0.95)
  percentile_99_HQ <- quantile(sums_HQ, 0.99)
  
  # Calculate statistics for CR sums.
  mean_sum_CR <- mean(sums_CR)
  percentile_95_CR <- quantile(sums_CR, 0.95)
  percentile_99_CR <- quantile(sums_CR, 0.99)
  
  # Add results to dataframe.
  OCP_risks <- rbind(OCP_risks, list(num_samples_per_iteration,mean_sum_HQ, percentile_95_HQ, percentile_99_HQ, mean_sum_CR, percentile_95_CR, percentile_99_CR))
  
} # End Loop over different values of num_samples_per_iteration.

# Name columns.
colnames(OCP_risks) <- c("num_samples", "mean_HQ", "percentile_95_HQ", "percentile_99_HQ", "mean_CR", "percentile_95_CR", "percentile_99_CR")

# Display results.
view(OCP_risks)

#--------------------------------------------------------------------------
# Fitting Sum_PFAS_HQ -----------------------------------------------------

# Load the data.
Sum_PFAS_HQ=c(na.omit(risk_data$Sum_PFAS_HQ))

# Divide the data by the number of meals to get the by-meal vector.
Sum_PFAS_HQ_MEAL=Sum_PFAS_HQ/24

# Look at the shape of the data.
summary(Sum_PFAS_HQ_MEAL)
hist(Sum_PFAS_HQ_MEAL)

# See if we can transform the J shape to normal. 
# Go up or down the roots to see if we can get the best fit RE diagnostics.
Sum_PFAS_HQ_MEAL_S=Sum_PFAS_HQ_MEAL^(1/5)
hist(Sum_PFAS_HQ_MEAL_S)

# Fit the transformed data to a normal distribution.
fit_PFAS_HQ=fitdist(Sum_PFAS_HQ_MEAL_S, "norm")

# Look at the fit metrics.
gofstat(fit_PFAS_HQ)
summary(fit_PFAS_HQ)

# Look at the diagnostics.
plot(fit_PFAS_HQ)

# Save the plot as a PNG with 1000 dpi and 6in x 6in dimensions
png("PFAS_HI_diagnostics.png", width = 6.5, height = 6, units = "in", res = 1000)
plot(fit_PFAS_HQ)  # Plot the diagnostics
dev.off()  # Close the PNG device

#--------------------------------------------------------------------------
# Simulating PFAS Risk -----------------------------------------------------

# Number of iterations
num_iterations <- 10000

# Initialize a dataframe to store results
PFAS_risks <- data.frame(
  num_samples = numeric(),
  mean_HQ = numeric(),
  percentile_95_HQ = numeric(),
  percentile_99_HQ = numeric(),
  stringsAsFactors = FALSE)

# Loop over different values of num_samples_per_iteration.
for (num_samples_per_iteration in 1:60) {
  # Initialize vectors to store sums of HQ values
  sums_HQ <- numeric(num_iterations)
  # Perform simulations
  for (i in 1:num_iterations) {
    # Simulate HQ values
    samples_PFAS_HQ <- rnorm(num_samples_per_iteration, mean = fit_PFAS_HQ$estimate[1], sd = fit_PFAS_HQ$estimate[2])
    samples_PFAS_HQ <- samples_PFAS_HQ^5
    sum_PFAS_HQ <- sum(samples_PFAS_HQ)
    # Store sums
    sums_HQ[i] <- sum_PFAS_HQ
  } # End iterations.
  
  # Calculate statistics for HQ sums
  mean_sum_HQ <- mean(sums_HQ)
  percentile_95_HQ <- quantile(sums_HQ, 0.95)
  percentile_99_HQ <- quantile(sums_HQ, 0.99)
  
  # Add results to dataframe
  PFAS_risks <- rbind(PFAS_risks, list(num_samples_per_iteration,mean_sum_HQ, percentile_95_HQ, percentile_99_HQ))
  
}  # End loop over different values of num samples per iteration.

# Name columns.
colnames(PFAS_risks) <- c("num_samples", "mean_HQ", "percentile_95_HQ", "percentile_99_HQ")

# Display results.
view(PFAS_risks)

#-------------------------------------------------------------------------
# Total_HI ----------------------------------------------------------------

# Load the data.
Total_HI=c(na.omit(risk_data$Total_HI))

# Divide the data by the number of meals to get the by-meal vector.
Total_HI_MEAL=Total_HI/24

# Look at the shapeof the data.
summary(Total_HI_MEAL)
hist(Total_HI_MEAL)

# See if we can transform the J shape to normal. 
# Go up or down the roots to see if we can get the best fit RE diagnostics.
Total_HI_MEAL_S=Total_HI_MEAL^(1/10)
hist(Total_HI_MEAL_S)

# Fit the transformed data to a normal distribution.
fit_Total_HI=fitdist(Total_HI_MEAL_S, "norm")

# Look at the metrics of the fit. 
gofstat(fit_Total_HI)
summary(fit_Total_HI)

# Look at the diagnostics.
plot(fit_Total_HI)

# Save the plot as a PNG with 1000 dpi and 6in x 6in dimensions
png("Total_HI_diagnostics.png", width = 6.5, height = 6, units = "in", res = 1000)
plot(fit_Total_HI)  # Plot the diagnostics
dev.off()  # Close the PNG device

#--------------------------------------------------------------------------
# Total_CR -------------------------------------------------------------

# Load the data.
Total_CR=c(na.omit(risk_data$Total_CR))

# Divide the data by the number of meals to get the by-meal vector.
Total_CR_MEAL=Total_CR/24

# Look at the shape of the data.
summary(Total_CR_MEAL)
hist(Total_CR_MEAL)

# See if we can transform the J shape to normal. 
# Go up or down the roots to see if we can get the best fit RE diagnostics.
Total_CR_MEAL_S=Total_CR_MEAL^(1/100)
hist(Total_CR_MEAL_S)

# Fit the transformed data to a normal distribution.
fit_total_CR=fitdist(Total_CR_MEAL_S, "norm")

# Look at the metrics of the fit.
gofstat(fit_total_CR)
summary(fit_total_CR)

# Look at the diagnostics.
plot(fit_total_CR)

# Save the plot as a PNG with 1000 dpi and 6in x 6in dimensions
png("Total_CR_diagnostics.png", width = 6.5, height = 6, units = "in", res = 1000)
plot(fit_total_CR)  # Plot the diagnostics
dev.off()  # Close the PNG device

#--------------------------------------------------------------------------
# Simulating Total Risk ---------------------------------------------------

# Number of iterations.
num_iterations <- 10000

# Initialize a dataframe to store results.
Total_risks <- data.frame(
  num_samples = numeric(),
  mean_HQ = numeric(),
  percentile_95_HQ = numeric(),
  percentile_99_HQ = numeric(),
  mean_CR = numeric(),
  percentile_95_CR = numeric(),
  percentile_99_CR = numeric(),
  stringsAsFactors = FALSE)

# Loop over different values of num_samples_per_iteration.
for (num_samples_per_iteration in 1:60) {
  # Initialize vectors to store sums of HQ and CR values.
  sums_HQ <- numeric(num_iterations)
  sums_CR <- numeric(num_iterations)
  # Perform simulations.
  for (i in 1:num_iterations) {
    # Simulate HQ values.
    samples_total_HQ <- rnorm(num_samples_per_iteration, mean = fit_Total_HI$estimate[1], sd = fit_Total_HI$estimate[2])
    samples_total_HQ <- samples_total_HQ^10
    sum_total_HQ <- sum(samples_total_HQ)
    # Simulate CR values.
    samples_total_CR <- rnorm(num_samples_per_iteration, mean = fit_total_CR$estimate[1], sd = fit_total_CR$estimate[2])
    samples_total_CR <- samples_total_CR^100
    sum_total_CR <- sum(samples_total_CR)
    # Store sums.
    sums_HQ[i] <- sum_total_HQ
    sums_CR[i] <- sum_total_CR
  } # End number iterations loop.
  
  # Calculate statistics for HQ sums.
  mean_sum_HQ <- mean(sums_HQ)
  percentile_95_HQ <- quantile(sums_HQ, 0.95)
  percentile_99_HQ <- quantile(sums_HQ, 0.99)
  
  # Calculate statistics for CR sums.
  mean_sum_CR <- mean(sums_CR)
  percentile_95_CR <- quantile(sums_CR, 0.95)
  percentile_99_CR <- quantile(sums_CR, 0.99)
  
  # Add results to dataframe.
  Total_risks <- rbind(Total_risks, list(num_samples_per_iteration,
                                         mean_sum_HQ, percentile_95_HQ, percentile_99_HQ, 
                                         mean_sum_CR, percentile_95_CR, percentile_99_CR))
  
} # End loop over different values of num_samples_per_iteration. 

# Name columns.
colnames(Total_risks) <- c("num_samples", "mean_HQ", "percentile_95_HQ", "percentile_99_HQ", "mean_CR", "percentile_95_CR", "percentile_99_CR")

# Display results.
view(Total_risks)
