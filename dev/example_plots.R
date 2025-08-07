# Load required libraries
library(ggplot2)
library(robustbase)
library(expm)
library(stats)
library(MASS)
library(cellWise)
library(fMRIscrub)
library(isotree)
library(tibble)
library(here)
library(fastICA)

devtools::load_all()
source(here::here("dev", "ICA_extract_kurt.R"))

# ABIDE
data_matrix <- fMRIscrub::Dat1
kurt_data <- ICA_extract_kurt(time_series = data_matrix)
plot_RD(kurt_data$hk, log = TRUE, show_f_density = FALSE)


# Random
norm_dat = mvrnorm(500, mu = rep(0, 5), Sigma = diag(5))
plot_RD(norm_dat)
