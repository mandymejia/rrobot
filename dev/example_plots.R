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

set.seed(2025)
plot_RD(kurt_data$hk, log = TRUE, show_f_density = TRUE)


# Random
norm_dat = mvrnorm(500, mu = rep(0, 5), Sigma = diag(5))
plot_RD(norm_dat)


###############################################################################
# RE-FACTORED FOR S3                                                          #
###############################################################################
#########################################
# Method: F                             #
#########################################
data_matrix <- fMRIscrub::Dat1
kurt_data <- ICA_extract_kurt(time_series = data_matrix)

result_F <- RD(x = kurt_data$hk,
               method = "F",
               mode = "auto",
               quantile = 0.01)

set.seed(2025)
RD_obj <- compute_RD(x = kurt_data$hk, mode = "auto")
result_F_thresh <- threshold_RD(x = kurt_data$hk,
                                method = "F",
                                RD_obj = RD_obj,
                                quantile = 0.01)


plot(result_F_thresh)

