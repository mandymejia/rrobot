library(MASS)
library(expm)
library(cellWise)
library(forecast)
library(fMRIscrub)
library(isotree)
library(future.apply)
library(future)
library(tibble)
library(here)
library(fastICA)

devtools::load_all()
source(here::here("dev", "ICA_extract_kurt.R"))

set.seed(2025)

data_matrix <- fMRIscrub::Dat1
kurt_data <- ICA_extract_kurt(time_series = data_matrix)
hk_data <- kurt_data$hk


out_result <- univOut(x = hk_data, cutoff = 4, method = "SHASH")
RD_obj_shash <- compute_RD(x = out_result$x_norm, mode = "auto", dist = TRUE)

