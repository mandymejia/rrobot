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
RD_org_obj <- compute_RD(x = hk_data, mode = "auto")






SHASH_F_result <- thresh_SASH(x = hk_data, cutoff = 4, quantile = 0.01)
SHASH_F_result$threshold




HR_result <- thresh_F(Q = ncol(hk_data), n = nrow(hk_data), h = RD_org_obj$h, quantile = 0.01)
HR_result$threshold
summary(HR_result)
