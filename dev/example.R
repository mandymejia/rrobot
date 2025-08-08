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

# Load the data (example using fMRIscrub::Dat1)
data_matrix <- fMRIscrub::Dat1

###------ STEP 0: Extract the low and high kurtosis of the given data------------
kurt_data <- ICA_extract_kurt(time_series = data_matrix)

# Access outputs
dim(kurt_data$hk)        # High-kurtosis ICA components: 193 x 25
dim(kurt_data$lk)        # Low-kurtosis ICA components: 193 x 4
kurt_data$highkurt       # Logical vector indicating high-kurtosis components

###------ STEP 1: to detect univariate outliers.---------------------------------
hk_data <- kurt_data$hk
out_result <- univOut(hk_data = hk_data, cutoff = 4, trans = "SHASH")

###------ STEP 2: to impute univariate outliers.--------------------------------
imp_result <- impTemp_univOut(x = hk_data, outlier_mask = out_result$outliers)

###------ STEP 3: to compute RD and related objects of kurt_data----------------
RD_org_obj <- RD(data_matrix = hk_data, mode = "auto")

str(RD_org_obj)
RD_org <- RD_org_obj$RD
best_org <- RD_org_obj$ind_incld

###------ STEP 4: to obtain the threshold value for detecting outliers----------
#--1---SI-----------------------------------------------------------------------
SI_results <- thresh_SI(RD_org_obj = RD_org_obj
                 ,imp_data = imp_result$imp_data
                 ,alpha = 0.01)

SI_results$SI_threshold
summary(SI_results)

#--2---SIBoot-------------------------------------------------------------------
SI_boot_results <- thresh_SI_boot( RD_org_obj = RD_org_obj
                            ,imp_data = imp_result$imp_data
                            , B = 500, alpha = 0.01, boot_quant = 0.95,
                            verbose = TRUE)
SI_boot_results$LB_CI
SI_boot_results$UB_CI
summary(SI_boot_results)

#-----Multiple Imputation-------------------------------------------------------
set.seed(2025) # or result <- future_lapply(1:5, function(i) rnorm(1), future.seed = 2025 or TRUE)
multiple_imp <- MImpute( x = imp_result$imp_data,
                         w = kurt_data$lk,
                         outlier_matrix = out_result$outliers,
                         M = 5, k = 10)

#--3---MI-----------------------------------------------------------------------
MI_results <- thresh_MI(RD_org_obj = RD_org_obj
                 , imp_datasets = multiple_imp$imp_datasets
                 , alpha = 0.01)

MI_results$thresholds        # vector of 99th percentiles (length M)
MI_results$voted_outliers    # logical vector: TRUE if outlier in > M/2 imputations
summary(MI_results)

#--4---MI Boot-----------------------------------------------------------------------
MI_boot_results <- thresh_MI_boot(
  RD_org_obj   = RD_org_obj,                  # Output from RD(hk_data, mode = "auto")
  imp_datasets = multiple_imp$imp_datasets,   # Multiply imputed datasets (from MImpute)
  B            = 500,                         # Number of bootstrap samples per imputation
  alpha        = 0.01,                        # 99th percentile cutoff
  boot_quant   = 0.95,                        # 95% CI -> use 5% quantile as lower bound
  verbose      = TRUE                         # Print progress per imputation
)

# Print the final threshold
MI_boot_results$final_threshold

# See which time points were flagged as outliers
which(MI_boot_results$flagged_outliers)
summary(MI_boot_results)

#----Hardin and Rocke-----------------------------------------------------------
HR_result <- thresh_F(Q = ncol(hk_data), n = nrow(hk_data), h = RD_org_obj$h, quantile = 0.01)

HR_result$threshold
summary(HR_result)


###################
## NEW WRAPPER   #
##################
data_matrix <- fMRIscrub::Dat1
kurt_data <- ICA_extract_kurt(time_series = data_matrix)
hk_data <- kurt_data$hk
all_results <- threshold_RD(x = hk_data,
                            w = kurt_data$lk,
                            threshold_method = "all",
                            M = 5,
                            k = 10,
                            B = 500,
                            alpha = 0.01,
                            boot_quant = 0.95,
                            verbose = TRUE,
                            cutoff = 4,
                            trans = "SHASH",
                            quantile = 0.01)
