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

# Use smaller subset for faster tests
data_matrix <- fMRIscrub::Dat1[1:50, 1:100]  # Smaller for speed

# Run pipeline up to the point we need
kurt_data <- ICA_extract_kurt(time_series = data_matrix)
hk_data <- kurt_data$hk
out_result <- univOut(hk_data = hk_data, cutoff = 4, trans = "SHASH")
imp_result <- impTemp_univOut(x = hk_data, outlier_mask = out_result$outliers)
RD_org_obj <- comp_RD(data_matrix = hk_data, mode = "auto")

# Save the setup data needed for all tests
saveRDS(list(
  hk_data = hk_data,
  kurt_data = kurt_data,
  out_result = out_result,
  imp_result = imp_result,
  RD_org_obj = RD_org_obj
), "tests/testthat/fixtures/test_setup_data.rds")

# Test 1: SI method
SI_results <- SI(RD_org_obj = RD_org_obj,
                 imp_data = imp_result$imp_data,
                 alpha = 0.01)
saveRDS(SI_results, "tests/testthat/fixtures/SI_reference.rds")

# Test 2: SI_boot method
SI_boot_results <- SI_boot(RD_org_obj = RD_org_obj,
                           imp_data = imp_result$imp_data,
                           B = 50, alpha = 0.01, boot_quant = 0.95,  # Smaller B for speed
                           verbose = FALSE)
saveRDS(SI_boot_results, "tests/testthat/fixtures/SI_boot_reference.rds")

# Test 3: MI method (with smaller M for speed)
multiple_imp <- MImpute(x = imp_result$imp_data,
                        w = kurt_data$lk,
                        outlier_matrix = out_result$outliers,
                        M = 3, k = 5)  # Smaller M, k for speed
MI_results <- MI(RD_org_obj = RD_org_obj,
                 imp_datasets = multiple_imp$imp_datasets,
                 alpha = 0.01)
saveRDS(list(multiple_imp = multiple_imp, MI_results = MI_results),
        "tests/testthat/fixtures/MI_reference.rds")

# Test 4: Hardin-Rocke method
HR_result <- Fit_F(Q = ncol(hk_data), n = nrow(hk_data),
                   h = RD_org_obj$h, quantile = 0.01)
saveRDS(HR_result, "tests/testthat/fixtures/HR_reference.rds")

cat("All reference data saved!\n")
