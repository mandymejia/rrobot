library(fMRIscrub)
data_matrix <- fMRIscrub::Dat1
# source("point to scripts read the funct)"
kurt_data <- ICA_extract_kurt(time_series = data_matrix) # use fMRIscrub
hk_abide1 <- kurt_data$hk
lk_abide1 <- kurt_data$lk
dim(hk_abide1)

# want to test if we need to save RDs for each methods:

t <- system.time({
  hk1_obj <- rrobot::compute_RD(hk_abide1, mode = "auto")
  hk1_thrs <- rrobot::threshold_RD(
                  x = hk_abide1, w = lk_abide1,
                  method = "all", RD_obj = hk1_obj,
                  impute_method = "interp",
                  M = 50, k = 5, B = 1000)
})

print(t)              # user / system / elapsed
t["elapsed"]          # numeric seconds


# Base R profiler
Rprof("profile.out")
hk1_thrs <- rrobot::threshold_RD(x = hk_abide1, w = lk_abide1,
                                 method = "all", RD_obj = hk1_obj,
                                 impute_method = "interp",
                                 M = 50, k = 5, B = 1000)
Rprof(NULL)
summaryRprof("profile.out")  # shows top functions by time



out_result <- rrobot::univOut(x = hk_abide1, cutoff = 4, method = "SHASH")

MImpute_obj <- rrobot::MImpute(x = hk_abide1, w = lk_abide1,
                               outlier_matrix = out_result$outliers,
                               k = 5, M = 50)

MImpute_fast_obj <- MImpute_fast(x = hk_abide1, w = lk_abide1,
                               outlier_matrix = out_result$outliers,
                               k = 5, M = 50)


##---- speed benchmark---------------------------------------------------------
set.seed(2025)

# Quick wall-clock timing
t_old  <- system.time({
  MImpute_obj <- rrobot::MImpute(
    x = hk_abide1, w = lk_abide1,
    outlier_matrix = out_result$outliers,
    k = 5, M = 50
  )
})

t_fast <- system.time({
  MImpute_fast_obj <- MImpute_fast(
    x = hk_abide1, w = lk_abide1,
    outlier_matrix = out_result$outliers,
    k = 5, M = 50
  )
})

print(t_old["elapsed"]); print(t_fast["elapsed"])

# Optional: richer benchmark
# install.packages("bench")
library(bench)
bench::mark(
  MImpute      = rrobot::MImpute(hk_abide1, lk_abide1, out_result$outliers, k = 5, M = 50),
  MImpute_fast = MImpute_fast   (hk_abide1, lk_abide1, out_result$outliers, k = 5, M = 50),
  iterations = 3, check = FALSE
)
