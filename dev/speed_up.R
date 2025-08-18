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

# Or interactive flamegraph:
# install.packages("profvis")
# profvis::profvis({ ... your call ... })
