#' Dots parameter documentation
#' @param ... Additional arguments to to method-specific functions.
#' @name ...
NULL

#' Alpha parameter documentation
#' @param alpha Significance level used to compute RD threshold (default = 0.01 for 99th percentile).
#' @name alpha
NULL

#' B parameter documentation
#' @param B Integer; number of bootstrap samples per imputed dataset (default = 1000).
#' @name B
NULL

#' Boot_quant parameter documentation
#' @param boot_quant Numeric; confidence level for bootstrap confidence intervals (default = 0.95, for 95% CI).
#' @name boot_quant
NULL

#' Cov_mcd parameter documentation
#' @param cov_mcd Optional covariance matrix (p × p); required in "manual" mode.
#' @name cov_mcd
NULL

#' Cutoff parameter documentation
#' @param cutoff A numeric value indicating how many MADs away from the median to flag as outliers. The default value is set to be 4.
#' @name cutoff
NULL

#' Dist parameter documentation
#' @param dist Logical; if TRUE, compute squared robust Mahalanobis distances for all observations.
#' @name dist
NULL

#' Imp_data parameter documentation
#' @param imp_data A numeric matrix (T × p) of single-imputed data.
#' @name imp_data
NULL

#' Imp_datasets parameter documentation
#' @param imp_datasets A list of M numeric matrices (T × p); multiply imputed datasets.
#' @name imp_datasets
NULL

#' Impute_method parameter documentation
#' @param impute_method Character string; imputation method for univariate outliers.
#' @name impute_method
NULL

#' Ind_incld parameter documentation
#' @param ind_incld Optional vector of row indices used to compute the robust mean; required in "manual" mode.
#' @name ind_incld
NULL

#' K parameter documentation
#' @param k Integer; number of perturbation cycles per imputation (default = 10).
#' @name k
NULL

#' Quantile parameter documentation
#' @param quantile Numeric in (0,1) specifying the upper quantile for thresholding; the expected False Positive Rate for the chosen threshold.
#' @name quantile
NULL

#' M parameter documentation
#' @param M Integer; number of multiply imputed datasets (default = 5).
#' @name M
NULL

#' Mode parameter documentation
#' @param mode Character string; either "auto" (default) to compute MCD internally or "manual" to use user-supplied values.
#' @name mode
NULL

#' RD_obj parameter documentation
#' @param RD_obj Pre-computed RD_result object from \code{\link{compute_RD}}.
#' @name RD_obj
NULL

#' RD_org_obj parameter documentation
#' @param RD_org_obj Output list from \code{\link{compute_RD}} on the original data. Must contain $RD, $S_star, and $ind_incld.
#' @name RD_org_obj
NULL

#' Symmetric parameter documentation
#' @param symmetric Logical. Use symmetric bounds for empirical rule.
#' @name symmetric
NULL

#' Thr parameter documentation
#' @param thr Threshold multiplier for outlier detection (default = 4).
#' @name thr
NULL

#' Threshold_method parameter documentation
#' @param threshold_method Character string; one of "all","SI","SI_boot","MI","MI_boot","F", "SHASH".
#' @name threshold_method
NULL

#' Trans parameter documentation
#' @param trans Character string; transformation method, one of "SHASH" or "robZ".
#' @name trans
NULL

#' Upper_only parameter documentation
#' @param upper_only Logical. If TRUE, use right-tail only.
#' @name upper_only
NULL

#' Use_huber parameter documentation
#' @param use_huber Logical. If TRUE, apply Huber-like soft rejection.
#' @name use_huber
NULL

#' Verbose parameter documentation
#' @param verbose Logical; if TRUE, print progress messages.
#' @name verbose
NULL

#' W parameter documentation
#' @param w A numeric matrix (n_time × L) of low-kurtosis ICA components used as predictors (optional).
#' @name w
NULL

#' X parameter documentation
#' @param x A numeric matrix or data frame of dimensions T × p (observations × variables).
#' @name x
NULL
