#' @importFrom ggplot2 ggplot aes geom_histogram geom_line geom_vline geom_text geom_point
#' @importFrom ggplot2 scale_fill_manual scale_color_manual scale_x_log10 theme_minimal theme
#' @importFrom ggplot2 element_text element_blank after_stat ylab xlab labs geom_tile
#' @importFrom stats density df
#' @importFrom reshape2 melt
#' @importFrom rlang .data
NULL

# ==========================================================================
# Global Variables Declaration for R CMD Check
# ==========================================================================
#
# This file suppresses R CMD check NOTEs about "no visible binding for global variable"
# for variables used in ggplot2 non-standard evaluation (NSE) contexts.
#
# When using ggplot2's aes() function, variable names are evaluated using NSE,
# meaning they don't need to be explicitly defined as objects in the global
# environment. However, R CMD check doesn't understand this and flags these
# as potential problems.
#
# The variables listed below are:
# - Column names from data frames used in ggplot2 aesthetics
# - Variables created through data manipulation (e.g., reshape2::melt)
# - Computed variables in ggplot2 (e.g., after_stat(density))
#
# This is a standard practice recommended by:
# - Writing R Extensions manual
# - ggplot2 documentation
# - Most major R packages on CRAN
#
# References:
# - https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Registering-S3-methods
# - https://ggplot2.tidyverse.org/articles/ggplot2-in-packages.html
# ==========================================================================

if(getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    # Variables from plot_RD_histogram_multi() and related functions
    "thr_log",        # Log-transformed threshold values
    "method",         # Threshold method names (SI, MI, etc.)
    "z",              # Log-transformed robust distances
    "label",          # Observation inclusion/exclusion labels

    # Variables from plot_univOut() - outlier heatmap visualization
    "Time",           # Time/observation index
    "Variable",       # Variable/component index
    "Outlier",        # Logical outlier indicator

    # Variables from plot_imputations() - multiple imputation visualization
    "Original",       # Original time series values
    "Temporal",       # Temporally imputed values
    "Imputed",        # Multiple imputation values

    # Variables from plot_RD_histogram() and plot_F_histogram()
    "distance",       # Robust distance values
    "observation",    # Included/excluded observation status
    "x",              # Generic x-axis values
    "y",              # Generic y-axis values (densities, etc.)

    # Generic plotting variables
    ".data"           # rlang pronoun (if using .data$ syntax)
  ))
}

