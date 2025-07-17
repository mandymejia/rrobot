#' Extract high and low kurtosis ICA components using fMRIscrub
#'
#' This function performs ICA projection on a masked fMRI dataset using `pscrub()` from the
#' fMRIscrub package and separates high- and low-kurtosis components based on the
#' `highkurt` indicator returned by the method.
#'
#' @param data_matrix A numeric matrix of dimension T × V (time × voxels), with one time series per voxel.
#' @param mask_path File path to a binary brain mask in NIfTI format (.nii or .nii.gz).
#' @param save_prefix (Optional) Prefix for the saved `.RData` files (e.g., "abide1"). Required if `save_dir` is specified.
#' @param save_dir (Optional) Directory to save `hk_*.RData` and `lk_*.RData` matrices. If not specified, nothing is saved.
#'
#' @return A list with elements:
#' \describe{
#'   \item{hk}{High-kurtosis ICA component matrix (T × K)}
#'   \item{lk}{Low-kurtosis ICA component matrix (T × L)}
#'   \item{highkurt}{Logical vector indicating which components are high-kurtosis}
#' }
#'
#' @importFrom fMRIscrub pscrub
#' @importFrom oro.nifti readNIfTI
#' @export
#'
#' @examples
#' \dontrun{
#' result <- ICA_extract_kurt(
#'   data_matrix = fMRIscrub::Dat1,
#'   mask_path = system.file("extdata", "Dat1_mask.nii.gz", package = "fMRIscrub"),
#'   save_prefix = "abide1",
#'   save_dir = "data"
#' )
#' }

ICA_extract_kurt <- function(data_matrix, mask_path, save_prefix = NULL, save_dir = NULL) {
  stopifnot(is.matrix(data_matrix), file.exists(mask_path))

  # Load and apply mask
  mask_img <- oro.nifti::readNIfTI(mask_path, reorient = FALSE)
  mask_vec <- as.logical(mask_img)

  if (ncol(data_matrix) != sum(mask_vec)) {
    stop("Number of voxels in mask does not match number of columns in data_matrix.")
  }

  # Run pscrub with ICA
  message("Running pscrub ICA projection...")
  ICA_result <- fMRIscrub::pscrub(data_matrix, projection = c("ICA"))

  M <- ICA_result$ICA$M                 # Time series matrix (T × num_components)
  highkurt_idx <- ICA_result$ICA$highkurt  # Logical vector

  # Subset components
  hk <- M[, highkurt_idx, drop = FALSE]
  lk <- M[, !highkurt_idx, drop = FALSE]

  # Save if requested
  if (!is.null(save_prefix) && !is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

    save(hk, file = file.path(save_dir, paste0("hk_", save_prefix, "_ICA.RData")))
    save(lk, file = file.path(save_dir, paste0("lk_", save_prefix, "_ICA.RData")))

    message("Saved to: ", save_dir)
  }

  return(list(
    hk = hk,
    lk = lk,
    highkurt = highkurt_idx
  ))
}
