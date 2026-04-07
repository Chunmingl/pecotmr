#' Create an LD loader for on-demand block-wise LD retrieval
#'
#' Constructs a loader function that retrieves per-block LD matrices on
#' demand. This avoids loading all blocks into memory simultaneously,
#' which is critical for genome-wide analyses with hundreds of blocks.
#'
#' Three modes are supported:
#'
#' \describe{
#'   \item{list mode (R)}{Pre-loaded list of LD correlation matrices.
#'     Simple but uses more memory. Set \code{R_list}.}
#'   \item{list mode (X)}{Pre-loaded list of genotype matrices (n x p_g).
#'     Set \code{X_list}.}
#'   \item{file mode}{Loads LD from a pecotmr metadata TSV file on the fly
#'     via \code{\link{load_LD_matrix}}. Memory-efficient for large datasets.
#'     Set \code{ld_meta_path} and \code{regions}.}
#' }
#'
#' @param R_list List of G precomputed LD correlation matrices (p_g x p_g).
#'   Mutually exclusive with \code{X_list} and \code{ld_meta_path}.
#' @param X_list List of G genotype matrices (n x p_g). Mutually exclusive
#'   with \code{R_list} and \code{ld_meta_path}.
#' @param ld_meta_path Path to a pecotmr LD metadata TSV file (as used by
#'   \code{\link{load_LD_matrix}}). Mutually exclusive with \code{R_list}
#'   and \code{X_list}.
#' @param regions Character vector of G region strings (e.g.,
#'   \code{"chr22:17238266-19744294"}). Required when \code{ld_meta_path}
#'   is used.
#' @param return_genotype Logical. When using file mode, return the
#'   genotype matrix X (\code{TRUE}) or LD correlation R (\code{FALSE},
#'   default). Returning X enables low-rank RSS paths in downstream methods.
#' @param max_variants Integer or \code{NULL}. If set, randomly subsample
#'   blocks larger than this to control memory usage.
#'
#' @return A function \code{loader(g)} that, given a block index \code{g},
#'   returns the corresponding LD matrix or genotype matrix.
#'
#' @examples
#' # List mode with pre-computed LD
#' R1 <- diag(10)
#' R2 <- diag(15)
#' loader <- ld_loader(R_list = list(R1, R2))
#' loader(1)  # returns R1
#' loader(2)  # returns R2
#'
#' @export
ld_loader <- function(R_list = NULL, X_list = NULL,
                      ld_meta_path = NULL, regions = NULL,
                      return_genotype = FALSE,
                      max_variants = NULL) {
  # Validate: exactly one source
  n_sources <- sum(!is.null(R_list), !is.null(X_list), !is.null(ld_meta_path))
  if (n_sources != 1)
    stop("Provide exactly one of R_list, X_list, or ld_meta_path.")

  if (!is.null(R_list)) {
    # List mode (R matrices)
    loader <- function(g) {
      R <- R_list[[g]]
      if (!is.null(max_variants) && ncol(R) > max_variants) {
        keep <- sort(sample(ncol(R), max_variants))
        R <- R[keep, keep]
      }
      R
    }
  } else if (!is.null(X_list)) {
    # List mode (genotype matrices)
    loader <- function(g) {
      X <- X_list[[g]]
      if (!is.null(max_variants) && ncol(X) > max_variants) {
        keep <- sort(sample(ncol(X), max_variants))
        X <- X[, keep]
      }
      X
    }
  } else {
    # File mode: load on the fly via load_LD_matrix()
    if (is.null(regions))
      stop("'regions' is required when using ld_meta_path.")

    loader <- function(g) {
      ld <- load_LD_matrix(ld_meta_path, region = regions[g],
                           return_genotype = return_genotype)
      mat <- ld$LD_matrix
      if (!is.null(max_variants) && ncol(mat) > max_variants) {
        keep <- sort(sample(ncol(mat), max_variants))
        if (return_genotype || nrow(mat) > ncol(mat)) {
          mat <- mat[, keep]
        } else {
          mat <- mat[keep, keep]
        }
      }
      # Center and scale genotype matrices
      if (return_genotype || nrow(mat) > ncol(mat)) {
        mat <- scale(mat)
        mat[is.na(mat)] <- 0
      }
      mat
    }
  }

  loader
}
