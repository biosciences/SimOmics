#' Simulate Latent Factors
#'
#' Generate latent variables (shared or block-specific) to represent hidden biological variation.
#' These can be used to drive signal across omics blocks (e.g., transcriptome, proteome).
#'
#' @param n Integer. Number of samples.
#' @param n_factors Integer. Number of latent factors.
#' @param block_dims Named list. Each element is an integer representing the number of features in each omics block (e.g., list(X1 = 500, X2 = 200)).
#' @param shared Logical. If TRUE, the same latent factors affect all blocks. If FALSE, each block gets its own independent latent matrix.
#' @param loadings_sd Numeric. Standard deviation of the random loadings matrix.
#' @param seed Optional integer. Random seed for reproducibility.
#'
#' @return A list with elements:
#' \describe{
#'   \item{Z}{Latent variable matrix (n x n_factors)}
#'   \item{X_blocks}{List of simulated omics blocks, each matrix of dimension n x p}
#'   \item{W}{List of loading matrices, each of size n_factors x p}
#' }
#'
#' @examples
#' block_dims <- list(transcriptome = 1000, proteome = 200)
#' result <- simulate_latent_factors(n = 100, n_factors = 3, block_dims = block_dims)
#'
#' @export
simulate_latent_factors <- function(n, n_factors, block_dims, shared = TRUE, loadings_sd = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  if (!is.list(block_dims) || is.null(names(block_dims))) {
    stop("block_dims must be a named list of feature counts for each omics block.")
  }
  
  # Generate latent variable matrix Z (n x n_factors)
  Z <- matrix(rnorm(n * n_factors), nrow = n, ncol = n_factors)
  
  X_blocks <- list()
  W_list <- list()
  
  for (block_name in names(block_dims)) {
    p <- block_dims[[block_name]]
    
    # Random loadings matrix W: n_factors x p
    W <- matrix(rnorm(n_factors * p, sd = loadings_sd), nrow = n_factors, ncol = p)
    
    # Block-specific latent variables if not shared
    Z_block <- if (shared) Z else matrix(rnorm(n * n_factors), nrow = n)
    
    # Simulate data block: X = Z * W + noise
    X <- Z_block %*% W
    
    X_blocks[[block_name]] <- X
    W_list[[block_name]] <- W
  }
  
  return(list(Z = Z, X_blocks = X_blocks, W = W_list))
}