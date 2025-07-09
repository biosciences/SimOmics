#' Simulate Multi-Omics Dataset
#'
#' Generate a synthetic multi-omics dataset driven by latent factors,
#' with optional block-wise structure, correlation, and Gaussian noise.
#'
#' @param n Integer. Number of samples.
#' @param block_dims Named list. Number of features in each omics block (e.g., list(transcriptome = 1000, proteome = 200)).
#' @param n_factors Integer. Number of latent factors.
#' @param shared Logical. Whether latent factors are shared across all blocks. Defaults to TRUE.
#' @param noise_sd Numeric. Standard deviation of Gaussian noise added to each block. Default is 1.
#' @param block_corr Optional numeric. Inter-block correlation (used to build block covariance matrix).
#' @param seed Optional integer. Random seed for reproducibility.
#'
#' @return A list containing:
#' \describe{
#'   \item{X_blocks}{List of omics data matrices (each n x p)}
#'   \item{Z}{Latent factor matrix (n x n_factors)}
#'   \item{W}{List of loading matrices}
#'   \item{params}{List of parameters used to generate data}
#' }
#'
#' @export
simulate_multiomics <- function(
  n,
  block_dims,
  n_factors,
  shared = TRUE,
  noise_sd = 1,
  block_corr = NULL,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  # Step 1: Simulate latent factor-driven signal
  latent_res <- simulate_latent_factors(
    n = n,
    n_factors = n_factors,
    block_dims = block_dims,
    shared = shared,
    loadings_sd = 1,
    seed = seed
  )
  
  X_blocks <- latent_res$X_blocks
  W <- latent_res$W
  Z <- latent_res$Z
  
  # Step 2: Add inter-block correlation if specified
  if (!is.null(block_corr)) {
    Sigma <- simulate_block_structure(block_dims, block_corr = block_corr)
    merged <- do.call(cbind, X_blocks)
    
    # Simulate new block-structured correlated data
    correlated <- MASS::mvrnorm(n = n, mu = rep(0, ncol(merged)), Sigma = Sigma)
    
    # Split back into omics blocks
    block_lengths <- unlist(block_dims)
    end_points <- cumsum(block_lengths)
    start_points <- c(1, head(end_points + 1, -1))
    
    X_blocks <- mapply(function(start, end) {
      as.matrix(correlated[, start:end])
    }, start_points, end_points, SIMPLIFY = FALSE)
    
    names(X_blocks) <- names(block_dims)
  }
  
  # Step 3: Add Gaussian noise to each block
  X_blocks <- lapply(X_blocks, function(X) {
    X + matrix(rnorm(n * ncol(X), sd = noise_sd), nrow = n)
  })
  
  return(list(
    X_blocks = X_blocks,
    Z = Z,
    W = W,
    params = list(
      n = n,
      block_dims = block_dims,
      n_factors = n_factors,
      shared = shared,
      noise_sd = noise_sd,
      block_corr = block_corr,
      seed = seed
    )
  ))
}