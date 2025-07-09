#' Simulate Block Structure for Multi-Omics Covariance
#'
#' Generate block-wise covariance matrix with optional inter-block correlation.
#'
#' @param block_dims Named list. Each element gives number of features in each omics block.
#' @param block_corr Numeric (0-1). Constant correlation between blocks.
#' @param within_block_sd Numeric. Standard deviation within each block.
#' @param seed Optional integer for reproducibility.
#'
#' @return A symmetric positive-definite covariance matrix of size total_p x total_p.
#'
#' @export
simulate_block_structure <- function(block_dims, block_corr = 0.1, within_block_sd = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  if (!is.list(block_dims) || is.null(names(block_dims))) {
    stop("block_dims must be a named list of integers.")
  }
  
  block_names <- names(block_dims)
  block_sizes <- unlist(block_dims)
  total_p <- sum(block_sizes)
  
  Sigma <- matrix(0, nrow = total_p, ncol = total_p)
  rownames(Sigma) <- colnames(Sigma) <- unlist(
    lapply(seq_along(block_names), function(i) {
      paste0(block_names[i], "_", seq_len(block_sizes[i]))
    })
  )
  
  # 1. Fill within-block variances (diagonal blocks)
  current_index <- 1
  block_indices <- list()
  
  for (i in seq_along(block_sizes)) {
    p <- block_sizes[i]
    idx <- current_index:(current_index + p - 1)
    Sigma[idx, idx] <- diag(within_block_sd^2, p)
    block_indices[[block_names[i]]] <- idx
    current_index <- current_index + p
  }
  
  # 2. Fill between-block correlations (off-diagonal blocks)
  if (block_corr != 0) {
    for (i in 1:(length(block_names) - 1)) {
      for (j in (i + 1):length(block_names)) {
        idx_i <- block_indices[[block_names[i]]]
        idx_j <- block_indices[[block_names[j]]]
        
        Sigma[idx_i, idx_j] <- block_corr * within_block_sd^2
        Sigma[idx_j, idx_i] <- t(Sigma[idx_i, idx_j])
      }
    }
  }
  
  # 3. Ensure positive definiteness
  eigenvalues <- eigen(Sigma, symmetric = TRUE)$values
  if (any(eigenvalues <= 0)) {
    Sigma <- Sigma + diag(abs(min(eigenvalues)) + 1e-6, total_p)
  }
  
  return(Sigma)
}