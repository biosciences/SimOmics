#' Simulate Multivariate Data
#'
#' Generate multivariate normal data using a specified covariance matrix or identity.
#'
#' @param n Integer. Number of samples.
#' @param p Integer. Number of variables (features).
#' @param Sigma Optional covariance matrix (p x p). If NULL, identity matrix is used.
#' @param mean Optional numeric vector of length p. Defaults to zero vector.
#' @param seed Optional integer. Seed for reproducibility.
#'
#' @return A matrix of dimension n x p, drawn from MVN(mean, Sigma).
#'
#' @examples
#' set.seed(42)
#' sim_data <- simulate_multivariate(n = 100, p = 10)
#'
#' # With custom covariance
#' Sigma <- diag(1, 5); Sigma[1, 2] <- Sigma[2, 1] <- 0.5
#' sim_data2 <- simulate_multivariate(n = 50, p = 5, Sigma = Sigma)
#'
#' @export
simulate_multivariate <- function(n, p, Sigma = NULL, mean = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  if (is.null(Sigma)) {
    Sigma <- diag(1, p)
  } else {
    if (!is.matrix(Sigma) || ncol(Sigma) != p || nrow(Sigma) != p) {
      stop("Sigma must be a p x p matrix.")
    }
  }
  
  if (is.null(mean)) {
    mean <- rep(0, p)
  } else {
    if (length(mean) != p) {
      stop("Mean vector must be of length p.")
    }
  }
  
  # Simulate MVN samples
  mvtnorm::rmvnorm(n = n, mean = mean, sigma = Sigma)
}