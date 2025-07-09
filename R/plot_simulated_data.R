#' Plot Simulated Data Structures
#'
#' Visualize simulated multi-omics data blocks or latent variables using heatmaps, PCA, or correlation plots.
#'
#' @param sim_data List. Output from simulate_latent_factors() or simulate_multiomics().
#' @param type Character. Type of plot: "heatmap", "pca", or "correlation".
#' @param block Character. Optional. Which block to plot (for heatmap or PCA). Defaults to the first block.
#' @param center Logical. Whether to center the data before PCA.
#' @param scale Logical. Whether to scale the data before PCA.
#'
#' @return A ggplot2 object or a base plot.
#'
#' @examples
#' sim <- simulate_latent_factors(n = 100, n_factors = 3, block_dims = list(X1 = 200, X2 = 100))
#' plot_simulated_data(sim, type = "pca", block = "X1")
#'
#' @export
plot_simulated_data <- function(sim_data, type = c("pca", "heatmap", "correlation"),
                                block = NULL, center = TRUE, scale = TRUE) {
  type <- match.arg(type)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2.")
  }
  
  if (type == "pca") {
    if (is.null(block)) block <- names(sim_data$X_blocks)[1]
    X <- sim_data$X_blocks[[block]]
    pca <- prcomp(X, center = center, scale. = scale)
    df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2])
    
    ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2)) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::labs(title = paste("PCA of block:", block), x = "PC1", y = "PC2") +
      ggplot2::theme_minimal()
    
  } else if (type == "heatmap") {
    if (is.null(block)) block <- names(sim_data$X_blocks)[1]
    X <- sim_data$X_blocks[[block]]
    mat <- scale(X)
    heatmap(mat, Rowv = NA, Colv = NA, scale = "none",
            main = paste("Heatmap of block:", block),
            xlab = "Features", ylab = "Samples")
    return(invisible(NULL))
    
  } else if (type == "correlation") {
    all_blocks <- sim_data$X_blocks
    if (length(all_blocks) < 2) stop("Need at least two blocks for correlation plot.")
    
    merged <- do.call(cbind, all_blocks)
    corrmat <- cor(merged)
    
    if (!requireNamespace("corrplot", quietly = TRUE)) {
      stop("Please install the 'corrplot' package.")
    }
    
    corrplot::corrplot(corrmat, method = "color", order = "hclust",
                       tl.cex = 0.6, main = "Correlation matrix across all omics blocks")
    return(invisible(NULL))
  }
}