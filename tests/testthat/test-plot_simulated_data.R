test_that("plot_simulated_data runs PCA without error", {
  sim <- simulate_multiomics(n = 30, block_dims = list(x = 50), n_factors = 2)
  p <- plot_simulated_data(sim, type = "pca", block = "x")
  expect_s3_class(p, "ggplot")
})