test_that("simulate_multiomics returns proper structure", {
  sim <- simulate_multiomics(
    n = 50,
    block_dims = list(t1 = 100, p1 = 50),
    n_factors = 2,
    seed = 42
  )
  
  expect_type(sim, "list")
  expect_named(sim, c("X_blocks", "Z", "W", "params"))
  expect_true(all(sapply(sim$X_blocks, is.matrix)))
  expect_equal(nrow(sim$X_blocks[[1]]), 50)
})