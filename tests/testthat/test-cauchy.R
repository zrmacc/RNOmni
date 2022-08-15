test_that("Test Cauchy combination.", {
  
  # Only 1 p-value: should remain unchanged.
  p <- seq(from = 0.1, to = 0.9, by = 0.1)
  obs <- sapply(p, OmniP)
  expect_equal(obs, p)
  
  # Repeated p-value: should remain unchanged.
  p <- c(0.05, 0.05, 0.05)
  obs <- OmniP(p)
  expect_equal(obs, 0.05)
  
  obs <- OmniP(p, w = c(1, 2, 1))  # Weighting should not affect this.
  expect_equal(obs, 0.05)
  
  # Check case of uniform weights.
  p <- c(0.1, 0.2, 0.3)
  obs1 <- OmniP(p)
  obs2 <- OmniP(p, w = c(1, 1, 1))
  obs3 <- OmniP(p, w = c(2, 2, 2))
  expect_true(all.equal(obs1, obs2, obs3))
  
  # Check non-uniform weights.
  p <- c(0.1, 0.05)
  obs <- OmniP(p, w = c(1, 2))
  expect_lt(obs, 0.1)
  
  obs <- OmniP(p, w = c(2, 1))
  expect_gt(obs, 0.05)
  
})
