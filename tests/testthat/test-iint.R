test_that("Test basic association test.", {

  # Type 1 error.
  n <- 1e3
  withr::local_seed(101)
  x <- cbind(1, rnorm(n))
  g <- replicate(1e3, rbinom(n = n, size = 2, prob = 0.25))
  storage.mode(g) <- "numeric"
  y <- as.numeric(x %*% c(1, 1)) + stats::rnorm(n)
  p <- IINT(y = y, G = g, X = x, simple = TRUE)
  
  # Approximately correct type 1 error.
  t1e <- mean(p <= 0.05)
  expect_equal(t1e, 0.05, tolerance = 0.05)
  
})
