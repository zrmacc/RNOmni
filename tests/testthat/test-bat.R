test_that("Test basic association test.", {
  
  # Run checks.
  n <- 5
  withr::local_seed(100)
  y <- stats::rnorm(n)
  x <- matrix(rep(1, n))
  g <- matrix(stats::rbinom(n, 2, 0.25))
  g_null <- array(0, dim = c(n, 1))
  
  expect_error(BAT(y, g, x), NA)
  expect_error(BAT(y, g, NULL), NA)
  expect_error(BAT(y, g_null), NA)
  
  # Type 1 error.
  n <- 1e3
  withr::local_seed(101)
  x <- cbind(1, rnorm(n))
  g <- replicate(1e3, rbinom(n = n, size = 2, prob = 0.25))
  storage.mode(g) <- "numeric"
  y <- as.numeric(x %*% c(1, 1)) + stats::rnorm(n)
  p_score <- BAT(y = y, G = g, X = x, test = "Score", simple = TRUE)
  p_wald <- BAT(y = y, G = g, X = x, test = "Wald", simple = TRUE)
  
  # Approximately correct type 1 error.
  t1e_score <- mean(p_score <= 0.05)
  t1e_wald <- mean(p_wald <= 0.05)
  expect_equal(t1e_score, 0.05, tolerance = 0.05)
  expect_equal(t1e_wald, 0.05, tolerance = 0.05)
  
})
