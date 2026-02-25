test_that("DINT runs with X=NULL (marginal test)", {
  n <- 15
  withr::local_seed(104)
  y <- stats::rnorm(n)
  G <- matrix(stats::rbinom(n * 2, 2, 0.25), nrow = n, ncol = 2)
  expect_error(DINT(y = y, G = G, X = NULL, simple = TRUE), NA)
  p <- DINT(y = y, G = G, X = NULL, simple = TRUE)
  expect_length(p, 2)
  expect_true(all(p >= 0 & p <= 1))
})

test_that("DINT simple=FALSE returns matrix with statistic, SE, Z, P", {
  n <- 20
  withr::local_seed(105)
  X <- cbind(1, stats::rnorm(n))
  G <- matrix(stats::rbinom(n * 2, 2, 0.25), nrow = n, ncol = 2)
  y <- exp(as.numeric(X %*% c(1, 1)) + stats::rnorm(n))
  out <- DINT(y = y, G = G, X = X, test = "Score", simple = FALSE)
  expect_true(is.matrix(out))
  expect_equal(nrow(out), 2)
  expect_equal(colnames(out), c("Score", "SE", "Z", "P"))
})

test_that("DINT type I error approximately 0.05 under null", {
  n <- 1e3
  withr::local_seed(101)
  x <- cbind(1, rnorm(n))
  g <- replicate(1e3, rbinom(n = n, size = 2, prob = 0.25))
  storage.mode(g) <- "numeric"
  y <- as.numeric(x %*% c(1, 1)) + stats::rnorm(n)
  p_score <- DINT(y = y, G = g, X = x, test = "Score", simple = TRUE)
  p_wald <- DINT(y = y, G = g, X = x, test = "Wald", simple = TRUE)
  expect_equal(mean(p_score <= 0.05), 0.05, tolerance = 0.05)
  expect_equal(mean(p_wald <= 0.05), 0.05, tolerance = 0.05)
})
