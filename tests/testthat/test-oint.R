test_that("OINT runs with X=NULL (marginal test)", {
  n <- 15
  withr::local_seed(108)
  y <- stats::rnorm(n)
  G <- matrix(stats::rbinom(n * 2, 2, 0.25), nrow = n, ncol = 2)
  expect_error(OINT(y = y, G = G, X = NULL, simple = TRUE), NA)
  p <- OINT(y = y, G = G, X = NULL, simple = TRUE)
  expect_length(p, 2)
  expect_true(all(p >= 0 & p <= 1))
})

test_that("OINT simple=FALSE returns DINT-p, IINT-p, OINT-p columns", {
  n <- 20
  withr::local_seed(109)
  X <- cbind(1, stats::rnorm(n))
  G <- matrix(stats::rbinom(n * 2, 2, 0.25), nrow = n, ncol = 2)
  y <- exp(as.numeric(X %*% c(1, 1)) + stats::rnorm(n))
  out <- OINT(y = y, G = G, X = X, simple = FALSE)
  expect_true(is.matrix(out))
  expect_equal(nrow(out), 2)
  expect_equal(colnames(out), c("DINT-p", "IINT-p", "OINT-p"))
  expect_true(all(out >= 0 & out <= 1))
})

test_that("OINT accepts custom weights", {
  n <- 30
  withr::local_seed(110)
  X <- cbind(1, stats::rnorm(n))
  G <- matrix(stats::rbinom(n * 2, 2, 0.25), nrow = n, ncol = 2)
  y <- as.numeric(X %*% c(1, 0)) + stats::rnorm(n)
  p_default <- OINT(y = y, G = G, X = X, weights = c(1, 1), simple = TRUE)
  p_custom <- OINT(y = y, G = G, X = X, weights = c(2, 1), simple = TRUE)
  expect_length(p_default, 2)
  expect_length(p_custom, 2)
})

test_that("OINT type I error approximately 0.05 under null", {
  n <- 1e3
  withr::local_seed(101)
  x <- cbind(1, rnorm(n))
  g <- replicate(1e3, rbinom(n = n, size = 2, prob = 0.25))
  storage.mode(g) <- "numeric"
  y <- as.numeric(x %*% c(1, 1)) + stats::rnorm(n)
  p <- OINT(y = y, G = g, X = x, simple = TRUE)
  expect_equal(mean(p <= 0.05), 0.05, tolerance = 0.05)
})
