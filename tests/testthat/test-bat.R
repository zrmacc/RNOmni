library(RNOmni)

test_that("BAT runs without error for valid inputs", {
  n <- 5
  withr::local_seed(100)
  y <- stats::rnorm(n)
  x <- matrix(rep(1, n))
  g <- matrix(stats::rbinom(n, 2, 0.25))
  g_null <- array(0, dim = c(n, 1))
  expect_error(BAT(y, g, x), NA)
  expect_error(BAT(y, g, NULL), NA)
  expect_error(BAT(y, g_null), NA)
})

test_that("BAT rejects invalid test type", {
  n <- 10
  y <- stats::rnorm(n)
  G <- matrix(stats::rbinom(n, 2, 0.25))
  expect_error(BAT(y, G, test = "LR"))
  expect_error(BAT(y, G, test = "invalid"))
})

test_that("BAT simple=TRUE returns named vector of p-values", {
  n <- 20
  withr::local_seed(102)
  X <- cbind(1, stats::rnorm(n))
  G <- matrix(stats::rbinom(n * 3, 2, 0.25), nrow = n, ncol = 3)
  colnames(G) <- c("snp1", "snp2", "snp3")
  y <- as.numeric(X %*% c(1, 0)) + stats::rnorm(n)
  p <- BAT(y = y, G = G, X = X, simple = TRUE)
  expect_equal(length(p), 3)
  expect_named(p, c("snp1", "snp2", "snp3"))
  expect_true(all(p >= 0 & p <= 1))
})

test_that("BAT simple=FALSE returns matrix with Score/Wald, SE, Z, P", {
  n <- 20
  withr::local_seed(103)
  X <- cbind(1, stats::rnorm(n))
  G <- matrix(stats::rbinom(n * 2, 2, 0.25), nrow = n, ncol = 2)
  y <- as.numeric(X %*% c(1, 0)) + stats::rnorm(n)
  out <- BAT(y = y, G = G, X = X, test = "Score", simple = FALSE)
  expect_true(is.matrix(out))
  expect_equal(nrow(out), 2)
  expect_equal(colnames(out), c("Score", "SE", "Z", "P"))
  out_wald <- BAT(y = y, G = G, X = X, test = "Wald", simple = FALSE)
  expect_equal(colnames(out_wald), c("Wald", "SE", "Z", "P"))
})

test_that("BAT type I error approximately 0.05 under null", {
  n <- 1e3
  withr::local_seed(101)
  x <- cbind(1, rnorm(n))
  g <- replicate(1e3, rbinom(n = n, size = 2, prob = 0.25))
  storage.mode(g) <- "numeric"
  y <- as.numeric(x %*% c(1, 1)) + stats::rnorm(n)
  p_score <- BAT(y = y, G = g, X = x, test = "Score", simple = TRUE)
  p_wald <- BAT(y = y, G = g, X = x, test = "Wald", simple = TRUE)
  expect_equal(mean(p_score <= 0.05), 0.05, tolerance = 0.05)
  expect_equal(mean(p_wald <= 0.05), 0.05, tolerance = 0.05)
})
