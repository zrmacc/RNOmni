test_that("input checks pass for valid inputs", {
  n <- 3
  y <- seq_len(n)
  g <- array(0, dim = c(n, 2))
  x <- array(0, dim = c(n, 3))
  expect_error(RNOmni:::BasicInputChecks(y, g, x), NA)
})

test_that("input checks fail when y is invalid", {
  n <- 3
  y <- seq_len(n)
  g <- array(0, dim = c(n, 2))
  x <- array(0, dim = c(n, 3))
  # y contains NA
  y_na <- y
  y_na[1] <- NA
  expect_error(RNOmni:::BasicInputChecks(y_na, g, x))
  # y is a matrix (must be vector)
  y_mat <- matrix(y, ncol = 1)
  expect_error(RNOmni:::BasicInputChecks(y_mat, g, x))
})

test_that("input checks fail when X is invalid", {
  n <- 3
  y <- seq_len(n)
  g <- array(0, dim = c(n, 2))
  x <- array(0, dim = c(n, 3))
  # X contains NA
  x_na <- x
  x_na[1, 1] <- NA
  expect_error(RNOmni:::BasicInputChecks(y, g, x_na))
  # X is a data.frame
  x_df <- data.frame(x)
  expect_error(RNOmni:::BasicInputChecks(y, g, x_df))
})

test_that("input checks fail when G is not a matrix", {
  n <- 3
  y <- seq_len(n)
  x <- array(0, dim = c(n, 3))
  g_df <- data.frame(g1 = rep(0, n), g2 = rep(0, n))
  expect_error(RNOmni:::BasicInputChecks(y, g_df, x))
})
