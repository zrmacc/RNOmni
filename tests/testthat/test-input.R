test_that("Test input check.", {
  
  n <- 3
  y <- seq_len(n)
  g <- array(0, dim = c(n, 2))
  x <- array(0, dim = c(n, 3))
  
  # No formatting errors.
  expect_error(BasicInputChecks(y, g, x), NA)
  
  # y contains NA.
  y_na <- y
  y_na[1] <- NA
  expect_error(BasicInputChecks(y_na, g, x))
  
  # X contains NA.
  x_na <- x
  x_na[1, 1] <- NA
  expect_error(BasicInputChecks(y, g, x_na))
  
  # X is a data.frame.
  x_df <- data.frame(x)
  expect_error(BasicInputChecks(y, g, x_df))
  
})
