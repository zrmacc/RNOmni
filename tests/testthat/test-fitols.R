test_that("FitOLS returns list with Beta, V, Ibb, Resid", {
  n <- 10
  X <- cbind(1, stats::rnorm(n))
  y <- as.numeric(X %*% c(1, 0.5)) + stats::rnorm(n)
  fit <- FitOLS(y, X)
  expect_type(fit, "list")
  expect_named(fit, c("Beta", "V", "Ibb", "Resid"))
  expect_equal(length(fit$Beta), ncol(X))
  expect_length(fit$Resid, n)
  expect_equal(length(fit$V), 1)
  expect_equal(dim(fit$Ibb), c(ncol(X), ncol(X)))
})

test_that("FitOLS agrees with lm on coefficients and residuals", {
  withr::local_seed(42)
  n <- 50
  X <- cbind(1, stats::rnorm(n), stats::rnorm(n))
  y <- as.numeric(X %*% c(1, 0.3, -0.2)) + stats::rnorm(n)
  fit_ols <- FitOLS(y, X)
  fit_lm <- stats::lm(y ~ X - 1)
  expect_equal(as.numeric(fit_ols$Beta), as.numeric(stats::coef(fit_lm)), tolerance = 1e-10)
  expect_equal(as.numeric(fit_ols$Resid), as.numeric(stats::resid(fit_lm)), tolerance = 1e-10)
  # Residual variance: sigma^2 from lm
  sigma2_lm <- summary(fit_lm)$sigma^2
  expect_equal(fit_ols$V, sigma2_lm, tolerance = 1e-10)
})

test_that("FitOLS intercept-only gives mean and residual variance", {
  n <- 5
  y <- c(1, 2, 3, 4, 5)
  X <- matrix(1, nrow = n, ncol = 1)
  fit <- FitOLS(y, X)
  expect_equal(as.numeric(fit$Beta), mean(y))
  expect_equal(fit$V, stats::var(y), tolerance = 1e-10)
})
