test_that("OmniP single p-value unchanged", {
  p <- seq(from = 0.1, to = 0.9, by = 0.1)
  obs <- sapply(p, OmniP)
  expect_equal(obs, p)
})

test_that("OmniP repeated p-value unchanged", {
  p <- c(0.05, 0.05, 0.05)
  expect_equal(OmniP(p), 0.05)
  expect_equal(OmniP(p, w = c(1, 2, 1)), 0.05)
})

test_that("OmniP uniform weights equivalent", {
  p <- c(0.1, 0.2, 0.3)
  obs1 <- OmniP(p)
  obs2 <- OmniP(p, w = c(1, 1, 1))
  obs3 <- OmniP(p, w = c(2, 2, 2))
  expect_equal(obs1, obs2)
  expect_equal(obs1, obs3)
})

test_that("OmniP non-uniform weights change result", {
  p <- c(0.1, 0.05)
  obs <- OmniP(p, w = c(1, 2))
  expect_lt(obs, 0.1)
  obs <- OmniP(p, w = c(2, 1))
  expect_gt(obs, 0.05)
})

test_that("OmniP edge p=0 returns 0", {
  expect_equal(OmniP(c(0, 0.5)), 0)
  expect_equal(OmniP(0), 0)
})

test_that("OmniP edge p=1 returns 1", {
  expect_equal(OmniP(c(0.5, 1)), 1)
  expect_equal(OmniP(1), 1)
})

test_that("OmniP errors on invalid p or weight length", {
  expect_error(OmniP(c(0.1, 0.2), w = c(1, 2, 3)))
  expect_error(OmniP(c(-0.1, 0.5)))
  expect_error(OmniP(c(0.5, 1.1)))
  expect_error(OmniP(c(0, 1)))  # both 0 and 1
})
