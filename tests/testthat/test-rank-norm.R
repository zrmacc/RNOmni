test_that("RankNorm single observation gives 0 (Blom default)", {
  u <- 1
  obs <- RankNorm(u)
  expect_equal(obs, 0)
})

test_that("RankNorm with k=0 matches (r - 0)/(n - 0 + 1)", {
  u <- c(1, 2)
  obs <- RankNorm(u, k = 0)
  exp <- c(stats::qnorm(1/3), stats::qnorm(2/3))
  expect_equal(obs, exp)
})

test_that("RankNorm default k=0.375 is Blom", {
  u <- c(1, 2, 3)
  obs <- RankNorm(u, k = 0.375)
  n <- 3
  r <- rank(u, ties.method = "average")
  exp <- stats::qnorm((r - 0.375) / (n - 2 * 0.375 + 1))
  expect_equal(obs, exp)
})

test_that("RankNorm handles ties with ties.method", {
  u <- c(1, 1, 2, 2)
  obs_avg <- RankNorm(u, ties.method = "average")
  obs_first <- RankNorm(u, ties.method = "first")
  expect_length(obs_avg, 4)
  expect_length(obs_first, 4)
  expect_equal(obs_avg[1], obs_avg[2])
  expect_equal(obs_avg[3], obs_avg[4])
})

test_that("RankNorm rejects invalid k", {
  u <- c(1, 2, 3)
  expect_error(RankNorm(u, k = -0.1))
  expect_error(RankNorm(u, k = 0.6))
})

test_that("RankNorm rejects NA in input", {
  u <- c(1, 2, NA, 4)
  expect_error(RankNorm(u))
})

test_that("RankNorm rejects non-vector", {
  u <- matrix(1:4, ncol = 1)
  expect_error(RankNorm(u))
})
