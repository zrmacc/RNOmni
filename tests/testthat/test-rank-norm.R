test_that("Test rank normalization.", {
  
  # Case 1.
  u <- 1
  obs <- RankNorm(u)
  expect_equal(obs, 0)
  
  # Case 2.
  u <- c(1, 2)
  obs <- RankNorm(u, k = 0)
  exp <- c(stats::qnorm(1/3), stats::qnorm(2/3))
  expect_equal(obs, exp)
  
})
