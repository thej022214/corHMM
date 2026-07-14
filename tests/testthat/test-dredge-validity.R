library(testthat)
library(corHMM)

context("corHMMDredge model-space validity")

# rows are "from", columns are "to", diagonal is NA
ard <- matrix(c(NA, 2, 1, NA), 2, 2)
er <- matrix(c(NA, 1, 1, NA), 2, 2)
dollo <- matrix(c(NA, NA, 1, NA), 2, 2)  # 0 -> 1 only, state 1 absorbing
empty <- matrix(NA_real_, 2, 2)

test_that("is_valid_index_mat accepts models with at least one free rate", {
  expect_true(corHMM:::is_valid_index_mat(ard))
  expect_true(corHMM:::is_valid_index_mat(er))
  # irreversible models must stay in the search space
  expect_true(corHMM:::is_valid_index_mat(dollo))
  expect_true(corHMM:::is_valid_index_mat(t(dollo)))
})

test_that("is_valid_index_mat rejects an empty matrix", {
  expect_false(corHMM:::is_valid_index_mat(empty))
  # zeros are treated as absent rates, not as parameter 0
  expect_false(corHMM:::is_valid_index_mat(matrix(0, 2, 2)))
})

test_that("is_valid_index_mat rejects an isolated state", {
  # 3 states, state 3 has no rates in or out
  isolated <- matrix(NA_real_, 3, 3)
  isolated[1, 2] <- 1
  isolated[2, 1] <- 2
  expect_false(corHMM:::is_valid_index_mat(isolated))
  # connecting it in one direction only is enough
  isolated[2, 3] <- 3
  expect_true(corHMM:::is_valid_index_mat(isolated))
})

test_that("propose_stochastic_drop will not drop the last free rate", {
  # a one-parameter model: no drop is possible, so no matrix is proposed
  one_par <- list(
    index.mat = dollo,
    solution = matrix(c(NA, NA, 0.5, NA), 2, 2)
  )
  expect_null(corHMM:::propose_stochastic_drop(one_par, drop.threshold = 1e-7))
})

test_that("propose_stochastic_drop never proposes an invalid matrix", {
  # both rates sit below drop.threshold, so the proposer may try to drop both
  tiny <- list(
    index.mat = ard,
    solution = matrix(c(NA, 1e-10, 1e-10, NA), 2, 2)
  )
  set.seed(1)
  for (i in 1:100) {
    proposal <- corHMM:::propose_stochastic_drop(tiny, drop.threshold = 1e-7)
    if (!is.null(proposal)) {
      expect_true(corHMM:::is_valid_index_mat(proposal))
    }
  }
})
