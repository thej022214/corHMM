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

test_that("dropping the last free rate is rejected by name", {
  one_par <- list(
    rate.cat = 1,
    index.mat = dollo,
    solution = matrix(c(NA, NA, 0.5, NA), 2, 2)
  )
  proposal <- corHMM:::propose_stochastic_drop(one_par, drop.threshold = 1e-7)$new_index_mat
  expect_equal(corHMM:::validate_index_mat(proposal), "no free parameters")
})

test_that("every drop proposal is a model or names why it is not", {
  # both rates sit below drop.threshold, so the proposer may try to drop both
  tiny <- list(
    rate.cat = 1,
    index.mat = ard,
    solution = matrix(c(NA, 1e-10, 1e-10, NA), 2, 2)
  )
  set.seed(1)
  for (i in 1:100) {
    proposal <- corHMM:::propose_stochastic_drop(tiny, drop.threshold = 1e-7)$new_index_mat
    reason <- corHMM:::validate_index_mat(proposal)
    expect_true(is.null(reason) ||
      reason %in% c("no proposal", "no free parameters", "isolated state"))
  }
})

test_that("validate_index_mat names each way a matrix fails", {
  expect_null(corHMM:::validate_index_mat(ard))
  expect_equal(corHMM:::validate_index_mat(NULL), "no proposal")
  expect_equal(corHMM:::validate_index_mat(empty), "no free parameters")
  isolated <- matrix(NA_real_, 3, 3)
  isolated[1, 2] <- 1
  isolated[2, 1] <- 2
  expect_equal(corHMM:::validate_index_mat(isolated), "isolated state")
  expect_equal(corHMM:::validate_index_mat(dollo, require.entry = TRUE), "unenterable state")
  hmm <- getFullMat(replicate(2, ard, simplify = FALSE), getStateMat(2))
  hmm[hmm == 0] <- NA
  expect_null(corHMM:::validate_index_mat(hmm, 2))
  # a rate class with no rates at all is not a model, whatever the dimnames say
  dimnames(hmm) <- NULL
  expect_null(corHMM:::validate_index_mat(hmm, 2))
  hmm[3:4, ] <- hmm[, 3:4] <- NA
  expect_equal(corHMM:::validate_index_mat(hmm, 2), "isolated state")
})

test_that("diagonal parameters are never valid", {
  bad <- ard
  bad[1, 1] <- 3
  expect_equal(corHMM:::validate_index_mat(bad), "diagonal parameter")
  expect_false(corHMM:::is_valid_index_mat(bad))
})

test_that("models have one source component", {
  two_sources <- matrix(NA_real_, 4, 4)
  two_sources[1, 2] <- 1
  two_sources[3, 1] <- 2
  two_sources[4, 2] <- 2
  expect_equal(corHMM:::validate_index_mat(two_sources, 2),
    "multiple source components")
  precursor <- matrix(NA_real_, 3, 3)
  precursor[1, 2] <- 1
  precursor[2, 3] <- 2
  expect_null(corHMM:::validate_index_mat(precursor))
})

test_that("transitions outside the legal mask are never valid", {
  allowed <- matrix(NA_real_, 4, 4)
  allowed[1, 2] <- allowed[2, 1] <- 1
  allowed[2, 4] <- allowed[4, 2] <- 1
  allowed[4, 3] <- allowed[3, 4] <- 1
  candidate <- allowed
  candidate[1, 4] <- 2
  expect_equal(corHMM:::validate_index_mat(candidate, allowed_mat = allowed),
    "transition outside model space")
})

test_that("a user mask can permit a simultaneous transition", {
  allowed <- getFullMat(replicate(2, ard, simplify = FALSE), getStateMat(2))
  allowed[allowed == 0] <- NA
  candidate <- allowed
  candidate[3, 2] <- max(allowed, na.rm = TRUE) + 1
  expect_equal(corHMM:::validate_index_mat(candidate, 2,
    allowed_mat = allowed), "transition outside model space")
  allowed[3, 2] <- candidate[3, 2]
  expect_null(corHMM:::validate_index_mat(candidate, 2,
    allowed_mat = allowed))
  fit <- list(index.mat = candidate)
  fit$index.mat[3, 2] <- NA
  proposal <- corHMM:::propose_stochastic_free(fit, allowed)
  expect_false(is.na(proposal[3, 2]))
})

test_that("free moves preserve the legal transition mask", {
  allowed <- matrix(c(NA, 2, 1, NA), 2, 2)
  fit <- list(index.mat = matrix(c(NA, NA, 1, NA), 2, 2))
  for(i in 1:100) {
    proposal <- corHMM:::propose_stochastic_free(fit, allowed)
    expect_true(is.null(proposal) || all(is.na(proposal[is.na(allowed)])))
    if(!is.null(proposal)) expect_true(all(is.na(diag(proposal))))
  }
})

test_that("a drop that strands a hidden class collapses it", {
  hmm <- getFullMat(replicate(2, ard, simplify = FALSE), getStateMat(2))
  hmm[hmm == 0] <- NA
  stranded <- corHMM:::dropStateMatPars(hmm, c(3, 4, 5, 6))
  pruned <- corHMM:::prune_isolated_states(stranded, 2)
  expect_equal(pruned$rate_cat, 1)
  expect_equal(dim(pruned$new_index_mat), c(2, 2))
  # stranding a single state collapses nothing and is rejected by name
  one_state <- hmm
  one_state[4, ] <- one_state[, 4] <- NA
  pruned <- corHMM:::prune_isolated_states(one_state, 2)
  expect_null(pruned$rate_cat)
  expect_equal(corHMM:::validate_index_mat(pruned$new_index_mat, 2), "isolated state")
})

test_that("rates at the floor are not merge candidates", {
  expect_equal(sort(corHMM:::stochastic_merge_pars(c(1e-10, 1e-10, 0.5, 0.52), 0)), c(3, 4))
  expect_null(corHMM:::stochastic_merge_pars(c(1e-10, 1e-10, 0.5), 0))
})

test_that("printing preserves free rates at the optimizer lower bound", {
  data(primates)
  index_mat <- matrix(NA_real_, 4, 4)
  index_mat[1, 2] <- 1
  solution <- matrix(NA_real_, 4, 4)
  solution[1, 2] <- 2e-12
  fit <- list(loglik = -1, AIC = 4, AICc = 5, rate.cat = 1,
    phy = primates[[1]], data = primates[[2]], index.mat = index_mat,
    solution = solution, lower.bound = 1e-10)
  class(fit) <- "corhmm"
  dredge <- list(fit)
  class(dredge) <- "corhmm.dredge"
  attr(dredge, "dredge_history") <- list(list(rate_category = 1,
    iterations = 1, acceptance_rate = 0, restart_count = 0,
    best_fit = fit, best_score = fit$AIC, unique_structures = 1,
    stop_reason = "stalled", rejections = character()))

  printed <- capture.output(print(dredge))
  expect_true(any(grepl("2e-12", printed, fixed = TRUE)))
  expect_false(any(grepl("shown as NA", printed, fixed = TRUE)))
})

test_that("dredge restart fits report rates on the natural scale", {
  data(primates)
  phy <- ape::multi2di(primates[[1]])
  phy$edge.length <- phy$edge.length + 1e-6
  dat <- primates[[2]][, 1:2]
  index_mat <- matrix(c(NA, 1, NA, NA), 2, 2, byrow = TRUE)
  set.seed(101)
  fit <- corHMM:::corHMMDredgeBase(phy, dat, 1, root.p = "maddfitz",
    pen.type = "l1", lambda = 0, rate.mat = index_mat,
    node.states = "none", nstarts = 1, use_RTMB = TRUE)
  reproduced <- corHMM:::corHMMDredgeBase(phy, dat, 1,
    root.p = "maddfitz", pen.type = "l1", lambda = 0,
    rate.mat = fit$index.mat, node.states = "none",
    p = unname(MatrixToPars(fit)), use_RTMB = TRUE)
  expect_equal(reproduced$loglik, fit$loglik, tolerance = 1e-6)
  expect_equal(reproduced$AIC, fit$AIC, tolerance = 1e-6)
})

test_that("annealing temperature spans the iteration budget", {
  temps <- vapply(1:1000, corHMM:::annealing_temperature, numeric(1),
    max.iterations = 1000, initial.temp = 2, cooling.rate = 0.95,
    temp.schedule = "exponential", epoch.start = 1)
  expect_equal(temps[1], 2)
  expect_equal(temps[1000], 0.001, tolerance = 1e-12)
  expect_true(all(diff(temps) < 0))
})

test_that("BIC uses tip count and parameter count", {
  scores <- corHMM:::information_criteria(-100, 3, 60)
  expect_equal(scores$AIC, 206)
  expect_equal(scores$BIC, 200 + 3 * log(60))
  expect_equal(scores$AICc, 200 + 6 * 60 / 56)
})

test_that("an unproductive restart interval triggers polishing", {
  expect_false(corHMM:::polish_due(19, 19, 500, 20))
  expect_true(corHMM:::polish_due(20, 20, 500, 20))
  expect_true(corHMM:::polish_due(500, 0, 500, 20))
})

test_that("polishing follows improving neighbors to a local optimum", {
  initial <- list(AIC = 10, state = 0)
  neighbors <- function(fit) {
    if (fit$state == 0) return(list(list(AIC = 9, state = 1), list(AIC = 11, state = 3)))
    if (fit$state == 1) return(list(list(AIC = 8, state = 2), list(AIC = 9.5, state = 4)))
    list()
  }
  result <- corHMM:::polish_neighborhood(initial, neighbors, identity, "AIC")
  expect_equal(result$best_fit$AIC, 8)
  expect_equal(length(result$path), 2)
})

test_that("polishing bounds uncached evaluations", {
  initial <- list(AIC = 10)
  proposals <- lapply(1:10, function(i) list(AIC = 10 + i))
  evaluated <- 0L
  result <- corHMM:::polish_neighborhood(initial, function(fit) proposals,
    function(candidate) {
      evaluated <<- evaluated + 1L
      candidate
    }, "AIC", max.new.evaluations = 3,
    is.cached = function(candidate) FALSE)
  expect_equal(evaluated, 3)
  expect_equal(result$new_evaluations, 3)
  expect_equal(result$deferred, 7)
  expect_false(result$complete)
})

test_that("deterministic merge neighbors exclude floor parameters", {
  index_mat <- matrix(c(NA, 2, 1, NA), 2, 2)
  solution <- matrix(c(NA, 0.5, 1e-10, NA), 2, 2)
  dimnames(index_mat) <- dimnames(solution) <- list(c("0", "1"), c("0", "1"))
  fit <- list(rate.cat = 1, index.mat = index_mat, solution = solution)
  neighbors <- corHMM:::deterministic_neighbors(fit, index_mat, 1)
  expect_false(any(vapply(neighbors, function(x) x$move_type == "merge", logical(1))))
})

test_that("polishing can relocate a floor parameter", {
  max_mat <- matrix(NA_real_, 3, 3)
  max_mat[1, 2] <- max_mat[2, 1] <- max_mat[2, 3] <- max_mat[3, 2] <- 1
  index_mat <- matrix(NA_real_, 3, 3)
  index_mat[1, 2] <- 1
  index_mat[2, 3] <- 2
  solution <- matrix(NA_real_, 3, 3)
  solution[1, 2] <- 1e-10
  solution[2, 3] <- 0.5
  dimnames(index_mat) <- dimnames(solution) <- list(as.character(1:3), as.character(1:3))
  fit <- list(rate.cat = 1, index.mat = index_mat, solution = solution)
  neighbors <- corHMM:::deterministic_neighbors(fit, max_mat, 1)
  relocated <- neighbors[vapply(neighbors, function(x)
    x$move_type == "polish_relocate", logical(1))]
  expect_true(any(vapply(relocated, function(x)
    is.na(x$new_index_mat[1, 2]) && !is.na(x$new_index_mat[2, 1]), logical(1))))
})
