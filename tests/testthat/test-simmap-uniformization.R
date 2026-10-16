test_that("uniformization bridges have the analytic symmetric event distribution", {
  Q <- matrix(c(-0.8, 0.8, 0.8, -0.8), 2, 2)
  time <- 1.1
  lambda <- 0.8 * time
  uni <- corHMM:::makeUniformizationTable(Q, time, min(expm::expm(Q * time)))
  set.seed(210)
  same <- replicate(4000, length(corHMM:::sampleBranchBridge(1, 1, time, uni)) - 1L)
  different <- replicate(4000, length(corHMM:::sampleBranchBridge(1, 2, time, uni)) - 1L)
  expect_true(all(same %% 2 == 0))
  expect_true(all(different %% 2 == 1))
  expect_equal(mean(same), lambda * tanh(lambda), tolerance = 0.05)
  expect_equal(mean(different), lambda / tanh(lambda), tolerance = 0.05)
  expect_equal(mean(same == 0), 1 / cosh(lambda), tolerance = 0.025)
})

test_that("irreversible bridge times follow the truncated exponential distribution", {
  rate <- 1.3
  time <- 0.7
  Q <- matrix(c(-rate, 0, rate, 0), 2, 2)
  uni <- corHMM:::makeUniformizationTable(Q, time, 1 - exp(-rate * time))
  set.seed(211)
  paths <- replicate(4000, corHMM:::sampleBranchBridge(1, 2, time, uni), simplify = FALSE)
  expect_true(all(vapply(paths, function(path) identical(names(path), c("1", "2")), logical(1))))
  shift.times <- vapply(paths, function(path) unname(path[1]), numeric(1))
  expect_equal(mean(shift.times), 1 / rate - time / expm1(rate * time), tolerance = 0.015)
  at <- time / 3
  expect_equal(mean(shift.times <= at),
               (1 - exp(-rate * at)) / (1 - exp(-rate * time)), tolerance = 0.025)
})

test_that("uniformization preserves endpoints, allowed transitions, and branch lengths", {
  Q <- matrix(c(-2, 0, 0, 2, -1, 0, 0, 1, 0), 3, 3)
  time <- 1e-8
  uni <- corHMM:::makeUniformizationTable(Q, time, 1e-17)
  set.seed(212)
  for (i in 1:30) {
    path <- corHMM:::sampleBranchBridge(1, 3, time, uni)
    states <- as.integer(names(path))
    expect_equal(states, 1:3)
    expect_true(all(path >= 0))
    expect_equal(sum(path), time, tolerance = 1e-20)
  }
  expect_equal(corHMM:::sampleBranchBridge(3, 3, time, uni), c("3" = time))
  expect_error(corHMM:::sampleBranchBridge(3, 1, time, uni), "zero probability")
  expect_error(corHMM:::sampleBranchBridge(1, 2, 0, uni), "zero probability")
  expect_equal(corHMM:::sampleBranchBridge(1, 1, 0, uni), c("1" = 0))
  zero <- corHMM:::makeUniformizationTable(matrix(0, 2, 2), 2, 1)
  expect_equal(corHMM:::sampleBranchBridge(1, 1, 2, zero), c("1" = 2))
  expect_error(corHMM:::sampleBranchBridge(1, 2, 2, zero), "zero probability")
})

test_that("uniformization reports resource limits instead of fabricating a path", {
  Q <- matrix(c(-1, 1, 1, -1), 2, 2)
  expect_error(corHMM:::makeUniformizationTable(Q, 1e6, 0.5), "too large.*rejection")
})

test_that("uniformization maps support fitted and matrix models with hidden classes", {
  phy <- ape::read.tree(text = "((a:0.2,b:0.3):0.4,(c:0.5,d:0.6):0.7);")
  data <- data.frame(sp = letters[1:4], A = c(0, 0, 1, 1), B = c(0, 1, 0, 1))
  for (rc in 1:2) {
    index <- getStateMat4Dat(data)$rate.mat
    if (rc > 1) index <- getFullMat(rep(list(index), rc))
    fit <- corHMM(phy, data, rate.cat = rc, rate.mat = index,
                  p = seq(0.1, 1, length.out = max(index, na.rm = TRUE)),
                  node.states = "none", verbose = FALSE)
    set.seed(213)
    from.fit <- makeSimmap(model = fit, nSim = 3, method = "uniformization")
    set.seed(213)
    from.matrix <- makeSimmap(fit$phy, fit$data, fit$solution, rc,
                             nSim = 3, method = "uniformization")
    expect_identical(from.fit, from.matrix)
    for (map in from.fit) {
      expect_equal(unname(rowSums(map$mapped.edge)), map$edge.length)
      endpoints <- vapply(seq_along(map$tip.label), function(i) {
        tail(names(map$maps[[which(map$edge[, 2] == i)]]), 1)
      }, character(1))
      indices <- (match(endpoints, colnames(fit$solution)) - 1L) %% 4L + 1L
      expect_equal(indices, c(1, 3, 2, 4))
      for (edge in seq_along(map$maps)) {
        path <- map$maps[[edge]]
        states <- match(names(path), colnames(map$Q))
        if (length(states) > 1) {
          expect_true(all(diff(states) != 0))
          expect_true(all(map$Q[cbind(head(states, -1), tail(states, -1))] > 0))
        }
        durations <- vapply(colnames(map$Q), function(state) sum(path[names(path) == state]), numeric(1))
        expect_equal(unname(map$mapped.edge[edge, ]), unname(durations))
      }
    }
    set.seed(214)
    original <- makeSimmap(model = fit, nSim = 2)
    set.seed(214)
    explicit <- makeSimmap(model = fit, nSim = 2, method = "rejection")
    expect_identical(original, explicit)
    expect_error(makeSimmap(model = fit, method = "invalid"), "arg")
  }
})
