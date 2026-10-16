simmap_fit_fixture <- function(rate.cat = 1L, states = c("0_0", "0_1", "1_0", "1_1"),
                               codes = as.character(seq_along(states))) {
  n <- length(states)
  k <- n * rate.cat
  Q <- matrix(0, k, k)
  traits <- do.call(rbind, strsplit(states, "_", fixed = TRUE))
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      observed.i <- (i - 1L) %% n + 1L
      observed.j <- (j - 1L) %% n + 1L
      same.class <- (i - 1L) %/% n == (j - 1L) %/% n
      if ((same.class && sum(traits[observed.i, ] != traits[observed.j, ]) == 1L) ||
          (!same.class && observed.i == observed.j)) {
        Q[i, j] <- 0.1 + (i + 2 * j) / (10 * k)
      }
    }
  }
  Q[Q == 0] <- NA_real_
  labels <- paste0("(", rep(seq_len(n), rate.cat), ",R",
                   rep(seq_len(rate.cat), each = n), ")")
  dimnames(Q) <- list(labels, labels)
  phy <- ape::stree(length(codes), type = "left")
  phy$edge.length <- seq(0.2, 0.8, length.out = nrow(phy$edge))
  structure(list(solution = Q, rate.cat = rate.cat, phy = phy,
                 data.legend = data.frame(sp = phy$tip.label, d = codes),
                 root.p = "yang"), class = "corhmm")
}

simmap_tip_indices <- function(map, fit) {
  endpoints <- vapply(seq_along(map$tip.label), function(i) {
    tail(names(map$maps[[which(map$edge[, 2] == i)]]), 1L)
  }, character(1))
  (match(endpoints, colnames(fit$solution)) - 1L) %%
    (nrow(fit$solution) / fit$rate.cat) + 1L
}

test_that("mapping a 2.8 correlated fit preserves its fitted state identities", {
  for (rc in 1:2) {
    fit <- simmap_fit_fixture(rc)
    set.seed(81)
    maps <- makeSimmap(model = fit, nSim = 3)
    for (map in maps) {
      expect_equal(simmap_tip_indices(map, fit), 1:4)
      Q <- fit$solution
      Q[is.na(Q)] <- 0
      diag(Q) <- -rowSums(Q)
      expect_equal(map$Q, Q)
      durations <- t(vapply(map$maps, function(edge) {
        vapply(colnames(map$Q), function(state) sum(edge[names(edge) == state]), numeric(1))
      }, numeric(ncol(map$Q))))
      dimnames(durations) <- dimnames(map$mapped.edge)
      expect_equal(map$mapped.edge, durations)
    }
  }
})

test_that("saved numeric root priors and fixed states retain their matrix positions", {
  fit <- simmap_fit_fixture(2L)
  fit$root.p <- c(0, 1, rep(0, 6))
  root_labels <- function(map) {
    root <- length(map$tip.label) + 1L
    vapply(which(map$edge[, 1] == root), function(i) names(map$maps[[i]])[1], character(1))
  }
  set.seed(82)
  map <- makeSimmap(model = fit)[[1]]
  expect_true(all(root_labels(map) == colnames(fit$solution)[2]))

  prior <- c(0, 0, 1, rep(0, 5))
  set.seed(82)
  map <- makeSimmap(model = fit, root.p = prior)[[1]]
  expect_true(all(root_labels(map) == colnames(fit$solution)[3]))

  set.seed(82)
  map <- makeSimmap(model = fit, root.p = "yang",
                    fix.node = length(fit$phy$tip.label) + 1L, fix.state = 6L)[[1]]
  expect_true(all(root_labels(map) == colnames(fit$solution)[6]))
})

test_that("fit mapping supports collapsed and unobserved combinations", {
  for (rc in 1:2) {
    collapsed <- simmap_fit_fixture(rc, states = c("0_0", "0_1", "1_1"))
    uncollapsed <- simmap_fit_fixture(rc, codes = c("1", "2", "4", "4"))
    set.seed(83)
    expect_equal(simmap_tip_indices(makeSimmap(model = collapsed)[[1]], collapsed), 1:3)
    set.seed(83)
    expect_equal(simmap_tip_indices(makeSimmap(model = uncollapsed)[[1]], uncollapsed),
                 c(1, 2, 4, 4))
  }
})

test_that("fit mapping preserves numeric indices above nine and ambiguous tips", {
  states <- sort(as.vector(outer(0:2, 0:3, paste, sep = "_")))
  fit <- simmap_fit_fixture(states = states)
  set.seed(84)
  expect_equal(simmap_tip_indices(makeSimmap(model = fit)[[1]], fit), 1:12)

  fit <- simmap_fit_fixture(codes = c("1", "2&3", "?", "4"))
  set.seed(85)
  indices <- simmap_tip_indices(makeSimmap(model = fit)[[1]], fit)
  expect_equal(indices[c(1, 4)], c(1, 4))
  expect_true(indices[2] %in% c(2, 3))
  expect_true(indices[3] %in% 1:4)
})

test_that("current fitted objects agree with the existing matrix interface", {
  phy <- simmap_fit_fixture()$phy
  data <- data.frame(sp = phy$tip.label,
                     A = factor(c(0, 0, 1, 1), levels = c(1, 0)),
                     B = c(0, 1, 0, 1))
  fit <- corHMM(phy, data, rate.cat = 1, p = seq(0.1, 0.8, by = 0.1),
                node.states = "none", verbose = FALSE)
  set.seed(86)
  from.fit <- makeSimmap(model = fit, nSim = 2)
  set.seed(86)
  from.matrix <- makeSimmap(tree = fit$phy, data = fit$data, model = fit$solution,
                           rate.cat = fit$rate.cat, root.p = fit$root.p,
                           collapse = fit$collapse, nSim = 2)
  expect_identical(from.fit, from.matrix)
})

test_that("fit mapping validates saved coding and conflicting arguments", {
  fit <- simmap_fit_fixture()
  expect_error(makeSimmap(model = fit, data = fit$data.legend), "omit.*data")
  expect_error(makeSimmap(model = fit, rate.cat = 2), "rate.cat.*fit")
  fit$data.legend$d[1] <- "5"
  expect_error(makeSimmap(model = fit), "data.legend.*state")
  fit$data.legend <- NULL
  expect_error(makeSimmap(model = fit), "data.legend")
  fit <- simmap_fit_fixture()
  fit$rate.cat <- 3L
  expect_error(makeSimmap(model = fit), "solution matrix.*rate.cat")
})
