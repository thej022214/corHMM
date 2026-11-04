test_that("collapsing independent characters keeps only active rate parameters", {
  phy <- ape::read.tree(text="(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);")
  dat <- data.frame(sp=letters[1:8], X=rep(c(0,0,1,1), 2), Y=rep(c("a","b","b","c"), 2))
  full <- getStateMat4Dat(dat, indep=TRUE, collapse=FALSE)$rate.mat
  collapsed <- getStateMat4Dat(dat, indep=TRUE, collapse=TRUE)$rate.mat
  subset <- full[rownames(collapsed), colnames(collapsed)]
  active <- sort(unique(subset[subset > 0]))
  expect_length(active, 6)
  expect_gt(max(active), length(active))
  want <- subset
  want[subset > 0] <- match(subset[subset > 0], active)
  expect_equal(collapsed, want)

  # Old saved or hand-specified matrices can have gaps. Their p slots must not shift.
  p <- seq_len(max(active)) / 10
  fit <- function(mat, rates) corHMM(phy, dat, 1, rate.mat=mat, p=rates,
      root.p=rep(.25,4), node.states="marginal", verbose=FALSE)
  compact <- fit(want, p[active])
  sparse <- fit(subset, p)
  expect_equal(sparse$solution, compact$solution)
  expect_equal(sparse$states, compact$states)
  expect_equal(sparse$loglik, compact$loglik)
  expect_equal(sparse$AIC, -2*compact$loglik + 2*6)
  expect_equal(sparse$AICc, sparse$AIC + 2*6*7/(8-6-1))
  expect_equal(attr(logLik(sparse),"df"),6L)
  expect_equal(AIC(sparse),sparse$AIC)
  expect_equal(AIC(compact),compact$AIC)
  expect_equal(sparse$index.mat[subset>0], subset[subset>0])

  # An actual optimizer call has no spare dimensions for generated matrices.
  est <- corHMM(phy, dat, 1, rate.mat=collapsed, ip=rep(.2,6),
      root.p=rep(.25,4), node.states="none", verbose=FALSE)
  expect_true(is.finite(est$loglik))
  expect_equal(est$AIC, -2*est$loglik + 12)
})

test_that("shared and forbidden transitions survive removal of empty slots", {
  dat <- data.frame(sp=letters[1:6], X=rep(0:1,each=3), Y=c("a","b","c","b","c","d"))
  m <- getStateMat4Dat(dat,indep=TRUE,collapse=TRUE)$rate.mat
  expect_equal(m["0|b","0|c"],m["1|b","1|c"])
  expect_gt(m["0|b","0|c"],0)
  expect_false(m["0|b","0|c"] == m["0|c","0|b"])
  expect_equal(m["0|a","1|d"],0)
  expect_equal(sort(unique(m[m>0])),seq_len(max(m)))
})
