## basic (quick?) corHMM fitting test with simulated data
## (via ape::rtree and phangorn::simSeq, processed)
library(corHMM)
library(testthat)
set.seed(7)
m3 <- ape::rtree(20)
g1 <- reorder(m3, "pruningwise")

Q <- matrix(
  c(0, 0.017, 0.124, 0, 0.234,
    0, 0, 0.26, 0.052,
    0, 0, 0.108, 0,
    0.136, 0.107, 0),
  nrow = 4L,
  ncol = 4L,
  dimnames = list(c("0|0", "1|0", "0|1", "1|1"), c("0|0", "1|0", "0|1", "1|1"))
)

traitMatrix <- data.frame(
  nm = c(
    "t15", "t8", "t12", "t19", "t2", "t20", "t6", "t14", "t11", "t5", "t4", "t13",
    "t16", "t7", "t17", "t10", "t9", "t3", "t18", "t1"
  ),
  V1 = rep(c(0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L), c(2L, 6L, 1L, 1L, 6L, 1L, 1L, 1L, 1L)),
  V2 = rep(c(1L, 0L, 1L), c(3L, 4L, 13L)),
  row.names = c(
    "t15", "t8", "t12", "t19", "t2", "t20", "t6", "t14", "t11", "t5", "t4", "t13",
    "t16", "t7", "t17", "t10", "t9", "t3", "t18", "t1"
  )
)

## snapshot: dput(fitted_pars) from below
best_pars <- c(-1.96836334206339, 1.90081570465979,
               -2.28245419352619, 1.12531062633252)
best_loglik <- -15.0219701261726
## suppress output
cfun  <- function(...) {
  invisible(capture.output(x <- corHMM(...)))
  x
}
get_pars <- function(x) c(log(na.omit(c(x$solution))))
                          
## log-likelihood calc for fixed parameters
dev_orig <- cfun(g1, traitMatrix, rate.cat = 1, p = exp(best_pars))
dev_RTMB <- cfun(g1, traitMatrix, rate.cat = 1, p = exp(best_pars), use_RTMB = TRUE)
expect_equal(dev_orig$loglik, dev_RTMB$loglik)

## reset seed for interactive messing-about convenience
set.seed(101)
t_orig <-
  system.time(fit_orig <- cfun(g1, traitMatrix, rate.cat = 1)
              )

fitted_pars_orig <- get_pars(fit_orig)
expect_equal(fitted_pars_orig, best_pars)
expect_equal(fit_orig$loglik, best_loglik)

t_RTMB <-
  system.time(
    fit_RTMB <- cfun(g1, traitMatrix, rate.cat = 1, use_RTMB = TRUE)
  )
fitted_pars_RTMB <- c(log(na.omit(c(fit_RTMB$solution))))
expect_equal(fitted_pars_RTMB, best_pars, tolerance = 5e-3)
expect_equal(fit_RTMB$loglik, best_loglik)

cbind(t_RTMB, t_orig, topt_RTMB = fit_RTMB$opt.time, topt_orig = fit_orig$opt.time)

## devtools::load_all()
## source("tests/testthat/test-compRTMB.R", echo = TRUE)
