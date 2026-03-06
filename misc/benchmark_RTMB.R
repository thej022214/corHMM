remotes::install_github("bbolker/corHMM", ref = "bolker_clean")
## source("~/R/pkgs/corHMM/R/RTMB_devcorhmm.R")
library(corHMM)
source("misc/simfun.R")
set.seed(101)
fitfun  <- function(dat, ..., rate.cat = 1) {
  tt <- system.time(
    invisible(capture.output(x <- with(dat, corHMM(phy, data, rate.cat = rate.cat, ...))))
  )
  attr(x, "time") <- tt
  x
}
parfun <- function(x) log(na.omit(c(x$solution)))

sumfun <- function(ntrait = 2, ntaxa = 200, model = "ARD", seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  seed <- seed %||% NA
  ss <- simfun(ntrait = ntrait, ntaxa = ntaxa)
  fit_orig <- fitfun(ss, model = model, ...)
  fit_RTMB <- fitfun(ss, use_RTMB = TRUE, model = model, ...)
  p_orig <- parfun(fit_orig)
  p_RTMB <- parfun(fit_RTMB)
  p_diff <- p_orig-p_RTMB
  ## truncate par values at -10 (these diffs are mostly irrelevant)
  p_orig_trunc <- pmax(p_orig, -10)
  p_RTMB_trunc <- pmax(p_RTMB, -10)
  p_diff_trunc <- p_orig_trunc - p_RTMB_trunc
  ## get RMSE vs true rates on truncated scale ...
  p_RTMB_rmse <- sqrt(mean((p_RTMB_trunc - ss$true_rates)^2))
  p_orig_rmse <- sqrt(mean((p_orig_trunc - ss$true_rates)^2))
  data.frame(seed, ntrait, ntaxa, model,
             RTMB_opt.time = fit_RTMB$opt.time[["elapsed"]],
             orig_opt.time = fit_orig$opt.time[["elapsed"]],
             RTMB_tot.time = attr(fit_RTMB, "time")[["elapsed"]],
             orig_tot.time = attr(fit_orig, "time")[["elapsed"]],
             RTMB_loglik = fit_RTMB$loglik,
             orig_loglik = fit_orig$loglik,
             par.rmsdiff = sqrt(mean(p_diff^2)),
             par.maxdiff = max(abs(p_diff)),
             par.rmsdiff.trunc = sqrt(mean(p_diff_trunc^2)),
             par.maxdiff.trunc = max(abs(p_diff_trunc)),
             RTMB_rmse = p_RTMB_rmse,
             orig_rmse = p_orig_rmse
             )
}

ntaxvec <- round(exp(seq(log(20), log(1000), length.out = 15)))
dd <- expand.grid(seed = 101:105, ntaxa = ntaxvec, ntrait = 2:3, model = c("ARD", "SYM")) |>
  transform(model = as.character(model)) ## factor messes up do.call, converts to numeric

## test
sumfun(seed = 105, ntrait = 3, ntaxa = 20)

res <- list()
for (i in 1:nrow(dd)) {
  print(dd[i,])
  res[[i]] <- do.call(sumfun, dd[i,])
  saveRDS(res, file = "benchmark2.rds")
}
res <- do.call(rbind, res)
saveRDS(res, file = "benchmark2.rds")
