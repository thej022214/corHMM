coef.corhmm <- function(object, log = TRUE, ...) {
  p <- object$args.list$p
  if (log) p else exp(p)
}

logLik.corhmm <- function(object, ...) {
  L <- object$loglik
  attr(L, "df") <- length(coef(object))
  ## nobs doesn't necessarily make sense - could use number of tips but ...
  class(L) <- "logLik"
  L
}

## SEQUENTIAL anova
anova.corhmm <- function(object, ...) {
  lklist <- lapply(c(list(object),list(...)), logLik)
  dvec <- -2*vapply(lklist, c, FUN.VALUE = numeric(1))
  nvec <- vapply(lklist, attr, "df", FUN.VALUE = integer(1))
  dd <- data.frame(deviance = dvec,
                   df = nvec)
  dd <- within(dd, {
    delta_dev <- c(NA, abs(diff(deviance)))
    delta_df <- c(NA, abs(diff(df)))
    pval <- pchisq(delta_dev, df = delta_df, lower.tail = FALSE)
  })
  return(dd[,c("deviance", "df", "delta_dev", "delta_df", "pval")])
}
