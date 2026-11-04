#' @export
coef.corhmm <- function(object, log = TRUE, ...) {
  p <- object$args.list$p
  if (log) p else exp(p)
}

nCorHMMParameters <- function(indices, fog) {
  length(unique(indices[!is.na(indices) & indices > 0])) + length(unique(fog))
}

#' @export
logLik.corhmm <- function(object, ...) {
  L <- object$loglik
  attr(L, "df") <- if (is.null(object$index.mat)) length(coef(object)) else
    nCorHMMParameters(object$index.mat, object$fog.vec)
  ## nobs doesn't necessarily make sense - could use number of tips but ...
  class(L) <- "logLik"
  L
}

## SEQUENTIAL anova
#' @export
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
