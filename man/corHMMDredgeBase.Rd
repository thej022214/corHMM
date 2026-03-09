\name{corHMMDredgeBase}
\alias{corHMMDredgeBase}
\title{Fit a Single Penalized Hidden Markov Model}
\description{
\code{corHMMDredgeBase} fits a single hidden Markov model (HMM) to a
phylogenetic tree and discrete character data with optional regularization.
It is the workhorse called repeatedly by \code{\link{corHMMDredge}} during
the simulated annealing search, but can also be called directly to fit a
specific model structure.
}
\usage{
corHMMDredgeBase(phy, data, rate.cat, root.p="maddfitz", tip.fog=NULL,
  fog.ip=0.01, pen.type="l1", lambda=1, rate.mat=NULL,
  node.states="marginal", fixed.nodes=FALSE, ip=NULL, nstarts=0,
  n.cores=1, get.tip.states=FALSE, lewis.asc.bias=FALSE, collapse=FALSE,
  lower.bound=1e-10, upper.bound=100, opts=NULL, p=NULL,
  use_RTMB=FALSE, prep=NULL)
}
\arguments{
  \item{phy}{A phylogenetic tree of class \code{"phylo"}.}
  \item{data}{A data frame with species names in the first column and
    discrete character states in subsequent columns.}
  \item{rate.cat}{An integer specifying the number of hidden rate
    categories.}
  \item{root.p}{Root state prior. Can be \code{"maddfitz"}, \code{"yang"},
    a numeric vector, or \code{NULL}. See \code{\link{corHMM}}.}
  \item{tip.fog}{Optional tip observation error. A fixed numeric value
    less than 1 is treated as a per-state misclassification probability;
    if \code{sum(tip.fog) >= 1} the probabilities are estimated.}
  \item{fog.ip}{Starting value for tip fog estimation. Defaults to 0.01.}
  \item{pen.type}{Regularization type: \code{"l1"} (lasso),
    \code{"l2"} (ridge), \code{"er"} (equal-rates penalty), or
    \code{"unreg"} (none). Defaults to \code{"l1"}.}
  \item{lambda}{Penalty strength from 0 (no penalty) to 1 (full
    penalization). Defaults to 1.}
  \item{rate.mat}{An optional index matrix specifying the parameter
    structure of the rate matrix. If \code{NULL} an all-rates-different
    model is used.}
  \item{node.states}{Ancestral state reconstruction method:
    \code{"marginal"}, \code{"joint"}, \code{"scaled"}, or
    \code{"none"}.}
  \item{fixed.nodes}{Logical. Whether internal node states are fixed via
    node labels in \code{phy}.}
  \item{ip}{Optional starting values for optimisation.}
  \item{nstarts}{Number of random restarts. Defaults to 0.}
  \item{n.cores}{Number of cores for parallel random restarts.}
  \item{get.tip.states}{Logical. Whether to return tip state
    reconstructions.}
  \item{lewis.asc.bias}{Logical. Whether to correct for ascertainment
    bias following Lewis (2001).}
  \item{collapse}{Logical. Whether to collapse unobserved state
    combinations.}
  \item{lower.bound, upper.bound}{Bounds for rate parameter estimates.}
  \item{opts}{A list of options passed to the \code{nloptr} optimiser.}
  \item{p}{Optional fixed parameter vector. If supplied the likelihood
    is evaluated at these values without optimisation.}
  \item{use_RTMB}{Logical. Whether to use the RTMB automatic
    differentiation backend for gradient-based optimisation. Defaults to
    \code{FALSE} when called directly; \code{\link{corHMMDredge}} sets
    this to \code{TRUE} by default.}
  \item{prep}{An optional precomputed prep object from
    \code{corHMMDredgePrep}. When supplied, expensive data processing
    steps are skipped. This is used internally by
    \code{\link{corHMMDredge}} to avoid redundant computation across
    many model fits within the same rate category.}
}
\details{
There are three penalty types available:

\code{"l1"} penalizes the likelihood by the mean transition rate (analogous
to lasso regression). The mean rather than sum is used to avoid
over-penalizing complex models simply because they have more parameters.

\code{"l2"} penalizes by the squared mean transition rate (analogous to
ridge regression). This is the most aggressive penalty, shrinking
parameters most rapidly.

\code{"er"} penalizes by the average distance between parameter values.
An equal-rates model incurs no penalty; models with one outlier rate
incur a penalty unless strongly supported by the data. In practice this
behaves similarly to an unpenalized model.

Setting \code{pen.type = "unreg"} or \code{lambda = 0} recovers the
standard unpenalized likelihood, identical to \code{\link{corHMM}}.

When \code{use_RTMB = TRUE} the gradient is computed analytically via
automatic differentiation, which substantially improves optimisation
speed and stability for models with many parameters.

Many arguments are shared with \code{\link{corHMM}}; see that
documentation for further detail.
}
\value{
An object of class \code{"corhmm"} with elements:
\item{$loglik}{Maximum log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{AIC corrected for small sample size.}
\item{$rate.cat}{Number of rate categories.}
\item{$solution}{Matrix of maximum likelihood rate estimates.}
\item{$index.mat}{Index matrix of estimated parameters.}
\item{$data}{Input dataset.}
\item{$data.legend}{Input dataset with corHMM trait coding appended.}
\item{$phy}{Input phylogeny.}
\item{$states}{Most likely states at internal nodes.}
\item{$tip.states}{Most likely states at tips (if requested).}
\item{$states.info}{Information content (bits) at each node.}
\item{$iterations}{Number of optimiser iterations.}
\item{$root.p}{Root prior used in estimation.}
\item{$tip.fog.p}{Estimated or fixed tip fog probabilities (if
  \code{tip.fog} was specified).}
}
\examples{
\donttest{
library(corHMM)
data(primates)
phy <- multi2di(primates[[1]])
data <- primates[[2]]

# fit a single model with L1 regularization
fit <- corHMMDredgeBase(phy = phy, data = data, rate.cat = 1,
  pen.type = "l1", root.p = "maddfitz", lambda = 1)
fit

# fit a specific model structure (note that dredge doesn't allow for collapsing)
index_mat <- getStateMat4Dat(data, collapse = FALSE)$rate.mat
fit2 <- corHMMDredgeBase(phy = phy, data = data, rate.cat = 1,
  rate.mat = index_mat, pen.type = "l1", lambda = 0)
fit2
}
}
\references{
Boyko, J.D. (2026). Automatic Discovery of Optimal Discrete Character
Models. Systematic Biology, In press.

Lewis, P.O. (2001). A likelihood approach to estimating phylogeny from
discrete morphological character data. Systematic Biology, 50(6),
913-925.

Zhou, Y., Gao, M., Chen, Y., Shi, X. (2024). Adaptive Penalized
Likelihood method for Markov Chains.
}
\author{James D. Boyko}
\seealso{\code{\link{corHMMDredge}}, \code{\link{corHMM}}}
\keyword{models}
