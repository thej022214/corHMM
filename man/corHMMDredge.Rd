\name{corHMMDredge}
\alias{corHMMDredge}
\title{Automatic Discovery of Optimal Discrete Character Models via Simulated Annealing}
\description{
This function automates the search for the best-fitting model of discrete character evolution. It uses a combination of simulated annealing and regularization to explore a vast space of possible model structures. This approach helps to identify optimal models without requiring the user to manually specify and compare a large set of candidate models, which is often impractical for complex datasets with multiple characters or hidden rate categories.
}
\usage{
corHMMDredge(phy, data, max.rate.cat=1, init.rate.cat=1,
  root.p="maddfitz", tip.fog=NULL, fog.ip=0.01, pen.type = "l1", lambda = 0,
  drop.threshold = 1e-7, criterion="AIC", merge.threshold=0, index_mat=NULL,
  node.states = "none", fixed.nodes=FALSE, ip=NULL, nstarts=1, n.cores=1,
  get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = FALSE,
  lower.bound = 1e-10, upper.bound = 100, opts=NULL, verbose=TRUE,
  p=NULL, rate.cat=NULL, use_RTMB=TRUE, max.iterations = 200,
  initial.temp = 2, cooling.rate = 0.95, temp.schedule = "exponential",
  seed = NULL, return.all=FALSE, checkpoint.file=NULL,
  checkpoint.interval=50)
}
\arguments{
  \item{phy}{A phylogenetic tree of class \code{"phylo"}.}
  \item{data}{A data frame with two or more columns: the first containing
    species names and the second onwards containing the discrete character
    states.}
  \item{max.rate.cat}{The maximum number of hidden rate categories to
    explore. The search runs sequentially from \code{init.rate.cat} up to
    this value.}
  \item{init.rate.cat}{The initial number of hidden rate categories to
    begin the search. Defaults to 1.}
  \item{root.p}{Specifies the prior on the root state probabilities. Can
    be \code{"maddfitz"}, \code{"yang"}, a numeric vector of
    probabilities, or \code{NULL} (see \code{\link{corHMM}}).}
  \item{tip.fog}{Optional specification of tip observation error
    (misclassification probabilities). If a numeric value less than 1 is
    supplied it is treated as a fixed per-state fog probability; if
    \code{sum(tip.fog) >= 1} the fog probabilities are estimated during
    optimisation.}
  \item{fog.ip}{Starting value for tip fog probability estimation.
    Defaults to 0.01.}
  \item{pen.type}{The type of regularization penalty to apply. Options
    are \code{"l1"} (lasso), \code{"l2"} (ridge), and \code{"er"} (equal
    rates penalty). Defaults to \code{"l1"}.}
  \item{lambda}{A hyper-parameter (0 to 1) controlling the strength of
    the regularization penalty. A value of 0 means no penalization.}
  \item{drop.threshold}{The rate value below which a parameter becomes a
    candidate for being dropped (fixed to zero) during a drop move.}
  \item{criterion}{The information criterion used to evaluate and compare
    models during the search. Either \code{"AIC"} or \code{"AICc"}.
    Defaults to \code{"AIC"}.}
  \item{merge.threshold}{The absolute difference between two rate
    estimates below which they become candidates for merging (constrained
    equal) during a merge move.}
  \item{index_mat}{An optional user-specified starting rate matrix. If
    \code{NULL} the search starts from a fully parameterized (all rates
    different) model.}
  \item{node.states}{The method for ancestral state reconstruction.
    Options are \code{"marginal"}, \code{"joint"}, \code{"scaled"}, or
    \code{"none"}.}
  \item{fixed.nodes}{Logical. Whether to fix the states of labelled
    internal nodes in the provided phylogeny.}
  \item{ip}{Optional numeric vector of initial parameter values for the
    optimiser.}
  \item{nstarts}{The number of random restarts to perform during
    optimisation for each candidate model. Higher values improve the
    chance of finding the global likelihood peak at the cost of
    computation time.}
  \item{n.cores}{The number of cores to use for parallel random restarts
    within a single model fit.}
  \item{get.tip.states}{Logical. Whether to estimate marginal state
    probabilities at the tips.}
  \item{lewis.asc.bias}{Logical. Whether to correct for ascertainment
    bias following Lewis (2001).}
  \item{collapse}{Logical. Whether to collapse unobserved character state
    combinations.}
  \item{lower.bound, upper.bound}{Lower and upper bounds for rate
    parameter estimates.}
  \item{opts}{A list of control options passed to the \code{nloptr}
    optimiser.}
  \item{verbose}{Logical. Whether to print search progress to the
    console.}
  \item{p, rate.cat}{Used together to evaluate a single fixed model. If
    both are provided the dredge search is skipped and a single model fit
    is returned.}
  \item{use_RTMB}{Logical. Whether to use the RTMB automatic
    differentiation backend for gradient-based optimisation. Defaults to
    \code{TRUE}. When \code{TRUE} the gradient is computed analytically,
    which is substantially faster for complex models.}
  \item{max.iterations}{The maximum number of simulated annealing
    iterations within each rate category.}
  \item{initial.temp}{The starting temperature for the simulated
    annealing algorithm. Higher values allow broader exploration by
    accepting worse models with higher probability early in the search.}
  \item{cooling.rate}{The multiplicative rate at which temperature
    decreases per iteration (\eqn{T_{new} = T_{old} \times
    \text{cooling.rate}}). Must be between 0 and 1. Defaults to 0.95.}
  \item{temp.schedule}{The cooling schedule. Options are
    \code{"exponential"} and \code{"linear"}.}
  \item{seed}{An integer random seed for reproducibility.}
  \item{return.all}{Logical. If \code{TRUE}, returns a list containing
    both the pruned model set and full simulated annealing diagnostics
    (required for \code{\link{plotDredgeTrace}}).}
  \item{checkpoint.file}{Optional path to an RDS file used for
    checkpointing. If the file already exists the search resumes from the
    saved state; otherwise a new checkpoint is created at regular
    intervals.}
  \item{checkpoint.interval}{Number of iterations between checkpoint
    saves. Defaults to 50.}
}
\details{
This function implements an automated model discovery framework for
phylogenetic comparative methods. The core of the method is a simulated
annealing heuristic that explores the space of discrete character
evolution models combined with optional regularization.

\strong{Proposal moves.} At each iteration the algorithm proposes one of
four modifications to the current model's rate matrix index structure:
\enumerate{
  \item \strong{Drop:} A rate parameter is removed from the model (fixed
    to zero). Parameters with smaller estimated values are preferentially
    targeted.
  \item \strong{Merge:} Two or more rate parameters are constrained to be
    equal. Pairs with similar estimated values are preferentially
    targeted.
  \item \strong{Free:} A currently constrained or dropped parameter is
    freed to be estimated independently, expanding model complexity.
  \item \strong{Eigen-merge:} Parameters are merged based on the
    eigenstructure of the current rate matrix rather than raw value
    proximity. Specifically, each free parameter is represented as a
    profile of its absolute loadings onto the eigenvectors of Q. Parameters
    whose profiles have high cosine similarity -- meaning they drive
    transitions along the same directions in state space -- are identified
    as functionally redundant and proposed for merging. This move can
    identify candidates that value-based merging would miss, and is
    particularly effective for complex multi-state or multi-rate-category
    models. For models with multiple rate categories the full block rate
    matrix is used, so eigen-merge can also propose merging parameters
    across rate categories when they are functionally equivalent.
}

\strong{Acceptance criterion.} A proposed model that improves the
information criterion score is always accepted. A worse model may still
be accepted with probability \eqn{\exp(-\Delta / T)}, where \eqn{\Delta}
is the score difference and \eqn{T} is the current temperature. As
temperature cools the algorithm transitions from broad exploration to
focused exploitation.

\strong{Restart mechanism.} To avoid entrapment in local optima the
algorithm can restart from a previously accepted model when the search
has stagnated. The restart strategy is controlled internally and the
restart event is recorded in the trace for diagnostic purposes.

\strong{Memory efficiency.} During the search, \code{phy}, \code{data},
and \code{data.legend} are stripped from stored model objects to reduce
RAM usage. These fields are restored on the small set of models that
survive \code{prune_redundant} before the function returns, so every
model in the output is a complete \code{corhmm} object.

\strong{Regularization.} An optional L1 or L2 penalty on the rate
parameters is incorporated directly into the likelihood during
optimisation. The \code{lambda} parameter controls the penalty strength.
Setting \code{lambda = 0} recovers standard maximum likelihood.
}
\value{
By default, an object of class \code{"corhmm.dredge"}, which is a list
of unique accepted models pruned to remove redundant entries (identical
likelihoods and parameter counts). Each element is a standard
\code{corhmm} object including \code{$phy}, \code{$data}, and
\code{$data.legend}. The list is ordered by the information criterion
score, best model first.

If \code{return.all = TRUE}, a list with two elements:
\describe{
  \item{all_models}{The pruned \code{"corhmm.dredge"} object described
    above.}
  \item{sa_fits}{A list of detailed simulated annealing outputs for each
    rate category explored, including per-iteration scores, move types,
    and acceptance decisions. Required by \code{\link{plotDredgeTrace}}.}
}
}
\author{
James D. Boyko
}
\seealso{
\code{\link{corHMM}}, \code{\link{plotDredgeTrace}},
\code{\link{kFoldCrossValidation}}, \code{\link{getModelTable}}
}
\examples{
\donttest{
library(corHMM)
data(primates)
phy <- multi2di(primates[[1]])
phy$edge.length <- phy$edge.length + 1e-6
data <- primates[[2]]

# Basic dredge search -- single rate category, L1 penalty
dredge_fits <- corHMMDredge(phy = phy,
  data = data,
  max.rate.cat = 1,
  pen.type = "l1",
  root.p = "maddfitz",
  lambda = 1,
  return.all = TRUE)

# return.all = TRUE is required for the trace plot
plotDredgeTrace(dredge_fits)

# Inspect the model table
mod_table <- getModelTable(dredge_fits$all_models)
print(mod_table)

# The best model for downstream analysis
best_fit <- dredge_fits$all_models[[which.min(mod_table$dAIC)]]

# Without penalization (lambda = 0) likelihoods are directly comparable
# with standard corHMM output
# dredge_fits_ml <- corHMMDredge(phy = phy,
#   data = data,
#   max.rate.cat = 1,
#   pen.type = "l1",
#   root.p = "maddfitz",
#   lambda = 0)

# Resume an interrupted run using checkpointing
# dredge_fits <- corHMMDredge(phy = phy,
#   data = data,
#   max.rate.cat = 1,
#   lambda = 1,
#   checkpoint.file = "dredge_checkpoint.rds",
#   checkpoint.interval = 25)
}
}
\references{
Boyko, J.D. (2026). Automatic Discovery of Optimal Discrete Character
Models. Systematic Biology, In press.

Lewis, P.O. (2001). A likelihood approach to estimating phylogeny from
discrete morphological character data. Systematic Biology, 50(6),
913-925.

Pagel, M., & Meade, A. (2006). Bayesian analysis of correlated evolution
of discrete characters by reversible-jump Markov chain Monte Carlo. The
American Naturalist, 167(6), 808-825.

Kalyaanamoorthy, S., Minh, B.Q., Wong, T.K.F., von Haeseler, A., &
Jermiin, L.S. (2017). ModelFinder: fast model selection for accurate
phylogenetic estimates. Nature Methods, 14, 587-589.

Burnham, K.P., & Anderson, D.R. (2002). Model selection and multimodel
inference: a practical information-theoretic approach. Springer-Verlag.

Zhou, Y., Gao, M., Chen, Y., Shi, X. (2024). Adaptive Penalized
Likelihood method for Markov Chains.
}
