\name{corHMMDredge}
\alias{corHMMDredge}
\title{Automatic Discovery of Optimal Discrete Character Models via Simulated Annealing}
\description{
This function automates the search for the best-fitting model of discrete
character evolution. It uses a combination of simulated annealing and
optional regularization to explore a vast space of possible model
structures, identifying optimal models without requiring the user to
manually specify and compare candidate models.
}
\usage{
corHMMDredge(phy, data, max.rate.cat=1, init.rate.cat=1,
  root.p="maddfitz", tip.fog=NULL, fog.ip=0.01, pen.type="l1", lambda=0,
  drop.threshold=1e-7, criterion="AIC", merge.threshold=0, index_mat=NULL,
  node.states="none", fixed.nodes=FALSE, ip=NULL, nstarts=1, n.cores=1,
  get.tip.states=FALSE, lewis.asc.bias=FALSE, collapse=FALSE,
  lower.bound=1e-10, upper.bound=100, opts=NULL, verbose=TRUE,
  p=NULL, rate.cat=NULL, use_RTMB=TRUE, max.iterations=1000,
  initial.temp=2, cooling.rate=0.95, temp.schedule="exponential",
  seed=NULL, checkpoint.file=NULL, checkpoint.interval=50)
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
    \code{"none"}. Defaults to \code{"none"} during the search for speed;
    set to \code{"marginal"} if ancestral states are needed on the
    returned models.}
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
    is returned immediately.}
  \item{use_RTMB}{Logical. Whether to use the RTMB automatic
    differentiation backend for gradient-based optimisation. Defaults to
    \code{TRUE}. When \code{TRUE} the gradient is computed via autodifferentiation,
    which is substantially faster for complex models.}
  \item{max.iterations}{The maximum number of simulated annealing
    iterations within each rate category.}
  \item{initial.temp}{The starting temperature for the simulated
    annealing algorithm. Higher values allow broader exploration by
    accepting worse models with higher probability early in the search.}
  \item{cooling.rate}{The multiplicative rate at which temperature
    decreases per iteration (\eqn{T_{new} = T_{old} \times
    \mathrm{cooling.rate}}). Must be between 0 and 1. Defaults to 0.95.}
  \item{temp.schedule}{The cooling schedule. Options are
    \code{"exponential"} and \code{"linear"}.}
  \item{seed}{An integer random seed for reproducibility.}
  \item{checkpoint.file}{Optional path to an RDS file used for
    checkpointing. If the file already exists the search resumes from the
    saved state; otherwise a new checkpoint is created at regular
    intervals. When \code{max.rate.cat > 1} a separate checkpoint file is
    written per rate category (e.g. \code{ckpt_rc1.rds},
    \code{ckpt_rc2.rds}).}
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
    freed to be estimated independently, expanding model complexity[1].
    
  \item \strong{Eigen-merge:} Parameters are merged based on the
    eigenstructure of the current rate matrix rather than raw value
    proximity[2]. 
    
\strong{[1]} Only positions permitted by the maximum index matrix are eligible,
    preventing proposals that violate the structural constraints of the
    model space.
    
\strong{[2]} Each free parameter is represented as a profile of its
    absolute loadings onto the eigenvectors of Q. Parameters whose
    profiles have high cosine similarity -- meaning they drive transitions
    along the same directions in state space -- are identified as
    functionally redundant and proposed for merging. This move can
    identify candidates that value-based merging would miss, and is
    particularly effective for complex multi-state or multi-rate-category
    models.
}

\strong{Acceptance criterion.} A proposed model that improves the
information criterion score is always accepted. A worse model may still
be accepted with probability \eqn{\exp(-\Delta / T)}, where \eqn{\Delta}
is the score difference and \eqn{T} is the current temperature. As
temperature cools the algorithm transitions from broad exploration to
focused exploitation.

\strong{Restart mechanism.} To avoid entrapment in local optima the
algorithm periodically restarts from a previously accepted model when the
search has stagnated (no improvement in \code{restart.interval}
iterations). The restart event is recorded in the trace for diagnostic
purposes.

\strong{Regularization.} An optional L1 or L2 penalty on the rate
parameters is incorporated directly into the likelihood during
optimisation. The \code{lambda} parameter controls the penalty strength.
Setting \code{lambda = 0} recovers standard maximum likelihood, making
scores directly comparable with \code{\link{corHMM}}.
}
\value{
An object of class \code{"corhmm.dredge"}, which is a list of unique
accepted models pruned to remove redundant entries (identical likelihoods
and parameter counts). Each element is a standard \code{corhmm} object
including \code{$phy}, \code{$data}, and \code{$data.legend}. The list
is ordered by the information criterion score, best model first.

Full simulated annealing diagnostics for each rate category explored
(per-iteration scores, move types, acceptance decisions) are stored as a
hidden attribute and can be accessed via
\code{attr(result, "dredge_history")}. This attribute is used
automatically by \code{\link{plotDredgeTrace}} and is displayed in the
\code{print} summary.
}
\author{
James D. Boyko
}
\seealso{
\code{\link{corHMM}}, \code{\link{plotDredgeTrace}},
\code{\link{kFoldCrossValidation}}, \code{\link{getModelTable}}
}
\examples{
\dontrun{
library(corHMM)
data(primates)
phy <- multi2di(primates[[1]])
phy$edge.length <- phy$edge.length + 1e-6
data <- primates[[2]]

# basic dredge search -- single rate category, L1 penalty
dredge_fits <- corHMMDredge(phy = phy,
  data = data,
  max.rate.cat = 1,
  pen.type = "l1",
  root.p = "maddfitz",
  lambda = 1)

# printing shows a model table and the best model
print(dredge_fits)

# trace plot uses the hidden dredge_history attribute automatically
plotDredgeTrace(dredge_fits)

# inspect the model table directly
mod_table <- getModelTable(dredge_fits)
print(mod_table)

# the best model for downstream analysis
best_fit <- dredge_fits[[which.min(mod_table$dAIC)]]

# access full SA diagnostics if needed
history <- attr(dredge_fits, "dredge_history")

# note: models inside dredge_history are stripped of $phy, $data, and
# $data.legend to save memory during the search. if you need complete
# objects, reattach the shared fields manually:
for(i in seq_along(history)){
  for(j in seq_along(history[[i]]$accepted_models)){
    history[[i]]$accepted_models[[j]]$phy         <- dredge_fits[[1]]$phy
    history[[i]]$accepted_models[[j]]$data        <- dredge_fits[[1]]$data
    history[[i]]$accepted_models[[j]]$data.legend <- dredge_fits[[1]]$data.legend
  }
}

for(i in seq_along(history)){
  for(j in seq_along(history[[i]]$every_model)){
    if(is.null(history[[i]]$every_model[[j]])){next} 
    history[[i]]$every_model[[j]]$phy         <- dredge_fits[[1]]$phy
    history[[i]]$every_model[[j]]$data        <- dredge_fits[[1]]$data
    history[[i]]$every_model[[j]]$data.legend <- dredge_fits[[1]]$data.legend
  }
}

# resume an interrupted run using checkpointing
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
