\name{corHMMDredge}
\alias{corHMMDredge}
\title{Automatic Discovery of Optimal Discrete Character Models via Simulated Annealing}
\description{
This function automates the search for the best-fitting model of discrete character evolution. It uses a combination of simulated annealing and regularization to explore a vast space of possible model structures. This approach helps to identify optimal models without requiring the user to manually specify and compare a large set of candidate models, which is often impractical for complex datasets with multiple characters or hidden rate categories.
}
\usage{
corHMMDredge(phy, data, max.rate.cat=1, init.rate.cat=1, 
  root.p="maddfitz", pen.type = "l1", lambda = 0, drop.threshold = 1e-7, 
  criterion="AIC", merge.threshold=0, index_mat=NULL, node.states = "marginal", 
  fixed.nodes=FALSE, ip=NULL, nstarts=0, n.cores=1, get.tip.states = FALSE, 
  lewis.asc.bias = FALSE, collapse = FALSE, lower.bound = 1e-10, 
  upper.bound = 100, opts=NULL, verbose=TRUE, p=NULL, rate.cat=NULL, grad=FALSE, 
  max.iterations = 200, initial.temp = 2, cooling.rate = 0.95, 
  temp.schedule = "exponential", seed = NULL, return.all=FALSE)
}
\arguments{
  \item{phy}{A phylogenetic tree of class 'phylo'.}
  \item{data}{A data frame with two or more columns: the first containing species names and the second onwards containing the discrete character states.}
  \item{max.rate.cat}{The maximum number of hidden rate categories to explore. The dredge will run sequentially from \code{init.rate.cat} up to this value.}
  \item{init.rate.cat}{The initial number of hidden rate categories to begin the search. Defaults to 1.}
  \item{root.p}{Specifies the prior on the root state probabilities. Can be "maddfitz", "yang", a vector of probabilities, or NULL (see \code{\link{corHMM}}).}
  \item{pen.type}{The type of regularization penalty to apply. Options are "l1" (lasso), "l2" (ridge), and "er" (equal rates penalty). Defaults to "l1".}
  \item{lambda}{A hyper-parameter (from 0 to 1) that controls the severity of the regularization penalty. A value of 0 means no penalization.}
  \item{drop.threshold}{The rate value below which a parameter becomes a candidate for being dropped (set to 0) during a "drop" move in the simulated annealing search.}
  \item{criterion}{The information criterion to use for model evaluation during the search. "AIC" (Akaike Information Criterion) or "AICc" (AIC corrected for small sample sizes). Defaults to "AIC".}
  \item{merge.threshold}{The absolute difference between two rate estimates below which they become candidates for being merged (constrained to be equal) during a "merge" move.}
  \item{index_mat}{An optional, user-specified rate matrix to begin the search. If NULL, the search starts with a fully parameterized (all rates different) model.}
  \item{node.states}{The method for ancestral state reconstruction. Options are "marginal", "joint", "scaled", or "none".}
  \item{fixed.nodes}{Logical, indicating whether to fix the states of nodes with labels in the provided phylogeny.}
  \item{ip}{Initial parameter values for the optimization routine.}
  \item{nstarts}{The number of random restarts to perform during optimization for each model visited. Helps in finding the global likelihood peak.}
  \item{n.cores}{The number of cores to use for parallel processing of random restarts.}
  \item{get.tip.states}{Logical, whether to estimate the marginal probabilities of states at the tips (useful for cross-validation).}
  \item{lewis.asc.bias}{Logical, whether to correct for ascertainment bias as described by Lewis (2001).}
  \item{collapse}{Logical, whether to collapse character states that are not observed in the data.}
  \item{lower.bound, upper.bound}{The lower and upper bounds for rate parameter estimates.}
  \item{opts}{A list of control options for the \code{nloptr} optimization routine.}
  \item{verbose}{Logical, whether to print progress of the search to the console.}
  \item{p, rate.cat}{Used to fit a single, fixed model. If \code{p} (a vector of parameters) and \code{rate.cat} are provided, the function will not perform the dredge search and instead will return a single model fit.}
  \item{grad}{Logical, whether to use the gradient-based optimizer "NLOPT_LD_MMA".}
  \item{max.iterations}{The maximum number of iterations for the simulated annealing search within each rate category.}
  \item{initial.temp}{The starting "temperature" for the simulated annealing algorithm. Higher values allow the search to more freely explore the model space, accepting worse models with a higher probability.}
  \item{cooling.rate}{The rate at which the temperature decreases in the simulated annealing search (e.g., T_new = T_old * cooling.rate). Must be between 0 and 1. Defaults to 0.95.}
  \item{temp.schedule}{The cooling schedule for the temperature. Currently only "exponential" is implemented.}
  \item{seed}{An integer to set the random seed for reproducibility.}
  \item{return.all}{Logical. If TRUE, returns a detailed list including all visited models and search diagnostics, not just the pruned final set.}
}
\details{
This function implements an automated model discovery framework for phylogenetic comparative methods based on the work of Boyko (2025). The core of the method is a simulated annealing heuristic search algorithm combined with regularization.

The algorithm explores the vast space of possible model structures by iteratively proposing modifications to the current model. Three types of moves are used:
\enumerate{
  \item \strong{Drop:} A rate parameter is removed from the model (i.e., fixed to zero). Smaller rates are more likely to be dropped.
  \item \strong{Merge:} Two rate parameters are constrained to be equal. Pairs with more similar estimated rates are more likely to be merged.
  \item \strong{Free:} A currently constrained (dropped or merged) parameter is freed, allowing it to be estimated independently.
}

The decision to accept a new model is probabilistic. A proposed model that improves the information criterion score (e.g., lower AIC) is always accepted. A model that is worse may still be accepted with a probability that depends on the current "temperature" of the system. The temperature starts high, allowing broad exploration, and is gradually lowered according to a cooling schedule. As the temperature approaches zero, the algorithm behaves like a greedy hill-climbing search, converging on an optimal solution. To avoid getting stuck in local optima, a restart mechanism is included.

Regularization (L1 or L2) is incorporated into the likelihood calculation. This adds a penalty term for model complexity (either the magnitude or squared magnitude of rate parameters), which helps to prevent overfitting and improve the generalizability of the final model. The \code{lambda} parameter tunes the strength of this penalty.
}
\value{
By default, an object of class \code{corhmm.dredge}, which is a list containing all unique models that were accepted during the simulated annealing search. The models are pruned to remove redundant entries (i.e., models with identical likelihoods and parameter counts but different rate structures). Each element of the list is a standard \code{corhmm} object. The list is sorted by the information criterion, with the best model first.

If \code{return.all=TRUE}, a list containing two elements:
  \item{all_models}{The pruned \code{corhmm.dredge} object described above.}
  \item{sa_fits}{A detailed list containing the raw output from the simulated annealing search for each rate category explored.}
}
\author{
James D. Boyko
}
\seealso{
\code{\link{corHMM}}, \code{\link{kFoldCrossValidation}}
}
\examples{
\donttest{
library(corHMM)
data(primates)
phy <- multi2di(primates[[1]])
phy$edge.length <- phy$edge.length + 1e-6
data <- primates[[2]]
# fit the models following the same input style as corHMM. 
# here we are NOT going to look for multiple rate classes (max.rate.cat=1)
dredge_fits <- corHMMDredge(phy = phy, 
  data = data,
  max.rate.cat = 1, 
  pen.type = "l1", 
  root.p = "maddfitz", 
  lambda = 1, 
  return.all=TRUE)
# return.all must be true to get the trace plot
plotDredgeTrace(dredge_fits)
# produce a model table
mod_table <- getModelTable(dredge_fits$all_models)
print(mod_table)
# which ever model is best is the one used for downstream analysis
best_fit <- dredge_fits[[which.min(mod_table$dAIC)]]
best_fit

# you can also fit dredge without any penalization. 
# this will make the likelihoods directly comparable with corHMM- just set lambda to 0
dredge_fits_no_penalty <- corHMMDredge(phy = phy, 
  data = data, 
  max.rate.cat = 1, 
  pen.type = "l1", 
  root.p = "maddfitz", 
  lambda = 0)
}
}
\references{
Boyko, J.D. (2025). Automatic Discovery of Optimal Discrete Character Models.

Pagel, M., & Meade, A. (2006). Bayesian analysis of correlated evolution of discrete characters by reversible-jump Markov chain Monte Carlo. The American Naturalist, 167(6), 808-825.

Burnham, K. P., & Anderson, D. R. (2002). Model selection and multimodel inference: a practical information-theoretic approach. Springer-Verlag.

Zhou, Y., Gao, M., Chen, Y., Shi, X., 2024. Adaptive Penalized Likelihood method for Markov Chains. 
}