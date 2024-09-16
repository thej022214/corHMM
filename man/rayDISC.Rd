\name{rayDISC}
\alias{rayDISC}
\alias{dev.raydisc}
\title{Evolution of categorical traits}
\description{Fits a model of evolution for categorical traits, allowing for multi-state characters, polymorphisms, missing data, and incompletely resolved trees}
\usage{
rayDISC(phy, data, charnum, lb, ub, ...)
}
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
\item{data}{a data matrix containing species information (see Details).}
\item{charnum}{integer index + 1 of the character to analyze.}
\item{lb}{lower bound for the likelihood search.}
\item{ub}{upper bound for the likelihood search.}
\item{...}{additional arguments passed to \code{corHMM}.}
}
\details{
This function serves as a wrapper for \code{corHMM} function with 
\code{rate.cat = 1}. All functionality of \code{rayDISC} is contained within 
the \code{corHMM} function (see _Generalized corHMM_ vignette), but we still 
provide it for those that are more comfortable using it. 

See \code{\link{corHMM}} for full parameter details.

\code{rayDISC} takes a tree and a trait file and estimates transition rates and ancestral states for binary or multistate characters. The first column of the trait file must contain the species labels to match to the tree, with the second, third, fourth, and so on, corresponding to the traits of interest. The input phylogeny need not be bifurcating as the algorithm is implemented to handle multifucations. Polytomies are allowed by generalizing Felsenstein's (1981) pruning algorithm to be the product of the probability of observing the tip states of n descendant nodes, rather than two, as in the completely bifurcating case.

The user can fix the root state probabilities by supplying a vector to the \code{root.p}. If the user supplies the flag \code{root.p}=\dQuote{yang}, then the estimated transition rates are used to set the weights at the root (see pg. 124 Yang 2006), whereas specifying \code{root.p}=\dQuote{maddfitz} employs the same procedure described by Maddison et al. (2007) and FitzJohn et al. (2009). Note that the default \code{root.p=NULL} assumes equal weighting among all possible states.

Ambiguities (polymorphic taxa or taxa missing data) are assigned likelihoods following Felsenstein (2004, p. 255). Taxa with missing data are coded \dQuote{?} with all states observed at a tip. Polymorphic taxa are coded with states separated by an \dQuote{&}. For example, if a trait has four states and taxonA is observed to be in state 1 and 3, the character would be coded as \dQuote{1&3}. corHMM then uses this information to assign a likelihood of 1.0 to both states. Missing data are treated as ambiguous for all states, thus all states for taxa missing data are assigned a likelihood of 1.0. For example, for a four-state character (i.e. DNA), a taxon missing data will have likelihoods of all four states equal to 1.0 [e.g. L(A)=1.0, L(C)=1.0, L(G)=1.0, L(T)=1.0].

In all ancestral state reconstruction implementations, the rates are first estimated, and subsequently, the MLE estimates of the rates are used to determine either the state probabilities (i.e., marginal or "scaled") or maximum likelihood states at nodes. This is the default -- i.e., the state.recon="subsequently" argument. However, for this function only, we also allow for both rates and states to be estimated jointly. This can be done with state.recon="estimate". We also allow for a hypothesis about states at all or even some nodes to help fixed, with the rates (and in some cases some of the states) being estimated. This is state.recon="given". For more information please see Vignette "Getting Likelihoods From Reconstructions".
}
\value{
An object of class \code{\link{corHMM}}.
}
\seealso{\code{\link{plotRECON}}}
\examples{
\donttest{
### Example 1
data(rayDISC.example)
## Perform ancestral state estimation, using an asymmetric model of evolution and marginal
## reconstruction of ancestral states
recon <- rayDISC(phy = rayDISC.example$tree, data = rayDISC.example$trait, model="ARD", node.states="marginal")

## Plot reconstructions on tree
plotRECON(rayDISC.example$tree, recon$states)

### Example 2
## Perform ancestral state estimation on the second character, using a single-rate model of
## evolution, marginal reconstruction of ancestral states, and setting the lower bound for
## parameter estimates to 0.01
recon <- rayDISC(rayDISC.example$tree, rayDISC.example$trait, charnum = 2, model="ER",
node.states="marginal",lb=0.01)

### Example 3
## Perform ancestral state estimation on the third character, using a single-rate model of
## evolution and joint reconstruction of ancestral states
recon <- rayDISC(rayDISC.example$tree,rayDISC.example$trait, charnum = 3, 
model="ER",node.states="joint")
}
}
\references{
Felsenstein, J. 1981. A likelihood approach to character weighting and what it tells us about parsimony and compatibility. Biological Journal of the Linnean Society 16: 183-196.

Felsenstein J. 2004. Inferring phylogenies. Sunderland MA: Sinauer Associates.

FitzJohn, R.G., W.P. Maddison, and S.P. Otto. 2009. Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Systematic Biology 58:595-611.

Lewis, P.O. 2001. A likelihood approach to estimating phylogeny from discrete morphological character data. Systematic Biology 50:913-925.

Maddison, W.P., P.E. Midford, and S.P. Otto. 2007. Estimating a binary characters effect on speciation and extinction. Systematic Biology 56:701-710.
}
\author{Jeffrey C. Oliver and Jeremy M. Beaulieu}
\keyword{models}
