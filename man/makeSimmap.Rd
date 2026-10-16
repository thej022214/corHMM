\name{makeSimmap}
\alias{makeSimmap}
\title{Simulate a character history}
\description{Produces a character history from a fitted corHMM object or a transition rate matrix.}
\usage{
makeSimmap(tree, data, model, rate.cat, root.p="yang", nSim=1, nCores=1, fix.node=NULL,
fix.state=NULL, parsimony = FALSE, max.attempt = 1e5, collapse=TRUE,
method=c("rejection", "uniformization"))
}       
\arguments{
\item{tree}{A phylogeny of class phylo. If omitted when \code{model} is a fitted corHMM object, uses the fit's phylogeny.}
\item{data}{a data.frame containing species information. The first column must be species names matching the phylogeny. Additional columns contain discrete character data. Omit this argument when \code{model} is a fitted corHMM object; its saved \code{data.legend} supplies the tip-state indices.}
\item{model}{A transition rate matrix or a fitted object of class \code{corhmm}. Passing the complete fit preserves the state coding used during fitting, including for fits from older package versions.}
\item{rate.cat}{The number of rate categories. Taken from \code{model} when it is a fitted corHMM object; an explicitly supplied value must agree with the fit.}
\item{root.p}{The root prior to begin the sampling at the root. If omitted when \code{model} is a fitted corHMM object, uses the fit's root prior.}
\item{nSim}{The number of simmaps to be simulated.}
\item{nCores}{The number of cores to be used.}
\item{fix.node}{A vector specifying node numbers to be fixed. Also possible to fix tips if using a hidden Markov model. Tips are in the order of tree$tip.label.}
\item{fix.state}{Specifies which states to fix the nodes. States are specified according to position in the rate matrix. E.g. If I had binary observed characters 0/1 and two hidden rate classes A/B and wanted to fix a node as 1B, I would set this to 4.}
\item{parsimony}{A boolean indicating whether node states should be based on conditional likelihood (per Bollback 2006), or if they should be consistent with a parsimonious model (if TRUE). Parsimony states are evaluted by dividing the rates present in the variable, \code{model}, by 1000 and evaluating the conditional likelihood of each state. However, by lowering the rates we can approximate a parsimony reconstruction (Steel and Penny 2000).}
\item{max.attempt}{For \code{method = "rejection"}, the maximum number of attempts to create a path between the states at the ends of a branch. When this value is reached, the Floyd-Walsh algorithm produces a shortest path with equally spaced changes. Ignored for \code{method = "uniformization"}.}
\item{collapse}{a boolean indicating whether to collapse multiple character combinations into only the observed states. For example, if true a two character dataset contained (0,0), (1,0), and (1,1), this would be collapsed into 1,2,3. However, if set to false it would 1,2,4. In combination with a custom rate matrix this allows for the estimation of transitions between the unobserved character combination. The default is \code{TRUE}. When \code{model} is a fitted corHMM object, this argument is not used: the fitted state space and saved indices are retained.}
\item{method}{The within-branch history sampler. \code{"rejection"} retains the existing behavior. \code{"uniformization"} samples paths conditional on their endpoints using a Poisson process and a discrete Markov chain, with shared precomputed matrix powers. It avoids rejected paths and the shortest-path fallback, and can be faster when the required endpoint is unlikely.}
}
\details{
This function will generate a character history given a model and dataset. It has a similar structure to the simmap generated in phytools and follows the methods of Bollback (2006). Conditional tip and node likelihoods are computed during mapping; saved ancestral reconstructions are not required. We chose not to implement any new plotting functions, instead \code{makeSimmap} produces a simmap object which is formatted so it can used with other R packages such as phytools (Revell, 2012). For additional capabilities, options, and biological examples we refer readers to the detailed _Generalized corHMM_ vignette.

For mapping an existing fit, use \code{makeSimmap(model = fit, nSim = 100)}.
The fit must retain its original \code{solution}, \code{rate.cat}, and
\code{data.legend}. Numeric root priors and \code{fix.state} continue to refer to
positions in the fitted matrix. Matrix state names are retained in the maps.

When passing a matrix instead, its rows and columns must be in the same state
order as the supplied data are assigned by the installed version of corHMM.
Column names label the output; they do not align tip states to matrix positions.
Version 2.8 sorted combined character states, while 2.10.x uses the order from
\code{expand.grid}. For two binary characters with all four combinations, these
orders are 00, 01, 10, 11 and 00, 10, 01, 11, respectively. Setting
\code{collapse = TRUE} selects the observed combinations but does not restore the
older ordering. Passing the complete fit avoids this mismatch without refitting.
Maps previously generated with mismatched state coding should be regenerated;
changing their labels alone does not correct the histories.

Uniformization follows Hobolth and Stone (2009). Virtual self-transitions are
removed from the returned histories. The Poisson expansion is truncated with
omitted conditional probability bounded by \eqn{10^{-12}}, using the computed
branch transition probabilities. Tables exceeding 10,000 virtual events or
20 million entries produce an error suggesting the rejection method. The two
methods sample the same target distribution before the rejection method's
fallback is reached; speed depends on the rates, branch lengths, and endpoints.
}
\value{
A list of simmaps.
}
\examples{
\donttest{
data(primates)
phy <- primates[[1]]
phy <- multi2di(phy)
data <- primates[[2]]

##run corhmm
MK <- corHMM(phy, data, 1)

## preserve the fitted state coding, including across package versions
simmaps <- makeSimmap(model=MK, nSim=100, nCores=1)
# To use uniformization instead, add method="uniformization".

## get the summary
simmap_summaries <- lapply(simmaps, summarize_single_simmap)
summary_df <- summarize_transition_stats(simmap_summaries)
print(summary_df)

## we import phytools plotSimmap for plotting
# library(phytools)
# plotSimmap(simmap[[1]])
}
}
\references{
Boyko, J. D., and J. M. Beaulieu. 2021. Generalized hidden Markov models for phylogenetic comparative datasets. Methods in Ecology and Evolution 12:468-478.

Bollback, J. P. 2006. SIMMAP: stochastic character mapping of discrete traits on phylogenies. BMC Bioinformatics 7:88.

Hobolth, A., and E. A. Stone. 2009. Simulation from endpoint-conditioned,
continuous-time Markov chains on a finite state space, with applications to
molecular evolution. Annals of Applied Statistics 3:1204-1231.

Revell, L. J. 2012. phytools: an R package for phylogenetic comparative biology (and other things). Methods in Ecology and Evolution, 3(2), 217-223.

Steel, M., and D. Penny. 2000. Parsimony, Likelihood, and the Role of Models in Molecular Phylogenetics. Molecular Biology and Evolution 17:839-850.
}
\author{James D. Boyko}
