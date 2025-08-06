\name{makeSimmap}
\alias{makeSimmap}
\title{Simulate a character history}
\description{Produces a character history given some of the outputs of a corHMM object.}
\usage{
makeSimmap(tree, data, model, rate.cat, root.p="yang", nSim=1, nCores=1, fix.node=NULL,
fix.state=NULL, parsimony = FALSE, max.attempt = 1e5, collapse=TRUE)
}       
\arguments{
\item{tree}{A phylogeny of class phylo.}
\item{data}{a data.frame containing species information. The first column must be species names matching the phylogeny. Additional columns contain discrete character data.}
\item{model}{The transition rate matrix.}
\item{rate.cat}{The number of rate categories.}
\item{root.p}{The root prior to begin the sampling at the root.}
\item{nSim}{The number of simmaps to be simulated.}
\item{nCores}{The number of cores to be used.}
\item{fix.node}{A vector specifying node numbers to be fixed. Also possible to fix tips if using a hidden Markov model. Tips are in the order of tree$tip.label.}
\item{fix.state}{Specifies which states to fix the nodes. States are specified according to position in the rate matrix. E.g. If I had binary observed characters 0/1 and two hidden rate classes A/B and wanted to fix a node as 1B, I would set this to 4.}
\item{parsimony}{A boolean indicating whether node states should be based on conditional likelihood (per Bollback 2006), or if they should be consistent with a parsimonious model (if TRUE). Parsimony states are evaluted by dividing the rates present in the variable, \code{model}, by 1000 and evaluating the conditional likelihood of each state. However, by lowering the rates we can approximate a parsimony reconstruction (Steel and Penny 2000).}
\item{max.attempt}{A numeric value indicating the maximum number of attempts to create a possible path between an initial and final state on a branch. When the maximum value is reached we use the Floyd-Walsh algorithm to produce the shortest path between the two states and divide the branch into equal segments.}
\item{collapse}{a boolean indicating whether to collapse multiple character combinations into only the observed states. For example, if true a two character dataset contained (0,0), (1,0), and (1,1), this would be collapsed into 1,2,3. However, if set to false it would 1,2,4. In combination with a custom rate matrix this allows for the estimation of transitions between the unobserved character combination. The default is \code{TRUE}}.
}
\details{
This function will generate a character history given a model and dataset. It has a similar structure to the simmap generated in phytools and follows the methods of Bollback (2006). If using hidden states, then it is necessary to reconstruct the tip probabilities as well as the node probabilities (i.e. get.tip.states must be \code{TRUE} when running \code{corHMM}). We chose not to implement any new plotting functions, instead \code{makeSimmap} produces a simmap object which is formatted so it can used with other R packages such as phytools (Revell, 2012). For additional capabilities, options, and biological examples we refer readers to the detailed _Generalized corHMM_ vignette.
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

##get simmap from corhmm solution
model <- MK$solution
simmaps <- makeSimmap(tree=phy, data=data, model=model, rate.cat=1, nSim=100, nCores=1)

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

Revell, L. J. 2012. phytools: an R package for phylogenetic comparative biology (and other things). Methods in Ecology and Evolution, 3(2), 217-223.

Steel, M., and D. Penny. 2000. Parsimony, Likelihood, and the Role of Models in Molecular Phylogenetics. Molecular Biology and Evolution 17:839-850.
}
\author{James D. Boyko}