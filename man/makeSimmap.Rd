\name{makeSimmap}
\alias{makeSimmap}
\title{Simulate a character history}
\description{Produces a character history given some of the outputs of a corHMM object.}
\usage{
makeSimmap(tree, tip.states, states, model, nSim=1, nCores=1)
}       
\arguments{
\item{tree}{A phylogeny of class phylo.}
\item{tip.states}{A data.frame where each column represents the presence or absence of a particular state and the rownames are the species names (see example below).}
\item{states}{A data.frame where each column represents the marginal probability of a particular state at a node and the rownames are the species names (see example below).}
\item{model}{The transition rate matrix.}
\item{nSim}{The number of simmaps to be simulated.}
\item{nCores}{The number of cores to be used.}
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

##get simmap inputs from corhmm outputs
phy <- MK$phy
data <- MK$data
model <- MK$solution
model[is.na(model)] <- 0
diag(model) <- -rowSums(model)
states <- MK$states
tip.states <- MK$tip.states

## run get simmap (can be plotted using phytools)
simmap <- makeSimmap(tree=phy, tip.states=tip.states, states=states, model=model, 
	nSim=1, nCores=1)

## we import phytools plotSimmap for plotting
plotSimmap(simmap[[1]])
}
}
\references{
Bollback, J. P. 2006. SIMMAP: stochastic character mapping of discrete traits on phylogenies. BMC Bioinformatics 7:88.

Revell, L. J. 2012. phytools: an R package for phylogenetic comparative biology (and other things). Methods in Ecology and Evolution, 3(2), 217-223.
}
\author{James D. Boyko}
\keyword{models}