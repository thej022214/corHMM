\name{makeSimmap}
\alias{makeSimmap}
\title{Simulate a character history}
\description{Produces a character history given some of the outputs of a corHMM object.}
\usage{
makeSimmap(tree, data, model, rate.cat, root.p="yang", nSim=1, nCores=1, fix.node=NULL, fix.state=NULL)
}       
\arguments{
\item{tree}{A phylogeny of class phylo.}
\item{data}{a data.frame containing species information. The first column must be species names matching the phylogeny. Additional columns contain discrete character data.}
\item{model}{The transition rate matrix.}
\item{rate.cat}{The number of rate categories.}
\item{root.p}{The root prior to begin the sampling at the root. Currently only "yang" allowed.}
\item{nSim}{The number of simmaps to be simulated.}
\item{nCores}{The number of cores to be used.}
\item{fix.node}{Specifies the node number to be fixed. Also possible to fix tips if using a hidden Markov model. Tips are in the order of tree$tip.label.}
\item{fix.state}{Specifies which state to fix the node. State must be specified according to position in the rate matrix. E.g. If I had binary observed characters 0/1 and two hidden rate classes A/B and wanted to fix a node as 1B, I would set this to 4.}
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
simmap <- makeSimmap(tree=phy, data=data, model=model, rate.cat=1, nSim=1, nCores=1)

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