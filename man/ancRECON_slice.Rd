\name{ancRECON_slice}
\alias{ancRECON_slice}
\title{Ancestral state reconstruction for a particular time}
\description{Infers marginal ancestral states based on a set of model parameters at a particular time}
\usage{
ancRECON_slice(corhmm.obj, time_slice, collapse=TRUE, ncores = 1)
}       
\arguments{
\item{corhmm.obj}{a corhmm object which is the output of the main corhmm function.}
\item{time_slice}{a vector of times to reconstruct (present = 0, root = max(branching.times(phy)))}
\item{collapse}{set collapse to be the same as it was during the corhmm run.}
\item{ncores}{number of cores to use during parallelization.}
}
\details{
This is a stand alone function for computing the marginal likelihoods at particular points along a branch for a given set of transition rates. ancRECON has the technical details of ancestral state reconstruction if you're interested. The time_slice argument will specify a time point where marginal reconstructions will be produced. You can imagine the time slice intersecting the branches of the phylogeny and doing a reconstruction there rather than at nodes as is typically done. 
}
\value{
a data.frame. Each row of time_slice is the time period that was specified. node is the closest tipward node to the slice on a particular branch. position is the amount of time (rootward) from the nearest node. Remaining columns are the marginal probabilities of each state at that particular node. There is also a plotting function that I don't currently export because it's unfinished. The example shows how to access it. 
}
\examples{
\donttest{
# library(corHMM)
# library(viridis)
# data(primates)
# phy <- multi2di(primates[[1]])
# data <- primates[[2]]
# corhmm.obj <- corHMM(phy = phy, data = data, rate.cat = 1)
# test <- ancRECON_slice(corhmm.obj, time_slice = c(1, 10, 20, 30, 40, 49), collapse = TRUE, ncores = 4)
# corHMM:::plot_slice_recon(phy, test, col = viridis(3))
}
}
\references{
Yang, Z. 2006. Computational Molecular Evolution. London:Oxford.

Boyko, J. D., and J. M. Beaulieu. 2021. Generalized hidden Markov models for phylogenetic comparative datasets. Methods in Ecology and Evolution 12:468-478.
}

\author{James D. Boyko}
\keyword{reconstructions}

