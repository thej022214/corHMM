\name{lewisMkv}
\alias{lewisMkv}
\title{Phylogenetic inference of morphology using mkv model}
\description{An implementation of Lewis' (2001) mkv model to infer topology and branch lengths based on morphological characters.}
\usage{
lewisMkv(phy, data, include.gamma=FALSE, ngammacats=4, include.beta=FALSE, 
exclude.sites=NULL, max.tol=.Machine$double.eps^0.25, ncores=NULL)
}       
\arguments{
\item{phy}{a fixed rooted or unrooted topology in \code{ape} \dQuote{phylo} format.}
\item{data}{a data frame containing character alignment (see Details).}
\item{include.gamma}{a logical indicating whether a discrete gamma should be used to allow rates to vary across sites.}
\item{ngammacats}{indicates how many discrete gamma categories to use.} 
\item{include.beta}{a logical indicating whether a beta distribution should be used to allow for unequal equilibrium frequencies among states. NOT YET IMPLEMENTED.}
\item{exclude.sites}{a vector indicating which sites to exclude from the analysis.}
\item{max.tol}{supplies the relative optimization tolerance to \code{nlopt}.}
\item{ncores}{specifies the number of independent processors to conduct the analysis.}
}
\details{
This function implements a maximum likelihood version of the morphological model described by Lewis (2001), including the correction for only including variable characters. I've also included an option to allow for variation across sites based on a discrete gamma distribution. At the moment the function only estimates branch lengths (and the shape parameter if include.gamma==TRUE). In future versions I will add in a topology estimation algorithm. In the meantime, this function can be used to test different topological hypotheses. For example, in my case, I wrote this function as a means of testing the position of a fossil taxon. 

A couple of things to note. First, the starting branch lengths are based on the Rogers-Swofford approach used by PAUP and described in Rogers and Swofford (1998). However, rather than using their maximum parsimony reconstruction approach, I rely on the acctran parsimony implementation in \code{phangorn}. Thus, I named the approach \dQuote{Rogers-Swoffordish}. Second, I implicitly assume that the input tree is unrooted, but this need not be the case. In the end, however, the outputted tree will be unrooted because \code{phangorn} automatically unroots a rooted tree when doing the acctran parsimony analysis. Finally, I've included a nexus reader that accompanies this function. Essentially, it reads in a nexus file of character alignment and organizes the data in such a way that it is readable by \code{corHMM}. For more information about this nexus reader see \dQuote{readNexusMorph}. But in general, the input format follows \code{rayDISC}, where ambiguous characters are given as \dQuote{?}, and partial ambiguous characters are concatenated with \dQuote{&} separating them. However, I assume that the taxon names are provided as row names, not as a separate column.
}
\value{
\code{lewisMkv} returns an object of class \code{lewis.mkv}. This is a list with elements:
\item{$loglik}{The maximum negative log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{gamma.shape}{The MLE of the discrete gamma shape parameter.}
\item{$phy}{The resulting phylogeny with inferred branch lengths.}
\item{$data}{User-supplied character set.}
\item{$opts}{Internal settings of the likelihood search}
}
\examples{
## Not run
# morph.data <- readNexusMorph("AsteralesPollen_MPfixed.nex")
# phy <- read.nexus("Asterales.tre")
# tree.set <- corHMM:::AddTaxonEverywhere(phy, "Tubulifloridites_lilliei")
# pp <- lewisMkv(tree.set[[1]], morph.data, exclude.sites=c(17,18), ncores=4)
}
\references{
Felsenstein, J. 1992. Phylogenies from restriction sites: A maximum-likelihood approach. Evolution 46:159-173.

Lewis, P.O. 2001. A likelihood approach to estimating phylogeny from discrete morphological character data. Systematic Biology 50: 913-925.

Rogers, J.S., and D.L. Swofford. 1998. A fast method for approximating maximum likelihoods of phylogenetic trees from nucleotide sequences. Systematic Biology 47: 77-89.
}
\author{Jeremy M. Beaulieu}
\keyword{models}
