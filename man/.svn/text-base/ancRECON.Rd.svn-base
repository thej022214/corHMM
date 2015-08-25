\name{ancRECON}
\alias{ancRECON}
\title{Ancestral state reconstruction}
\description{Infers ancestral states based on a set of model parameters}
\usage{
ancRECON(phy,data, p, method=c("joint", "marginal", "scaled"), hrm=FALSE, 
rate.cat, ntraits=NULL, charnum=NULL, rate.mat=NULL,
model=c("ER", "SYM", "ARD"), root.p=NULL)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
\item{data}{a data matrix containing species information (see Details).}
\item{p}{a vector of transition rates to be used to estimate ancestral states.}
\item{method}{method used to calculate ancestral states at internal nodes (see Details).}
\item{hrm}{a logical indicating whether the underlying model is the hidden rates model (HRM). The default is \code{FALSE}.} 
\item{rate.cat}{specifies the number of rate categories in the HRM.} 
\item{ntraits}{specifies the number of traits in the data file if the underlying model is not the HRM.}
\item{charnum}{specifies the number of characters in the data file used in rayDISC.}
\item{rate.mat}{a user-supplied rate matrix index of parameters to be optimized.}
\item{model}{if the model is not HRM, specifies the underlying model.}
\item{root.p}{a vector used to fix the probabilities at the root, but \dQuote{maddfitz} can also be supplied to use the method of Maddison et al (2007) and FitzJohn et al (2009) (see details).}

}
\details{
This is a stand alone function for computing the marginal, joint, or scaled likelihoods of internal nodes for a given set of transition rates. Like all other functions contained in corHMM, the tree does not have to be bifurcating in order for analyses to be carried out. IMPORTANT: If the corDISC, corHMM, and rayDISC functions are used they automatically provide a tree with the likeliest states as internal node labels. This function is intended for circumstances where the user would like to reconstruct states based on rates estimated elsewhere (e.g. BayesTraits, Mesquite, \code{ape}).

The algorithm based on Pupko et al. (2000, 2002) is used to calculate the \code{joint} estimates of ancestral states. The \code{marginal} method was originally based on a description of an algorithm by Yang (2006). The basic idea of Yang (2006) is that the tree is rerooted on each internal node, with the marginal likelihood being the probabilities of observing the tips states given that the focal node is the root. However, this takes a ton of time as the number of nodes increase. But more importantly, this does not work easily when the model contains asymmetric rates. Here we use the same dynamic programming algorithm as Mesquite (Maddison and Maddison, 2011) and \code{diversitree} (FitzJohn et al. 2009), which is time linear with the number of species and calculates the marginal probability at a node using an additional up and down pass of the tree.  If \code{scaled}, the function uses the same algorithm from ace(). Note that the \code{scaled} method of ace() is simply the conditional likelihoods of observing everything at or above the focal node and these should generally not be used for ancestral state estimation.

The user can fix the root state probabilities by supplying a vector to \code{root.p}. For example, in the two trait case, if the hypothesis is that the root is 00, then the root vector would be \code{root.p=c(1,0,0,0)} for state combinations 00, 01, 10, and 11, respectively. If the user supplies the flag \code{root.p}=\dQuote{maddfitz} the same procedure described by Maddison et al (2007) and FitzJohn et al (2009) is used. Note, however, that the default \code{root.p=NULL} assumes equal weighting among all possible states. 
}
\value{
For the \code{joint}, a vector of likeliest states at internal nodes and tips. For either \code{marginal} or \code{scaled}, a matrix of the probabilities of each state for each internal node are returned. 
}
\examples{
# Not run
## Load tree and trait
# data(primates)
## Obtain the marginal reconstruction for a set of parameters:
# param<-c(0.05,10,0.01,0.01,0.06,0,0.02,51.2)
# states<-ancRECON(primates$tree,primates$trait,p=param,method="marginal",
# hrm=FALSE,ntraits=2,model="ARD")
## Put likeliest states on the tree:
# pr<-apply(states$lik.anc.states,1,which.max)
# primates$tree$node.label <- pr
}
\references{
FitzJohn, R.G., W.P. Maddison, and S.P. Otto. 2009. Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Systematic Biology 58:595-611.

Maddison, W.P. and D.R. Maddison. 2011. Mesquite: a modular system for evolutionary analysis. Version 2.75 http://mesquiteproject.org

Pupko, T., I. Pe'er, R. Shamir, and D. Graur. 2000. A fast algorithm for joint reconstruction of ancestral amino-acid sequences. Molecular Biology and Evolution 17:890-896.

Pupko, T., I. Pe'er, D. Graur, M. Hasegawa, and N Friedman N. 2002. A branch-and-bound algorithm for the inference of ancestral amino-acid sequences when the replacement rate varies among sites: application to the evolution of five gene families. Bioinformatics 18:1116-1123.

Yang, Z. 2006. Computational Molecular Evolution. London:Oxford.
}

\author{Jeremy M. Beaulieu and Jeffrey C. Oliver}
\keyword{reconstructions}

