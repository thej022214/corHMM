\name{corHMM}
\alias{corHMM}
\title{Hidden Rates Model}
\description{Estimates hidden rates underlying the evolution of a binary character}
\usage{
corHMM(phy, data, rate.cat, rate.mat=NULL, node.states=c("joint", "marginal","scaled"), 
optim.method=c("subplex"), p=NULL, root.p=NULL, ip=NULL, nstarts=10, n.cores=NULL, 
lb=0, ub=100, diagn=FALSE)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
\item{data}{a data matrix containing species information (see Details).}
\item{rate.cat}{specifies the number of rate categories in the HRM.} 
\item{rate.mat}{a user-supplied rate matrix index of parameters to be optimized.}
\item{node.states}{method used to calculate ancestral states at internal nodes (see Details).}
\item{optim.method}{method used to perform optimization. The default is \code{subplex}.} 
\item{p}{a vector of transition rates. Allows the user to calculate the likelihood given a specified set of parameter values to specified as fixed and calculate the likelihood.}
\item{root.p}{a vector used to fix the probabilities at the root, but \dQuote{maddfitz} can also be supplied to use the method of Maddison et al (2007) and FitzJohn et al (2009) (see details).}
\item{ip}{initial values used for the likelihood search. Can be a single value or a vector of unique values for each parameter. The default is \code{ip=1}.}
\item{nstarts}{the number of random restarts to be performed. The default is \code{nstarts=10}.}
\item{n.cores}{the number of processor cores to spread out the random restarts.}
\item{lb}{lower bound for the likelihood search. The default is \code{lb=0}.}
\item{ub}{upper bound for the likelihood search. The default is \code{ub=100}.}
\item{diagn}{logical indicating whether diagnostic tests should be performed. The default is \code{FALSE}.}
}
\details{
The function takes a tree and a trait file and estimates transition rates and ancestral states for a single binary character using the hidden rates model (HRM). The HRM is a generalization of the covarion model that allows different rate classes to be treated as "hidden" states in reconstructing ancestral character states. For example, for a model with two rate classes, slow (S) and fast (F), underlie each observed state of 0 and 1. Since we only observe states, we treat each observation as having a probability of 1 for being either in the F and S categories. In other words, a character state 0 at a tip is assumed to have a probability of 1 for being 0_S and 0_F. The likelihood function is then maximized using the bounded subplex optimization routine (\code{optim.method=subplex}) implemented in the R package \code{nloptr}, which provides a common interface to NLopt, an open-source library for nonlinear optimization. Users can also set \code{optim.method=rgenoud} to specify that likelihood function is to maximized using a genetic algorithm. In the former case, however, it is recommended that \code{nstarts} is set to a large value (e.g. 100) to ensure that the maximum likelihood solution is found. Users can set \code{n.cores} to parse the random restarts onto multiple processors.

The input phylogeny need not be bifurcating as the algorithm is implemented to handle multifucations. Polytomies are allowed by generalizing Felsenstein's (1981) pruning algorithm to be the product of the probability of observing the tip states of n descendant nodes, rather than two, as in the completely bifurcating case. The first column of the trait file must contain the species labels to match to the tree, with the second corresponding to the binary trait of interest. Any variant of a model that assume either 1, 2, 3, 4, or 5 rate categories underlying the observed data can be evaluated. Note that for a given full model, the different rate classes are ordered from slowest (rate class R1) to fastest (rate class Rn) with respect to state 0. The user can fix the root state probabilities by supplying a vector to \code{root.p}. For example, if the hypothesis is that the root is 0_S in a model with two hidden rates, then the root vector would be \code{root.p=c(1,0,0,0)} for state combinations 0_S, 1_S, 0_F, and 1_F, respectively. If the user supplies the flag \code{root.p}=\dQuote{maddfitz} the same procedure described by Maddison et al. (2007) and FitzJohn et al. (2009) is used. Note that the default \code{root.p=NULL} assumes equal weighting among all possible states. 
}
\value{
\code{corHMM} returns an object of class \code{corHMM}. This is a list with elements:
\item{$loglik}{the maximum negative log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{Akaike information criterion corrected for sample size.}
\item{$rate.cat}{The number of rate categories specified.}
\item{$solution}{a matrix containing the maximum likelihood estimates of the transition rates. Note that the rate classes are ordered from slowest (R1) to fastest (Rn) with respect to state 0}
\item{$solution.se}{a matrix containing the approximate standard errors of the transition rates. The standard error is calculated as the square root of the diagonal of the inverse of the Hessian matrix.}
\item{$index.mat}{The indices of the parameters being estimated are returned. The numbers correspond to the row in the \code{eigvect} and can useful for identifying the parameters that are causing the objective function to be at a saddlepoint.}
\item{$opts}{Internal settings of the likelihood search}
\item{$data}{User-supplied dataset.}
\item{$phy}{User-supplied tree.}
\item{$states}{The likeliest states at each internal node. The state and rates reconstructed at internal nodes are in the order of the column headings of the rates matrix.}
\item{$tip.states}{NULL}
\item{$iterations}{The number of iterations used by the optimization routine.}
\item{$eigval}{The eigenvalues from the decomposition of the Hessian of the likelihood function. If any \code{eigval<0} then one or more parameters were not optimized during the likelihood search}
\item{$eigvect}{The eigenvectors from the decomposition of the Hessian of the likelihood function is returned}
}
\examples{
## Not run
# data(primates)
## Obtain the fit of second rate class underlying a binary character:
# pp<-corHMM(primates$tree,primates$trait[,c(1,2)],rate.cat=2,node.states="marginal")
# pp
}
\references{
Beaulieu J.M., B.C. O'Meara, and M.J. Donoghue. 2013. Identifying hidden rate changes in the evolution of a binary morphological character: the evolution of plant habit in campanulid angiosperms. Systematic Biology In press.

Felsenstein, J. 1981. A likelihood approach to character weighting and what it tells us about parsimony and compatibility. Biological Journal of the Linnean Society 16: 183-196.

Felsenstein J. 2004. Inferring phylogenies. Sunderland MA: Sinauer Associates. 

FitzJohn, R.G., W.P. Maddison, and S.P. Otto. 2009. Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Systematic Biology 58:595-611.

Maddison, W.P., P.E. Midford, and S.P. Otto. 2007. Estimating a binary characters effect on speciation and extinction. Systematic Biology 56:701-710.
}
\author{Jeremy M. Beaulieu}
\keyword{models}
