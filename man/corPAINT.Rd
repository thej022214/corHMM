\name{corPAINT}
\alias{corPAINT}
\title{Binary character evolution with tree painting}
\description{Fits multiple rate models of correlated evolution between one, two, or three binary traits to paintings on branches}
\usage{
corPAINT(phy,data, ntraits=2, rate.mat=NULL, model=c("ER","SYM","ARD"), 
node.states=c("joint", "marginal", "scaled"), p=NULL, root.p=NULL, ip=NULL, 
lb=0, ub=100,diagn=FALSE)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
\item{data}{a data matrix containing species information (see Details).}
\item{ntraits}{specifies the number of traits to included in the analysis.}
\item{rate.mat}{a user-supplied rate matrix index of parameters to be optimized.}
\item{model}{specifies the underlying model.} 
\item{node.states}{method used to calculate ancestral states at internal nodes (see Details).}
\item{p}{a vector of transition rates. Allows the user to calculate the likelihood given a specified set of parameter values to specified as fixed and calculate the likelihood.}
\item{root.p}{a vector used to fix the probabilities at the root.}
\item{ip}{initial values used for the likelihood search. Can be a single value or a vector of unique values for each parameter. The default is \code{ip=1}.}
\item{lb}{lower bound for the likelihood search. The default is \code{lb=0}.}
\item{ub}{upper bound for the likelihood search. The default is \code{ub=100}.}
\item{diagn}{logical indicating whether diagnostic tests should be performed. The default is \code{FALSE}.}
}
\details{
The function fits a model that applies different transition models between one, two or three binary characters based on the user-defined painting of branches on the tree of discrete "selective regimes". The trait file must be contructed in the following way: the first column of the trait file must contain the species labels to match to the tree, with the second, and so on, corresponding to the binary traits of interest. The last column in the trait file defines the current "selective regime" for each tip. The user can fix the root state probabilities by supplying a vector to the \code{root.p}, otherwise, the program assumes the marginal probability for the root. Also, like all other functions scoring information that is missing for a species can be incorporated in the analysis by including an NA for that particular trait. NOTE THAT ALTHOUGH THIS FUNCTION SHOULD WORK, IT IS CURRENTLY BEING DEVELOPED. SO USE AT YOUR OWN RISK. 
}
\value{
\code{corPAINT} returns an object of class \code{corPAINT}. This is a list with elements:
\item{$loglik}{the maximum negative log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{Akaike information criterion corrected for sample size.}
\item{$ntraits}{The number of traits specified.}
\item{$solution}{a matrix containing the maximum likelihood estimates of the transition rates.}
\item{$solution.se}{a matrix containing the approximate standard errors of the transition rates. The standard error is calculated as the square root of the diagonal of the inverse of the Hessian matrix.}
\item{$index.mat}{The indices of the parameters being estimated are returned. The numbers correspond to the row in the \code{eigvect} and can useful for identifying the parameters that are causing the objective function to be at a saddlepoint.}
\item{$opts}{Internal settings of the likelihood search}
\item{$data}{User-supplied dataset.}
\item{$phy}{User-supplied tree.}
\item{$states}{The likeliest states at each internal node.}
\item{$tip.states}{NULL}
\item{$iterations}{The number of iterations used by the optimization routine.}
\item{$eigval}{The eigenvalues from the decomposition of the Hessian of the likelihood function. If any \code{eigval<0} then one or more parameters were not optimized during the likelihood search}
\item{$eigvect}{The eigenvectors from the decomposition of the Hessian of the likelihood function is returned}
}
\examples{
# Not run
## Load tree and data
# data(primates.paint)
## Obtain the fit for two binary characters
# pp.null<-corDISC(primates.paint$tree,primates.paint$trait,ntraits=2,model="ER",
# node.states="marginal")
# pp.null
# pp.paint<-corPAINT(primates.paint$tree,primates.paint$trait,ntraits=2,model="ER",
# node.states="marginal")
# pp.paint
}
\author{Jeremy M. Beaulieu}
\keyword{models}
