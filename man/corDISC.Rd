\name{corDISC}
\alias{corDISC}
\title{Correlated evolution binary traits}
\description{Fits a model of correlated evolution between two or three binary traits}
\usage{
corDISC(phy,data, ntraits=2, rate.mat=NULL, model=c("ER","SYM","ARD"), 
node.states=c("joint", "marginal", "scaled"), p=NULL, root.p=NULL, ip=NULL, 
lb=0, ub=100, diagn=FALSE)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
\item{data}{a data matrix containing species information (see Details).}
\item{ntraits}{specifies the number of traits to be included in the analysis.}
\item{rate.mat}{a user-supplied rate matrix index of parameters to be optimized.}
\item{model}{specifies the underlying model.} 
\item{node.states}{method used to calculate ancestral states at internal nodes (see Details).}
\item{p}{a vector of transition rates. Allows the user to calculate the likelihood given a specified set of parameter values to specified as fixed and calculate the likelihood.}
\item{root.p}{a vector used to fix the probabilities at the root, but \dQuote{maddfitz} can also be supplied to use the method of Maddison et al (2007) and FitzJohn et al (2009) (see details).}
\item{ip}{initial values used for the likelihood search. Can be a single value or a vector of unique values for each parameter. The default is \code{ip=1}.}
\item{lb}{lower bound for the likelihood search. The default is \code{lb=0}.}
\item{ub}{upper bound for the likelihood search. The default is \code{ub=100}.}
\item{diagn}{logical indicating whether diagnostic tests should be performed. The default is \code{FALSE}.}
}
\details{
The function takes a tree and a trait file and estimates transition rates and ancestral states for two or three binary characters (see Pagel 1994). Note, however, that rayDISC can be used to evaluate the same models as in corDISC, with the major difference being that, with rayDISC, the rate matrix would have to be manipulated using \code{rate.mat.maker} in order to remove parameters associated with dual transitions. With corDISC, the input phylogeny need not be bifurcating as the algorithm is implemented to handle multifucations. Polytomies are allowed by generalizing Felsenstein's (1981) pruning algorithm to be the product of the probability of observing the tip states of n descendant nodes, rather than two, as in the completely bifurcating case. For the trait file, the first column of the trait file must contain the species labels to match to the tree, with the second column onwards corresponding to the binary traits of interest. The user can fix the root state probabilities by supplying a vector to \code{root.p}. For example, in the two trait case, if the hypothesis is that the root is 00, then the root vector would be \code{root.p=c(1,0,0,0)} for state combinations 00, 01, 10, and 11, respectively. If user supplies the flag \code{root.p}=\dQuote{maddfitz} the same procedure described by Maddison et al (2007) and FitzJohn et al (2009) is used. Note, however, that the default \code{root.p=NULL} assumes equal weighting among all possible states.

We also note that scoring information that is missing for a species can be incorporated in the analysis by including an NA for that particular trait. corDISC will then set the trait vector so that the tip vector will reflect the probabilities that are compatible with our observations. For example, if the scoring for trait 1 is missing, but trait 2 is scored as 0, then the tip vector would be (1,0,1,0), for state combinations 00, 01, 10, and 11 respectively, given our observation that trait 2 is scored 0 (for a good discussion see Felsenstein 2004, pg. 255).
}
\value{
\code{corDISC} returns an object of class \code{corDISC}. This is a list with elements:
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
## Not run
## Load tree and data
# data(primates)

## Obtain the fit for two binary characters
# pp<-corDISC(primates$tree,primates$trait,ntraits=2,model="ARD",
# node.states="marginal", diagn=FALSE)
# pp

## State combination three is not an observed state, so for fun, let's remove
## these transitions:
# new.mat <- rate.mat.maker(hrm=FALSE, ntraits=2, model="ARD")
# new.mat <- rate.par.drop(new.mat, c(2,8,5,6))
# pp<-corDISC(primates$tree,primates$trait,ntraits=2,rate.mat=new.mat,model="ARD",
# node.states="marginal", diagn=FALSE)
# pp
}
\references{
Beaulieu J.M., and M.J. Donoghue 2013. Fruit evolution and diversification in campanulid angiosperms. Evolution In press.

Felsenstein, J. 1981. A likelihood approach to character weighting and what it tells us about parsimony and compatibility. Biological Journal of the Linnean Society 16: 183-196.

Felsenstein J. 2004. Inferring phylogenies. Sunderland MA: Sinauer Associates.

FitzJohn, R.G., W.P. Maddison, and S.P. Otto. 2009. Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Systematic Biology 58:595-611.

Maddison, W.P., P.E. Midford, and S.P. Otto. 2007. Estimating a binary characters effect on speciation and extinction. Systematic Biology 56:701-710.

Pagel, M. 1994. Detecting correlated evolution on phylogenies: a general method for the comparative analysis of discrete characters. Proceedings of the Royal Society, B. 255:37-45.
}
\author{Jeremy M. Beaulieu}
\keyword{models}
