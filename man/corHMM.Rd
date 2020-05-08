\name{corHMM}
\alias{corHMM}
\title{Hidden Rates Model}
\description{Estimates hidden rates underlying the evolution of a binary character}
\usage{
corHMM(phy, data, rate.cat, rate.mat=NULL, model = "ARD", node.states = "marginal", 
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=1, sann.its=5000, 
get.tip.states = FALSE, lewis.asc.bias = FALSE)
}
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
\item{data}{a data matrix containing species information (see Details).}
\item{rate.cat}{specifies the number of rate categories.}
\item{rate.mat}{a user-supplied index of parameters to be optimized.}
\item{model}{One of "ARD", "SYM", or "ER". ARD: all rates differ. SYM: rates between any two states do not differ. ER: all rates are equal.} 
\item{node.states}{method used to calculate ancestral states at internal nodes (see Details).}
\item{fixed.nodes}{specifies that states for nodes in the phylogeny are assumed fixed. These are supplied as node labels in the \dQuote{phylo} object.}
\item{p}{a vector of transition rates. Allows the user to calculate the likelihood given a specified set of parameter values to specified as fixed and calculate the likelihood.}
\item{root.p}{a vector used to fix the probabilities at the root, but \dQuote{yang} and \dQuote{maddfitz} can also be supplied to use the method of Yang (2006) and FitzJohn et al (2009) respectively (see details).}
\item{ip}{initial values used for the likelihood search. Can be a single value or a vector of unique values for each parameter. The default is \code{ip=1}.}
\item{nstarts}{the number of random restarts to be performed. The default is \code{nstarts=0}.}
\item{n.cores}{the number of processor cores to spread out the random restarts.}
\item{sann.its}{a numeric indicating the number of times the simulated annealing algorithm should call the objective function.}
\item{get.tip.states}{a boolean indicating whether tip reconstructions should be output. The default is \code{FALSE}.}
\item{lewis.asc.bias}{a boolean indicating whether to correct for observing a dataset that is not univariate. The default is \code{FALSE}}.
}
\details{
The function takes a tree and a trait file and estimates transition rates and ancestral states for any number of discrete characters using a Markov model. These models describe evolution as discrete transitions between observed states. If rate.class > 1, then the model is a hidden Markov model (HMM; also known as a hidden rates model (HRM)). The HRM is a generalization of the covarion model that allows different rate classes to be treated as "hidden" states. Essentially a hidden Markov model allows for multiple processes to describe the evolution of your observed character. This could be another (hidden) state or a large group of them. Regardless of the reason, an HMM is saying that not all observed characters are expected to act the same way. 

From a more technical standpoint, let's say we have two rate classes, slow (S) and fast (F). They underlie each observed state of 0 and 1. Since we only observe states, we treat each observation as having a probability of 1 for being either in the F and S categories. In other words, a character state 0 at a tip is assumed to have a probability of 1 for being 0_S and 0_F. The likelihood function is then maximized using the bounded subplex optimization routine implemented in the R package \code{nloptr}, which provides a common interface to NLopt, an open-source library for nonlinear optimization. In the former case, however, it is recommended that \code{nstarts} is set to a large value (e.g. 100) to ensure that the maximum likelihood solution is found. Users can set \code{n.cores} to parse the random restarts onto multiple processors.

If more than a single column of data is given double transitions between characters are disallowed. For example, if character 1 is the presence or absence of limbs, and character 2 is the presence or absence of fingers, then the transition from absence of limbs and fingers to presence of limbs and fingers is automatically disallowed. This is consistent with Pagel's (1994) model of correlated character evolution. Additionally, the input phylogeny need not be bifurcating as the algorithm is implemented to handle multifucations. Polytomies are allowed by generalizing Felsenstein's (1981) pruning algorithm to be the product of the probability of observing the tip states of n descendant nodes, rather than two, as in the completely bifurcating case. The first column of the trait file must contain the species labels to match to the tree, with any number of columns following. Note that for a given full model, the different rate classes are ordered from slowest (rate class R1) to fastest (rate class Rn) with respect to state 0.

The user can fix the root state probabilities by supplying a vector to \code{root.p}. For example, if the hypothesis is that the root is 0_S in a model with two hidden rates, then the root vector would be \code{root.p=c(1,0,0,0)} for state combinations 0_S, 1_S, 0_F, and 1_F, respectively. If the user supplies the flag \code{root.p}=\dQuote{yang}, then the estimated transition rates are used to set the weights at the root (see pg. 124 Yang 2006), whereas specifying \code{root.p}=\dQuote{maddfitz} employs the same procedure described by Maddison et al. (2007) and FitzJohn et al. (2009). Note that the default \code{root.p=NULL} assumes equal weighting among all possible states.

Ancestral states can be estimated using marginal, joint, scaled, or none approaches. Marginal gives the likelihood of state at each node, integrating over the states at other nodes. Joint gives the optimal state at each node for the entire tree at once. Scaled is included for compatibility with ape's ace() function. None suppresses calculation of ancestral states, which can dramatically speed up calculations if you're comparing models but make plotting difficult.
}
\value{
\code{corHMM} returns an object of class \code{corHMM}. This is a list with elements:
\item{$loglik}{the maximum negative log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{Akaike information criterion corrected for sample size.}
\item{$rate.cat}{The number of rate categories specified.}
\item{$solution}{a matrix containing the maximum likelihood estimates of the transition rates. Note that the rate classes are ordered from slowest (R1) to fastest (Rn) with respect to state 0.}
\item{$index.mat}{The indices of the parameters being estimated are returned}
\item{$opts}{Internal settings of the likelihood search}
\item{$data}{User-supplied dataset.}
\item{$data.legend}{User-supplied dataset with an extra column of trait values corresponding to how corHMM calls the user data.}
\item{$phy}{User-supplied tree.}
\item{$states}{The likeliest states at each internal node. The state and rates reconstructed at internal nodes are in the order of the column headings of the rates matrix.}
\item{$tip.states}{The likeliest state at each tip. The state and rates reconstructed at the tips are in the order of the column headings of the rates matrix.}
\item{info.states}{a vector containing the amount of information (in bits) that the tip states and model gives to each node.}
\item{$iterations}{The number of iterations used by the optimization routine.}
\item{$root.p}{The root prior used in model estimation.}
}
\examples{
## Not run
# data(primates)
# phy <- multi2di(primates[[1]])
# data <- primates[[2]]
## Obtain the fit of second rate class underlying a binary character:
# MK_3state <- corHMM(phy = phy, data = data, rate.cat = 1)
# MK_3state
}
\references{
Beaulieu J.M., B.C. O'Meara, and M.J. Donoghue. 2013. Identifying hidden rate changes in the evolution of a binary morphological character: the evolution of plant habit in campanulid angiosperms. Systematic Biology 62:725-737.

Felsenstein, J. 1981. A likelihood approach to character weighting and what it tells us about parsimony and compatibility. Biological Journal of the Linnean Society 16: 183-196.

Felsenstein J. 2004. Inferring phylogenies. Sunderland MA: Sinauer Associates.

FitzJohn, R.G., W.P. Maddison, and S.P. Otto. 2009. Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Systematic Biology 58:595-611.

Maddison, W.P., P.E. Midford, and S.P. Otto. 2007. Estimating a binary characters effect on speciation and extinction. Systematic Biology 56:701-710.

Pagel, M. 1994. Detecting correlated evolution on phylogenies: a gneeral method for the comparative analysis of discrete characters. Proc. R. Soc. Lond. B 255:37-45.

Yang, Z. 2006. Computational Molecular Evolution. Oxford Press:London.
}
\author{Jeremy M. Beaulieu}
\keyword{models}
