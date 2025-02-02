\name{corHMMDredgeBase}
\alias{corHMMDredgeBase}
\title{Hidden Rates Model with regularization with Penalization Options}
\description{
\code{corHMMDredgeBase} fits a hidden Markov model (HMM) to a given phylogenetic tree and character data. It offers additional options for penalization and optimization compared to the standard \code{corHMM} function.
}
\usage{
corHMMDredgeBase(phy, data, rate.cat, root.p="maddfitz", pen.type = "l1", lambda = 1, 
rate.mat=NULL, node.states = "marginal", fixed.nodes=FALSE, ip=NULL, nstarts=0, n.cores=1, 
get.tip.states = FALSE,lewis.asc.bias = FALSE, collapse = FALSE, lower.bound = 1e-10, 
upper.bound = 100, opts=NULL, p=NULL, grad=FALSE)
}
\arguments{
  \item{phy}{
    An object of class \code{phylo}, representing the phylogenetic tree.
  }
  \item{data}{
    A data frame containing character states for the tips of the phylogeny. The first column should match the tip labels in \code{phy}, and the second onwards column should contain the observed states.
  }
  \item{rate.cat}{
    An integer specifying the number of rate categories in the model.
  }
  \item{root.p}{
    A vector of probabilities for the root states or a method to estimate them. The default is the \code{"maddfitz"} method.
  }
  \item{pen.type}{
    The type of penalization applied to the model. Options include \code{"l1"}, \code{"l2"}, \code{"er"}, and \code{"unreg"}. See Details.
  }
  \item{lambda}{
    A hyper-parameter that adjusts the severity of the penalty, ranging from 0 (no regularization) to 1 (full penalization). Default is 1.
  }
  \item{rate.mat}{
    A user-supplied matrix  containing indexes of parameters to be optimized. If NULL, an all rates different model is estimated.
  }
  \item{node.states}{
    A method for estimating node states. Options include \code{"marginal"}, \code{"joint"}, and \code{"none"}.
  }
  \item{fixed.nodes}{
    Specifies that states for nodes in the phylogeny are assumed fixed. These are supplied as node labels in the "phylo" object. Default is \code{FALSE}.
  }
  \item{ip}{
    Initial values used for the likelihood search. Can be a single value or a vector of unique values for each parameter. The default is ip=1.
  }
  \item{nstarts}{
    The number of random restarts to be performed. The default is nstarts=0.
  }
  \item{n.cores}{
    The number of processor cores to spread out the random restarts.
  }
  \item{get.tip.states}{
    Logical value indicating whether tip reconstructions should be output. The default is FALSE.
  }
  \item{lewis.asc.bias}{
  Logical value indicating whether to correct for observing a dataset that is not univariate. The default is FALSE
}
  \item{collapse}{
    A logical value indicating whether to collapse branches with no variation in states. Default is \code{FALSE}.
  }
  \item{lower.bound}{
    The lower bound for the rate parameters during optimization. Default is \code{1e-10}.
  }
  \item{upper.bound}{
    The upper bound for the rate parameters during optimization. Default is \code{100}.
  }
  \item{opts}{
    A list of options to pass to nloptr.
  }
  \item{p}{
    A vector of transition rates. Allows the user to calculate the likelihood given a specified set of parameter values to specified as fixed and calculate the likelihood.
  }
  \item{grad}{
    A logical value indicating whether to use gradient-based optimization. Default is \code{FALSE}.
  }
}
\details{
There are 3 penalty types available for users (though this may be expanded in the future). They are \code{l1}, \code{l2}, and \code{er}. The first two penalty types are analagous to lasso and ridge regression. Whereas the \code{er} penalization is analagous to Zhou et al.'s (2024) McLasso and penalizes the distance between rate values (i.e., it prefers rate matrices that are closer to an equal-rates model).

Under an \code{l1} regularization scheme, the likelihood is penalized by the mean transition rate. The mean is used instead of the sum because a sum would overly penalize more complex models by virtue of them having more parameters. This leads to scenarios where simpler models have much better likelihoods than more complex models

Under an \code{l2} regularization scheme, the likleihood is penalized by the squared mean transition rate. 

Under an \code{er} regularization scheme, the likleihood is penalized by the average distance between parameters. If all the parameters were the same (an equal rates model), the penalty would be 0. 

A user can also set \code{pen.type = 'unreg'} for an unpenalized likelihood which would be identical to the original implementation of \code{corHMM}. More work needs to be done to determine when to use each type of penalization, but generally using any penalization will be good for smaller datasets which tend to be high variance. \code{l2} is the most aggresive penalization, shrinking paramaters more quickly than other methods and leading to more dropped (and potentially finding more unecessary) parameters. \code{er} is the most similar to an unregularized model as it does not necessarily penalize high parameter values. It will however penalize a model that has one parameter that is much higher than the rest unless there is significant evidence that this outlier parameter is needed. In practice, \code{er}, behaves very similarly to unregularized models. \code{l1} regularization is an intermediate penalization between \code{l2} and \code{er}.

The \code{grad} option employs a numerical gradient for the optimization. This is a particularly inefficient way to find a gradient as it will require at least \code{k} iterations per likelihood calculation. However, I have found that this improves the search stability and speed as the number of iterations is greatly reduced when a gradient is provided. This is also important in cases where there are a lot of parameters (\code{k} is large). In these cases the parameter space is so highly dimensional that many search algorithms struggle. In the future I hope to implement a more efficient gradient calculation and combine a gradient based local optimizaiton with a global optimization scheme. 

*Note: Many details of \code{corHMMDredgeBase} and \code{corHMM} are the same and will not be repeated here. If an arguement is unclear check the Details section of \code{corHMM}. The focus of these details will be on the unique aspects of the dredge fitting approach.
}
\value{
\code{corHMM} returns an object of class \code{corHMM}. This is a list with elements:
\item{$loglik}{the maximum negative log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{Akaike information criterion corrected for sample size.}
\item{$rate.cat}{The number of rate categories specified.}
\item{$solution}{a matrix containing the maximum likelihood estimates of the transition rates. Note that the rate classes are ordered from slowest (R1) to fastest (Rn) with respect to state 0.}
\item{$index.mat}{The indices of the parameters being estimated are returned. This also is a way to allow the estimation of transition rates for parameters not oberved in the dataset. Say you have 2 traits X and Y, where the combinations 00, 01, and 11 are observed (10 is not). A 4 by 4 index matrix could be used to force 10 into the model.}
\item{$data}{User-supplied dataset.}
\item{$data.legend}{User-supplied dataset with an extra column of trait values corresponding to how corHMM calls the user data.}
\item{$phy}{User-supplied tree.}
\item{$states}{The likeliest states at each internal node. The state and rates reconstructed at internal nodes are in the order of the column headings of the rates matrix.}
\item{$tip.states}{The likeliest state at each tip. The state and rates reconstructed at the tips are in the order of the column headings of the rates matrix.}
\item{$states.info}{a vector containing the amount of information (in bits) that the tip states and model gives to each node.}
\item{$iterations}{The number of iterations used by the optimization routine.}
\item{$root.p}{The root prior used in model estimation.}
}
\examples{
\donttest{
data(primates)
phy <- multi2di(primates[[1]])
data <- primates[[2]]
MK_3state <- corHMMDredgeBase(phy = phy, data = data, rate.cat = 1, pen.type = "l1", 
	root.p = "maddfitz", lambda = 1)
MK_3state
}
}
\references{
Boyko, J. D. 2024. Automatic Discovery of Optimal Discrete Character Models through Regularization. In prep.

Zhou, Y., Gao, M., Chen, Y., Shi, X., 2024. Adaptive Penalized Likelihood method for Markov Chains. 
}
\author{James D. Boyko}
\keyword{models}
