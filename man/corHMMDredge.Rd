\name{corHMMDredge}
\alias{corHMMDredge}
\title{Automatic search for optimal discrete character model with regularization with Penalization Options}
\description{
\code{corHMMDredge} fits a series of hidden Markov models (HMMs) to a given phylogenetic tree and discrete character data to automatically find the optimal model structure. It offers additional options for penalization and optimization compared to the standard \code{corHMM} function.
}
\usage{
corHMMDredge(phy, data, max.rate.cat, root.p="maddfitz", pen.type = "l1", lambda = 1, 
drop.par = TRUE, drop.threshold = 1e-7,  info.threshold=2, criterion="AIC", 
merge.params=TRUE, merge.threshold=0, rate.mat=NULL, 
node.states = "marginal", fixed.nodes=FALSE, ip=NULL, nstarts=0, n.cores=1, 
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = FALSE, lower.bound = 1e-10, 
upper.bound = 100, opts=NULL, verbose=TRUE, p=NULL, rate.cat=NULL, grad=FALSE)
}
\arguments{
  \item{phy}{
    An object of class \code{phylo}, representing the phylogenetic tree.
  }
  \item{data}{
    A data frame containing character states for the tips of the phylogeny. The first column should match the tip labels in \code{phy}, and the second onwards column should contain the observed states.
  }
  \item{max.rate.cat}{
    The maximum number of rate categories to try.
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
  \item{drop.par}{
    Logical. Whether to drop parameters during optimization based on a threshold. Default is \code{TRUE}.
  }
  \item{drop.threshold}{
    A numeric value determining the threshold below which parameters should be dropped. Default is \code{1e-7}.
  }
  \item{info.threshold}{
    A numeric value specifying the threshold for the amount of information required for parameter estimation. Default is \code{2}.
  }
  \item{criterion}{
    The model selection criterion to use. Options are \code{"AIC"}, \code{"AICc"}, or \code{"BIC"}. Default is \code{"AIC"}.
  }
  \item{merge.params}{
    Logical. Whether to merge similar parameters during the model search. Default is \code{TRUE}.
  }
  \item{merge.threshold}{
    A numeric threshold to determine when parameters should be merged. Default is \code{0}.
  }
  \item{rate.mat}{
    A user-supplied matrix containing indexes of parameters to be optimized. If NULL, an all-rates-different model is estimated.
  }
  \item{node.states}{
    A method for estimating node states. Options include \code{"marginal"}, \code{"joint"}, and \code{"none"}.
  }
  \item{fixed.nodes}{
    Logical. Specifies whether the states for nodes in the phylogeny are assumed fixed. These are supplied as node labels in the \code{phylo} object. Default is \code{FALSE}.
  }
  \item{ip}{
    Initial values used for the likelihood search. Can be a single value or a vector of unique values for each parameter. The default is \code{ip = 1}.
  }
  \item{nstarts}{
    The number of random restarts to be performed. Default is \code{nstarts = 0}.
  }
  \item{n.cores}{
    The number of processor cores to spread out the random restarts. Default is \code{n.cores = 1}.
  }
  \item{get.tip.states}{
    Logical. Indicates whether tip reconstructions should be output. Default is \code{FALSE}.
  }
  \item{lewis.asc.bias}{
    Logical. Indicates whether to correct for observing a dataset that is not univariate. Default is \code{FALSE}.
  }
  \item{collapse}{
    Logical. Indicates whether to collapse branches with no variation in states. Default is \code{FALSE}.
  }
  \item{lower.bound}{
    The lower bound for the rate parameters during optimization. Default is \code{1e-10}.
  }
  \item{upper.bound}{
    The upper bound for the rate parameters during optimization. Default is \code{100}.
  }
  \item{opts}{
    A list of options to pass to \code{nloptr} for controlling optimization behavior.
  }
  \item{verbose}{
    Logical. If \code{TRUE}, detailed messages about the model fitting process will be printed. Default is \code{TRUE}.
  }
  \item{p}{
    A vector of transition rates. Allows the user to calculate the likelihood given a specified set of parameter values to be fixed and calculate the likelihood.
  }
  \item{rate.cat}{
    An integer specifying the number of rate categories in the model. Only useful if fitting a fixed value \code{p}.
  }
  \item{grad}{
    Logical. If \code{TRUE}, numerical gradient-based optimization will be used. Default is \code{FALSE}. This is useful for highly parameterized models, but because it is a numerical gradient it is slow. 
  }
}
\details{
\code{corHMMDredge} will automatically search different model structures dropping parameters which are estimated near 0 and/or equating parameter values which are near one another. This can be combined with a regularization term (see below) to encrouage lower rate values and thus lead to more parameters being dropped. It will do this iteratively until a stopping criterion is met. The stopping criteria is currently a dAIC of 2, meaning if the next step has made the model worse as indicated by dAIC > 2, the dredge will stop. No model averaging is conducted and only the best model should be used from a dredge search. I explain this in more detail in Boyko (2024), but by dredging we are not specifying a model set ourselves of distinct hypotheses. Model averaging is useful in that case, but not in the dredge case, because each model as it relates to a hypothesis provides unique information about the system, but the dredge model fits can be very similar to one another, differing in only one or two parameters. These do not really provide unique information and are just minor variations of essentially the same model. 

There are 3 penalty types available for users (though this may be expanded in the future). They are \code{l1}, \code{l2}, and \code{er}. The first two penalty types are analagous to lasso and ridge regression. Whereas the \code{er} penalization is analagous to Zhou et al.'s (2024) McLasso and penalizes the distance between rate values (i.e., it prefers rate matrices that are closer to an equal-rates model).

Under an \code{l1} regularization scheme, the likelihood is penalized by the mean transition rate. The mean is used instead of the sum because a sum would overly penalize more complex models by virtue of them having more parameters. This leads to scenarios where simpler models have much better likelihoods than more complex models

Under an \code{l2} regularization scheme, the likleihood is penalized by the squared mean transition rate. 

Under an \code{er} regularization scheme, the likleihood is penalized by the average distance between parameters. If all the parameters were the same (an equal rates model), the penalty would be 0. 

A user can also set \code{pen.type = 'unreg'} for an unpenalized likelihood which would be identical to the original implementation of \code{corHMM}. More work needs to be done to determine when to use each type of penalization, but generally using any penalization will be good for smaller datasets which tend to be high variance. \code{l2} is the most aggresive penalization, shrinking paramaters more quickly than other methods and leading to more dropped (and potentially finding more unecessary) parameters. \code{er} is the most similar to an unregularized model as it does not necessarily penalize high parameter values. It will however penalize a model that has one parameter that is much higher than the rest unless there is significant evidence that this outlier parameter is needed. In practice, \code{er}, behaves very similarly to unregularized models. \code{l1} regularization is an intermediate penalization between \code{l2} and \code{er}.

The \code{grad} option employs a numerical gradient for the optimization. This is a particularly inefficient way to find a gradient as it will require at least \code{k} iterations per likelihood calculation. However, I have found that this improves the search stability and speed as the number of iterations is greatly reduced when a gradient is provided. This is also important in cases where there are a lot of parameters (\code{k} is large). In these cases the parameter space is so highly dimensional that many search algorithms struggle. In the future I hope to implement a more efficient gradient calculation and combine a gradient based local optimizaiton with a global optimization scheme. 

*Note: Many details of \code{corHMMDredgeBase} and \code{corHMM} are the same and will not be repeated here. If an arguement is unclear check the Details section of \code{corHMM}. The focus of these details will be on the unique aspects of the dredge fitting approach.
}
\value{
\code{fit_set} returns an object of class \code{corhmm.dredge}. This is a list of the models being fit. Each element of that list is of class \code{corhmm}.
}
\examples{
\donttest{
data(primates)
phy <- multi2di(primates[[1]])
data <- primates[[2]]
dredge_fits <- corHMMDredge(phy = phy, data = data, max.rate.cat = 1, pen.type = "l1", 
	root.p = "maddfitz", lambda = 1)
getModelTable(dredge_fits)
}
}
\references{
Boyko, J. D. 2024. Automatic Discovery of Optimal Discrete Character Models through Regularization. In prep.

Zhou, Y., Gao, M., Chen, Y., Shi, X., 2024. Adaptive Penalized Likelihood method for Markov Chains. 
}
\author{James D. Boyko}
\keyword{models}
