\name{get_batch_profile_lik}
\alias{get_batch_profile_lik}
\title{Perform Batch Profile Likelihood Analysis for Multiple Parameters}
\description{
\code{get_batch_profile_lik} is a wrapper function that performs profile likelihood analysis for multiple parameters of a fitted \code{corHMM} object. It allows exploration of the likelihood surface by evaluating the likelihood at various points along the parameter space.
}
\usage{
get_batch_profile_lik(corhmm_obj, range_factor, n_points, verbose=FALSE, ncores = NULL, 
dredge = FALSE)
}
\arguments{
  \item{corhmm_obj}{
    A fitted \code{corHMM} object representing the model for which profile likelihood analysis is to be performed.
  }
  \item{range_factor}{
    A numeric factor determining the range over which to generate points for the profile likelihood analysis. This value is used to calculate the bounds around the maximum likelihood estimates (MLEs) for each parameter.
  }
  \item{n_points}{
    The number of points to generate for each parameter along its profile likelihood curve.
  }
    \item{verbose}{
    Whether to print messages about which parameter is being optimized. 
  }
  \item{ncores}{
    The number of processor cores to be used for parallel computation. Default is \code{NULL}, which uses a single core.
  }
  \item{dredge}{
    Logical value indicating whether to include model penalization factors such as \code{pen_type} and \code{lambda} from the fitted \code{corHMM} object. Default is \code{FALSE}.
  }
}
\details{
This function performs a profile likelihood analysis for each parameter in a fitted \code{corHMM} model. It evaluates the likelihood at logarithmically spaced points around the maximum likelihood estimates (MLEs) for the parameters. If \code{dredge} is set to \code{TRUE}, the function also considers model penalization terms.

The function works by first generating a set of points along a logarithmic scale for each parameter, then fixing one parameter at a time while optimizing over the others. The resulting profile likelihood values are returned for each parameter.

Parallel computation can be enabled using the \code{ncores} argument to speed up the analysis.
}
\value{
A list containing the profile likelihood results for each parameter. Each entry in the list corresponds to a parameter and contains the profile likelihood values across the range of points evaluated. The original \code{corHMM} object is also included in the output.
}
\examples{
\donttest{
# Assuming you have a fitted corHMM object:
data(primates)
phy <- multi2di(primates[[1]])
data <- primates[[2]]
corhmm_fit <- corHMM(phy = phy, data = data, rate.cat = 1)
profile_results <- get_batch_profile_lik(corhmm_fit, range_factor = 2, n_points = 50)
plot_batch_profile_lik(profile_results)
}
}
\seealso{
\code{\link{corHMM}} for fitting hidden Markov models to phylogenetic data.
}
\author{
Your Name
}
