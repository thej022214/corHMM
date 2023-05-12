\name{ComputeCI}
\alias{ComputeCI}
\title{Compute confidence around rate estimates}
\description{corHMM gives a single point estimate of rates, but this point estimate could be very uncertain. A traditional way to evaluate confidence is to vary each parameter until the log likelihood gets 2 units worse, while holding other parameters at their maximum likelihood estimates. That's fine and fast, but it can be misled by ridges. So, instead, we want all values that lead to a likelihood within two log likelihood units of the best. The range will be at least as wide as the univariate estimates but probably much larger. }
\usage{
ComputeCI(corhmm.object, desired.delta = 2, n.points=5000, verbose=TRUE,  
print_freq=50, ...)
}       
\arguments{
\item{corhmm.object}{The result of a corHMM search.}
\item{desired.delta}{How many log likelihood units to deviate from the optimal likelihood.}
\item{n.points}{How many points to use.}
\item{print_freq}{Output progress every print_freq steps.}
\item{verbose}{Other arguments to pass into the likelihood function.}
\item{...}{further arguments to be passed dentist.}
}
\details{
The algorithm tunes: if it is moving too far away from the desired likelihoods, it will decrease the proposal width; if it staying in areas better than the desired likelihood, it will increase the proposal width. It will also expand the proposal width for parameters where the extreme values still appear good enough to try to find out the full range for these values.

In general, the idea of this is not to give you a pleasingly narrow range of possible values -- it is to try to find the actual uncertainty, including finding any ridges that would not be seen in univariate space.
}
\value{
A dentist object containing results, the data.frame of negative log likelihoods and the parameters associated with them; acceptances, the vector of whether a proposed move was accepted each step; best_neglnL, the best value passed into the analysis; delta, the desired offset; all_ranges, a summary of the results.
}
\examples{
\donttest{
# data(primates)
# phy <- multi2di(primates[[1]])
# data <- primates[[2]]
# MK_3state <- corHMM(phy = phy, data = data, rate.cat = 1)
# confidence_results <- ComputeCI(MK_3state, desired.delta = 2, 200)
# print(confidence_results)
# plot.dentist(confidence_results)
}
}
\author{Brian O'Meara (see also the R package Dentist) & exported by James D. Boyko}
\keyword{models}