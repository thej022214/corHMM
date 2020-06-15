\name{getFullMat}
\alias{getFullMat}
\alias{getRateCatMat}
\alias{equateStateMatPars}
\alias{dropStateMatPars}
\title{Combines several rate class index matrices}
\description{Combines several index matrices which describe transitions between observed states into output a single index matrix for use in corHMM}
\usage{
getFullMat(StateMats, RateClassMat = NULL)
}       
\arguments{
\item{StateMats}{A list of index matrices describing transitions between observed states. Each unique number from 1 to n, will be independently estimated. Values of 0 are not estimated. Matrix entries of the same value are estimated to be the same rate.}
\item{RateClassMat}{An optional index matrix which describes how StateMats are related to one another. This will be a matrix of size: length(StateMats) by length(StateMats). By default, all transitions between StateMats are allowed and independently estimated.}
}
\details{
This function is the final step in creating a custom hidden Markov model. It takes a list of index matrices (StateMats) which describe different ways that the observed states are related to one another and creates a single matrix to describe the model. The matrices are combined following equation 3 of Tarasov (2019). getFullMat is part of several functions which help the user efficiently create custom index matrices. Often, it will be more practical to begin constructing a custom model with \code{getRateMat4Dat}. 

getStateMat will generate an index matrix of size n by n in which all transitions between the n states are allowed and independently estimated. That index matrix can then be manipulated by dropStateMatPars and equateStateMatPars. dropStateMatPars will drop specific rates from an index matrix. dropStateMatPars requires an index matrix and a vector of which rates should be dropped. equateStateMatPars will equates rates within an index matrix. equateStateMatPars requires an index matrix and a list of vectors each element of which should correspond to two or more rates to be equated.
}
\value{
Returns an index matrix.
}
\seealso{
\code{getRateMat4Dat}
}
\examples{
data(primates)
phy <- primates[[1]]
phy <- multi2di(phy)
data <- primates[[2]]
# create a legend and rate mat from a multi-character dataset.
LegendAndRateMat <- getStateMat4Dat(data)
rate.mat <- LegendAndRateMat$rate.mat
legend <- LegendAndRateMat$legend

# To create a hidden markov model first define your rate classes (state-dependent processes)
# R1 will be a manually created SYM model
R1 <- equateStateMatPars(rate.mat, c(1:6))
# R2 will only allow transitions between 1 and 2
R2 <- dropStateMatPars(rate.mat, c(3,4))
# R1 and R2 will transtion at equal rates (i.e. the parameter process will be ER)
P <- getRateCatMat(2)
P <- equateStateMatPars(P, c(1,2))
# combine our state-dependnet processes and parameter process
HMM <- getFullMat(list(R1, R2), P)

# This can now be used in a corHMM model
CustomModel <- corHMM(phy = phy, data = data, rate.cat = 2, rate.mat = HMM, node.states = "none")
}
\author{James D. Boyko}
\keyword{models}