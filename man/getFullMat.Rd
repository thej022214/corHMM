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
#Generate 2 matrices for 3 traits:
StateMatA <- getStateMat(3)
StateMatB <- getStateMat(3)
#Modify StateMatA such that allowed transitions are 1 to 2, 2 to 3, and 3 to 1
StateMatA_constrained <- dropStateMatPars(StateMatA, c(5, 1, 4))
#Modify StateMatB to be identical to the "SYM" model
StateMatB_constrained <- equateStateMatPars(StateMatB, list(c(1,3),c(2,5),c(4,6)))
#Make a list of our state matrices (each representing its own rate class)
StateMats <- list(StateMatA_constrained, StateMatB_constrained)
#Combine into a single index matrix
rate.mat <- getFullMat(StateMats)
}
\references{
Tarasov, S. 2019. Integration of Anatomy Ontologies and Evo-Devo Using Structured Markov Models Suggests a New Framework for Modeling Discrete Phenotypic Traits. Systematic Biology.
}
\author{James D. Boyko}
\keyword{models}
