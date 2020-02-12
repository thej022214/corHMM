\name{getFullMat}
\alias{getFullMat}
\title{Creates the index matrix for use in corHMM}
\description{Combines the rate class index matrices and how they are related and outputs a single index matrix for use in corHMM}
\usage{
getFullMat(StateMats, RateClassMat)
}       
\arguments{
\item{StateMats}{A list of index matrices of class matrix. These matrices describe individual rate classes}
\item{RateClassMat}{An index matrices of class matrix. This matrix describes how rate class matrices are related to one another.}
}
\details{
Outputs an index rate matrix to be used in corHMM. This function takes a list of rate class matrices (StateMats) and how they are related to one another (RateClassMat). getFullMat is part of several functions which help the user efficiently create custom index matrices. The use of this function comes in combination with other functions such as getStateMat, dropStateMatPars, equateStateMatPars, and updateStateMats.
}
\value{
Returns a list of index matrices.
}
\examples{
#Generate 2 matrices for 3 traits:
StateMatA <- getStateMat(3)
StateMatB <- getStateMat(3)
#Modify StateMatB to be identical to the "SYM" model
StateMatB_constrained <- equateStateMatPars(StateMatB, list(c(1,3),c(2,5),c(4,6)))
#Update StateMatA and StateMatB_constrained so their parameters are independently estimated
StateMats <- list(StateMatA, StateMatB_constrained)
StateMats_indep <- updateStateMats(StateMats)
#Determine how they rate classes are related to one another
RateClassMat <- getStateMat(2)
#Combine into a single index matrix
rate.mat <- getFullMat(StateMats_indep, RateClassMat)
}
\author{James D. Boyko}
\keyword{models}
