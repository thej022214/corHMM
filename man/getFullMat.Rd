\name{getFullMat}
\alias{getFullMat}
\title{Combines several rate class index matrices}
\description{Combines rate class index matrices by a RateClassMat to output a single index matrix for use in corHMM}
\usage{
getFullMat(StateMats, RateClassMat, UpdateRates = TRUE)
}       
\arguments{
\item{StateMats}{A list of index matrices of class matrix. These matrices describe individual rate classes.}
\item{RateClassMat}{An index matrix of class matrix. This matrix describes how rate class matrices are related to one another.}
\item{UpdateRates}{A logical constant indicating whether you want your distinct rate class matrices to have independent parameters.}
}
\details{
This function takes a list of rate class matrices (StateMats) and how they are related to one another (RateClassMat). getFullMat is part of several functions which help the user efficiently create custom index matrices. The use of this function comes in combination with other functions such as getStateMat, dropStateMatPars, equateStateMatPars.
}
\value{
Returns an index matrix.
}
\seealso{
\code{\link{getStateMat}},
\code{\link{dropStateMatPars}},
\code{\link{equateStateMatPars}}
}
\examples{
#Generate 2 matrices for 3 traits:
StateMatA <- getStateMat(3)
StateMatB <- getStateMat(3)
#Modify StateMatB to be identical to the "SYM" model
StateMatB_constrained <- equateStateMatPars(StateMatB, list(c(1,3),c(2,5),c(4,6)))
#Make a list of our state matrices (each representing its own rate class)
StateMats <- list(StateMatA, StateMatB_constrained)
#Determine how the rate classes are related to one another
RateClassMat <- getStateMat(2)
#Combine into a single index matrix
rate.mat <- getFullMat(StateMats, RateClassMat)
}
\author{James D. Boyko}
\keyword{models}
