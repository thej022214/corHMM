\name{updateStateMats}
\alias{updateStateMats}
\title{Updates a list of index matrices}
\description{Updates the indicies from a list of index matrix for use as custom hidden Markoc models}
\usage{
updateStateMats(StateMats)
}       
\arguments{
\item{StateMats}{A list of index matrices of class matrix.}
}
\details{
Outputs a list of rate matrices to be used in getFullMat. Each index matrix in the list represents its own rate class. This function automatically takes all index matrices and gives them independent index values. updateStateMats is part of several functions which help the user efficiently create custom index matrices. The use of this function comes in combination with other functions such as getStateMat, dropStateMatPars, equateStateMatPars, and getFullMat.
}
\value{
Returns a list of index matrices.
}
\examples{
#Generate 2 matrices for 3 traits:
StateMatA <- getStateMat(3)
StateMatB <- getStateMat(3)
StateMats <- list(StateMatA, StateMatB)
StateMats_indep <- updateStateMats(StateMats)
}
\author{James D. Boyko}
\keyword{models}
