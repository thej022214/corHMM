\name{dropStateMatPars}
\alias{dropStateMatPars}
\title{Drop specific rates from an index matrix}
\description{Drops rates from an index matrix for use as custom models}
\usage{
dropStateMatPars(StateMat, Pars)
}       
\arguments{
\item{StateMat}{The index matrix of class matrix to be modified.}
\item{Pars}{A vector of index values which will be dropped from the index matrix.}
}
\details{
dropStateMatPars is part of several functions which help the user efficiently create custom index matrices. The use of this function comes in combination with other functions such as getStateMat, equateStateMatPars, and getFullMat.
}
\value{
Returns a rate matrix index
}
\seealso{
\code{\link{getStateMat}},
\code{\link{equateStateMatPars}},
\code{\link{getFullMat}}
}
\examples{
#Generate a matrix for 3 traits:
StateMat <- getStateMat(3)
#Drop transitions from 1 -> 3, 2 -> 1, 3 -> 2:
StateMat_constrained <- dropStateMatPars(StateMat, c(5, 1, 4))
}
\author{James D. Boyko}
\keyword{models}
