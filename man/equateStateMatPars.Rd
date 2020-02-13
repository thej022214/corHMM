\name{equateStateMatPars}
\alias{equateStateMatPars}
\title{Equates specific rates within an index matrix}
\description{Equates rates within an index matrix for use in custom models}
\usage{
equateStateMatPars(StateMat, ParsList)
}       
\arguments{
\item{StateMat}{The index matrix to be modified.}
\item{ParsList}{A list of vectors containing index values which will be equated from the index matrix. }
}
\details{
The function will go through all the elements of the list equating the paramaters indicated within that element. equateStateMatPars is part of several functions which help the user efficiently create custom index matrices. The use of this function comes in combination with other functions such as getStateMat, dropStateMatPars, and getFullMat.
}
\value{
Returns a rate matrix index
}
\seealso{
\code{\link{getStateMat}},
\code{\link{dropStateMatPars}},
\code{\link{getFullMat}}
}
\examples{
#Generate a matrix for 3 traits:
StateMat <- getStateMat(3)
#Equate transitions 1 -> 2 & 2 -> 1; 1 -> 3 & 3 <- 1; 2 -> 3 & 3 -> 2:
#This is the same as setting model = "SYM" in corHMM
StateMat_constrained <- equateStateMatPars(StateMat, list(c(1,3),c(2,5),c(4,6)))
}
\author{James D. Boyko}
\keyword{models}
