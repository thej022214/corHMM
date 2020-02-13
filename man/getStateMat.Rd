\name{getStateMat}
\alias{getStateMat}
\title{Basic index matrix maker}
\description{Generates an index matrix that can be manipulated to create custom models}
\usage{
getStateMat(nState)
}       
\arguments{
\item{nState}{specifies the number of states.}
}
\details{
By itself, an nState index matrix is the same as running corHMM with model = "ARD" and a single rate category. This function is best used in combination with other functions: dropStateMatPars, equateStateMatPars, and getFullMat.
}
\value{
Returns a rate matrix index
}
\seealso{
\code{\link{dropStateMatPars}},
\code{\link{equateStateMatPars}},
\code{\link{getFullMat}}
}
\examples{
#Generate a matrix for 3 traits:
StateMat <- getStateMat(3)
}
\author{James D. Boyko}
\keyword{models}
