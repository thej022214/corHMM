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
Outputs an index of the rate parameters which are to be optimized in corHMM. By itself, an nState index matrix is the same as running corHMM with model = "ARD" and a single rate category. The use of this function comes in combination with other functions such as dropStateMatPars, equateStateMatPars, updateStateMats, and getFullMat.
}
\value{
Returns a rate matrix index
}
\examples{
#Generate a matrix for 3 traits:
StateMat <- getStateMat(3)
}
\author{James D. Boyko}
\keyword{models}
