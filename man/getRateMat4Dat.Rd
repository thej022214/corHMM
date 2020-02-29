\name{getRateMat4Dat}
\alias{getRateMat4Dat}
\title{Produce a index matrix and legend from a dataset}
\description{Takes a dataset and phylogeny to produce a legend and index matrix that can be manipulated and given to corHMM}
\usage{
getRateMat4Dat(data, model = "ARD")
}       
\arguments{
\item{data}{a data matrix containing species information.}
\item{model}{One of "ARD", "SYM", or "ER". ARD: all rates differ. SYM: rates between any two states do not differ. ER: all rates are equal.}
}
\details{
Outputs an index rate matrix to be used in corHMM. This function takes a list of rate class matrices (StateMats) and how they are related to one another (RateClassMat). getFullMat is part of several functions which help the user efficiently create custom index matrices. The use of this function comes in combination with other functions such as getStateMat, dropStateMatPars, equateStateMatPars, and updateStateMats.
}
\value{
\item{$legend}{A named vector. The elements of the vector are all the unique state combinations in the user data. The names of the vector are the state number assigned to each combination.}
\item{$rate.mat}{A rate index matrix.}
}
\examples{
data(primates)
phy <- primates[[1]]
phy <- multi2di(phy)
data <- primates[[2]]
getRateMat4Dat(data)
}
\author{James D. Boyko}
\keyword{models}
