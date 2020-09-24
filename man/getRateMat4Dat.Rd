\name{getStateMat4Dat}
\alias{getStateMat4Dat}
\title{Produce an index matrix and legend from a dataset}
\description{Takes a dataset to produce an index matrix that corresponds to a single state-dependent process (i.e. a single rate category) and a legend which matches input data to the rows and columns of the index matrix and corHMM solution.}
\usage{
getStateMat4Dat(data, model = "ARD", dual = FALSE)
}       
\arguments{
\item{data}{A data matrix containing species information in the same format as the main \code{corHMM} function: column[,1] is species names, column[,2:n] are the discrete states.}
\item{model}{One of "ARD", "SYM", or "ER". ARD: all rates differ. SYM: rates between any two states do not differ. ER: all rates are equal.}
\item{dual}{A boolean indicating whether or not to include dual transitions.}
}
\details{
This function will generate an index matrix based on user provided data. It provides a useful starting point for further modifications using \code{dropStateMatPars}, \code{equateStateMatPars}, and \code{getFullMat}. If more than a single column of data is given double transitions between characters are disallowed. For example, if character 1 is the presence or absence of limbs, and character 2 is the presence or absence of fingers, then the transition from absence of limbs and fingers to presence of limbs and fingers is automatically disallowed. This is consistent with Pagel's (1994) model of correlated character evolution. 
}
\value{
\item{$legend}{A named vector. The elements of the vector are all the unique state combinations in the user data. The names of the vector are the state number assigned to each combination.}
\item{$rate.mat}{A rate index matrix describing a single rate class.}
}
\seealso{
\code{getFullmat}
}
\examples{
\donttest{
data(primates)
phy <- primates[[1]]
phy <- multi2di(phy)
data <- primates[[2]]
# create a legend and rate mat from a multi-character dataset.
LegendAndRateMat <- getStateMat4Dat(data)
rate.mat <- LegendAndRateMat$rate.mat
legend <- LegendAndRateMat$legend
}
}
\references{
Pagel, M. 1994. Detecting correlated evolution on phylogenies: a gneeral method for the comparative analysis of discrete characters. Proc. R. Soc. Lond. B 255:37-45.
}
\author{James D. Boyko}
\keyword{models}