\name{plotMKmodel}
\alias{plotMKmodel}
\title{Plot a Markov model}
\description{Will plots a diagram of a Markov model from the output of corHMM or a custom index matrix}
\usage{
plotMKmodel(pp, rate.cat = NULL, col.func = plasma, arrow.scale = 2, text.scale = 3)
}
\arguments{
\item{pp}{an object of class corHMM or matrix.}
\item{rate.cat}{if using a custom matrix then the number of rate categories must be indicated.}
\item{col.func}{the type of color gradient to be used (see: ?viridis for all options).}
\item{arrow.scale}{size of the arrows for the Markov diagram.}
\item{text.scale}{size of the text for the plotted matrix.} 
}
\details{
Plots Markov models in a ball and stick type diagram next to matricies. If plotting a hidden rates model it will produce a compound plot describing how the different rate classes are related to one another. If the input is a corHMM result then arrows are colored by relative rate. If the input is a custom matrix arrows are colored by the paramater index.
}
\value{
Returns a ball and stick diagram of the input model.  
}
\examples{
## Not run
# data(primates)
## Obtain the fit of second rate class underlying a binary character:
# pp<-corHMM(primates$tree,primates$trait[,c(1,2)],rate.cat=2,node.states="none")
# plotMKmodel(pp)
}
\author{James D. Boyko}
\keyword{plot}
