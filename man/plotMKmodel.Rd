\name{plotMKmodel}
\alias{plotMKmodel}
\title{Plot a Markov model}
\description{Plots a diagram of a Markov model from the output of corHMM or a custom index matrix}
\usage{
plotMKmodel(pp, rate.cat = NULL, display = "column", color = c("blue", "red"), arrow.scale = 1, 
text.scale = 1, vertex.scale = 1)
}
\arguments{
\item{pp}{an object of class corHMM or matrix.}
\item{rate.cat}{if using a custom matrix then the number of rate categories must be indicated.}
\item{display}{the structure of the plot. one of "column", "square", or "row".}
\item{color}{Either, 1. a vector of 2 colors to create a gradient from low transition rates (first element) to high transition rates (second element), or 2. "col.blind" which will use the color pallete "plasma" from viridis.}
\item{arrow.scale}{determines the size of the arrows for the Markov diagram.}
\item{text.scale}{determines the size of the text for the plotted matrix.} 
\item{vertex.scale}{determines the size of the text for the Markov diagram.} 
}
\details{
Plots Markov models in a ball and stick type diagram next to its corresponding matrix. If plotting a hidden rates model it will produce a compound plot describing how the different rate classes are related to one another. If the input is a corHMM result then arrows are colored by relative rate. If the input is a custom matrix arrows are colored by the paramater index.
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
