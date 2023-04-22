\name{simMarkov}
\alias{simMarkov}
\title{Simulate a character on the tree}
\description{Simulates a character using an instanteous rate matrix, storing the states at nodes and at the tips of the tree.}
\usage{
simMarkov(phy, Q, root.freqs)
}       
\arguments{
\item{phy}{A phylogeny of class phylo.}
\item{Q}{an instantaneous rate matrix.}
\item{root.freqs}{A vector of root frequencies.}
}
\details{
This function will generate a character at tips and at interior nodes given a tree and Q matrix. 
}
\value{
A two element list containing `TipStates` and `NodeStates`
}
\author{James D. Boyko}
\keyword{models}