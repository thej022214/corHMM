\name{simCharacterHistory}
\title{Simulate a Markov model}
\description{Simulates a Markov or Hidden Markov model given a phylogeny, a rate matrix, and intial root frequency.}
\usage{
simCharacterHistory(phy, Q, root.freqs, Q2 = NA, NoI = NA)
}       
\arguments{
\item{phy}{A phylogeny of class phylo.}
\item{Q}{A rate matrix.}
\item{root.freqs}{Starting frequencies of the states at the root.}
\item{Q2}{An optional rate matrix which will be applied at NoI.}
\item{NoI}{The "node of interest". If you have given a Q2 that is different from Q1, starting at NoI Q2 will be applied at that node and every descendent branch.}
}
\details{
This function will simulate a dataset according to a given rate matrix. It will return a list of states at the tips and states at the nodes. This is a stochastic process so it will return different tip and node states each time it is run, even if the model specified is the same.
}
\value{
Returns a list of tip states and node states.
}
\author{James D. Boyko and Jeremy Beaulieu}
\keyword{models}
