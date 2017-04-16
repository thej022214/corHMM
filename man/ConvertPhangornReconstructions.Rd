\name{ConvertPhangornReconstructions}
\alias{ConvertPhangornReconstructions}
\title{Convert phangorn reconstruction to a vector}
\usage{
ConvertPhangornReconstructions(x, site = 1, best = TRUE)
}
\arguments{
\item{x}{The phyDat object that contains a character reconstruction from phangorn}

\item{site}{The character number to convert into a vector}

\item{best}{A logical indicating whether the state that maximizes some function (likelihood, parsimony, etc.) is to be returned.}
}
\description{
Converts a character reconstruction from phangorn into a vector
}
\details{
Creates a vector that contains the best tips and node state from a phangorn reconstruction.
}

