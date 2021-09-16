\name{getParamVector}
\alias{getParamVector}
\title{Returns a vector of rates from a fitted corhmm.object}
\description{This will return a vector with the transition rates between states from a fitted corhmm.object. This vector can be used to compare rates between runs or models as well as to use as input for the p argument in corHMM.}
\usage{
getParamVector(corhmm.object)
}       
\arguments{
\item{corhmm.object}{A fitted corhmm.object, output from corHMM().}
}
\author{Brian C. O'Meara}
