\name{getCVTable}
\alias{getCVTable}
\title{Print Method for corhmm.kfold Class Objects}
\description{
Summarize the scores for each fold and the average cross-validation score for objects of class \code{"corhmm.kfold"}.
}
\usage{
getCVTable(x)
}
\arguments{
  \item{x}{
    An object of class \code{"corhmm.kfold"} returned from \code{kFoldCrossValidation}.
  }
}
\value{
Returns a list of scores per fold for all lambdas and an average score for each lambda.  
}
\author{
James D. Boyko
}
