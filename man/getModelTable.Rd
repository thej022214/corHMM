\name{getModelTable}
\alias{getModelTable}
\title{Summarize Model Statistics for a List of corHMM Objects}
\description{
\code{getModelTable} extracts key statistics from a list of fitted \code{corHMM} models and returns a summary table. The table includes the number of parameters, log-likelihood, and model selection criteria such as AIC, delta AIC, and AIC weights.
}
\usage{
getModelTable(model_list, type = "AIC")
}
\arguments{
  \item{model_list}{
    A list of \code{corHMM} objects, each representing a fitted hidden Markov model to phylogenetic and character data.
  }
  \item{type}{
    The type of model selection criterion to use. Options are \code{"AIC"}, \code{"AICc"}, or \code{"BIC"}. Default is \code{"AIC"}.
  }
}
\details{
This function takes a list of models fitted using \code{corHMM} and calculates key statistics for comparison across models. Specifically, it calculates the number of parameters, log-likelihood, the chosen model selection criterion (e.g., AIC), the difference in the criterion relative to the best model (delta AIC), and the relative model weight based on the criterion.

\code{getModelTable} can handle different model selection criteria such as AIC, AICc, and BIC by specifying the \code{type} argument.
}
\value{
A data frame with the following columns:
\itemize{
  \item \code{np}: The number of parameters in the model.
  \item \code{lnLik}: The log-likelihood of the model.
  \item \code{AIC} (or the value of \code{type}): The model selection criterion value.
  \item \code{dAIC}: The difference in the criterion between the model and the best model (i.e., the model with the minimum criterion value).
  \item \code{AICwt}: The Akaike weights, representing the relative likelihood of the model given the data.
}
}
\examples{
\donttest{
# Assuming you have a list of fitted corHMM models:
models <- list(model1, model2, model3)
getModelTable(models)

# To use BIC instead of AIC:
getModelTable(models, type = "BIC")
}
}
\seealso{
\code{\link{corHMM}} for fitting hidden Markov models to phylogenetic data.
}
\author{
James D. Boyko
}
