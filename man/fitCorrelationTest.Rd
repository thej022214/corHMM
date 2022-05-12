\name{fitCorrelationTest}
\alias{fitCorrelationTest}
\title{Test for correlation}
\description{Automatically fits a set of independent and dependent models to test for correlation between characters.}
\usage{
  fitCorrelationTest(phy, data, simplified_models=FALSE)
}
\arguments{
  \item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
  \item{data}{a data.frame containing species information. The first column must be species names matching the phylogeny. Additional columns contain discrete character data.}
  \item{simplified_models}{A boolean which indicates whether to include simplified independent and dependent models (currently only works for two binary-state characters; see Details).}
}
\details{
This function automatically fit a set of multi-rate independent and dependent models (with default corHMM options) to drastically reduce false support for correlation. Currently, the simplified models are only available for two binary-state characters, but it is straightforward for users to use the tools available in corHMM to create model structures specific to their questions when the datasets are more complex. 

The correlation between two characters is often interpreted as evidence that there exists a significant and biologically important relationship between them. However, Maddison and FitzJohn (2015) recently pointed out that in certain situations find evidence of correlated evolution between two categorical characters is often spurious, particularly, when the dependent relationship stems from a single replicate deep in time. In Boyko and Beaulieu (in prep) we showed that there is, in fact, a statistical solution to the problem posed by Maddison and FitzJohn (2015) naturally embedded within the expanded model space afforded by the hidden Markov model (HMM) framework. 

The problem of single unreplicated evolutionary events manifests itself as rate heterogeneity within our models and that this is the source of the false correlation. Therefore, we argue that this problem is better understood as model misspecification rather than a failure of comparative methods to account for phylogenetic pseudoreplication. We utilize HMMs to develop a multi-rate independent model which, when implemented, drastically reduces support for correlation.
}
\value{
\code{fitCorrelationTest} returns an object of class \code{corhmm_list}. This is a list with elements:
\item{$independent_model_fit}{A corHMM object of the standard independent model ala Pagel (1994).}
\item{$correlated_model_fit}{A corHMM object of the standard dependent model ala Pagel (1994).}
\item{$hidden_Markov_independent_model_fit}{A corHMM object of the hidden Markov independent model which allows for rate heterogeneity independent of the focal character.}
\item{$hidden_Markov_correlated_model_fit.cat}{A corHMM object of the hidden Markov dependent model which allows for rate heterogeneity independent of the focal character as well as correlation between characters.}
\item{$simplified_*}{If simplified was set to TRUE, then the function will also return simplified versions of the above models. These models have fewer parameters than the above models while still being either dependent or independent models.}
}
\examples{
\donttest{
data(primates)
phy <- multi2di(primates[[1]])
data <- primates[[2]]
# not run because of long run times
#corr_test_fits <- fitCorrelationTest(phy = phy, data = data, simplified_models = TRUE)
#corr_test_fits
}
}
\references{
Maddison W.P., FitzJohn R.G. 2015. The Unsolved Challenge to Phylogenetic Correlation Tests for Categorical Characters. Syst Biol. 64:127-136.
}
\author{James D. Boyko}
