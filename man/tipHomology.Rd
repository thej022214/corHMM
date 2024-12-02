\name{tipHomology}
\alias{tipHomology}
\title{The Probability of a Tip State Being Homologous}
\description{
This function takes a fitted model (in the form of a corHMM object) and returns the probability that different tips are homologous with respect to a specified node. Two types of "homology" are included: the first is the probability that there was no change throughout the lineage's entire history, and the second is the probability that the starting and ending states were the same.
}
\usage{
tipHomology(corhmm_obj, type = "strict", node = NULL, return.likelihoods = FALSE)
}
\arguments{
  \item{corhmm_obj}{
    A corHMM fit returned from \code{corHMM}.
  }
  \item{type}{
    Either "strict" or "wagner". "strict" specifies the probability of no change throughout a lineage's history, while "wagner" returns the probability of starting and ending in the same state.
  }
  \item{node}{
    The node number from which to calculate the probability of homology. If NULL, the root is used. Only tips descending from the node are returned (a tip state cannot be homologous with a node it never descended from).
  }
  \item{return.likelihoods}{
    Whether to return the exact probability of the scenario or the probability relative to the likelihood of the dataset and model. 
  }
}
\value{
Returns a vector of homology probabilities.  
}
\details{
The probability of homology is defined as either the probability that a tip inherited its state without ever changing throughout its history ("type = strict") or that the specified ancestor and the tip had the same state ("type = wagner"). 

In the first case, the probability is calculated by fixing all branches between the tip and the specified node to a particular state. For example, if the tip state is 0, all branches between the tip and specified node are treated as being in state 0. Other branches not between the focal tip and the specified ancestral node are treated as unknown. 

In the "wagner" homology case, the probability is calculated as the likelihood of the tip and the ancestor sharing the same state. For example, if the tip state is 0, this corresponds to assuming that all nodes along the path between the tip and the ancestral node are also in state 0.

To build intuition, consider the probability of homology as being similar to ancestral state reconstruction. It integrates over all possibilities for branches and nodes not explicitly constrained. Another perspective is to interpret these probabilities as inversely proportional to the number of trait changes. This is particularly intuitive in the "strict" case, where the probability directly reflects the likelihood of no changes occurring along the lineage.

One could use a likelihood ratio test (LRT) to evaluate evidence for at least one change (using \code{chisq} with one degree of freedom). However, since the null hypothesis assumes no change, the LRT can only fail to reject strict homology; it cannot provide evidence for significant homology.
}
\examples{
data(primates)
phy <- multi2di(primates[[1]])
data <- primates[[2]]
MK_3state <- corHMM(phy = phy, data = data, rate.cat = 1)

# Calculate homology probabilities for all tips
homol_p <- tipHomology(MK_3state, type = "strict")

# Specify a particular node to calculate probabilities for its descendants
homol_p_node <- tipHomology(MK_3state, type = "strict", node = 63)
}
\author{
James D. Boyko
}
