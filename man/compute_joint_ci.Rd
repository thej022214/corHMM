\name{compute_joint_ci}
\alias{compute_joint_ci}
\title{Compute joint ancestral-state reconstructions by sampling}
\description{Compute joint ancestral state reconstructions (ASR) by sampling joint histories and returning a description of their likelihoods. The function first precomputes the necessary transition structure, then samples joint histories using one of several methods (simmap, marginal, or Pupko-based sampling), until a maximum number of samples is reached. The result is a list containing the best joint history, the corresponding log-likelihood, and the ensemble of sampled histories with their log-likelihoods, suitable for downstream uncertainty analyses.}
\usage{
compute_joint_ci(res, batch_size = 100, max_samples = 1000,
remove_outliers = TRUE, p_method = "simmap")
}
\arguments{
\item{res}{A corHMM results object. The function uses internal components of res, including:
\itemize{
\item \code{res$phy}: a phylogenetic tree (class \code{"phylo"}).
\item \code{res$data}: a data.frame of observed traits.
\item \code{res$collapse}: boolean indicating whether states were collapsed.
\item \code{res$rate.cat}: number of rate categories.
\item \code{res$index.mat}: index matrix describing parameter groupings for the model.
\item \code{res$root.p}: root prior (vector or method flag).}}
\item{batch_size}{Integer. Number of joint histories to sample per iteration. Default is 100.}
\item{max_samples}{Integer. Maximum total number of joint histories to sample. Default is 1000.}
\item{remove_outliers}{Logical. If \code{TRUE}, remove outliers in the log-likelihood distribution using a boxplot-based rule. Default is \code{TRUE}.}
\item{p_method}{Character. Sampling strategy. One of \code{"simmap"} (default), \code{"marginal"}, or \code{"pupko"}.}
}
\value{
A list with components:
\itemize{
\item{best_joint: }{The internal representation of the best joint history as returned by \code{ancRECON_internal(..., method = "joint", get.likelihood = FALSE)}.}
\item{best_lnlik: }{The log-likelihood corresponding to \code{best_joint} (obtained with \code{get.likelihood = TRUE}).}
\item{sampled_joints: }{A list of sampled joint histories. Each element is expected to contain at least \code{anc.states} (internal node states) and \code{tip.states} (tip states).}
\item{state_df:}{Data frame of internal node states across the sampled histories.}
\item{tip_state_df: }{Data frame of tip states across the sampled histories.}
\item{lnliks: }{Numeric vector of log-likelihoods for the sampled histories.}
}
}
\details{
This function implements a two-pass sampling approach to quantify the uncertainty of joint ancestral-state reconstructions. Key steps include:
\itemize{
\item{Computing p from the input result object and constructing the Q-matrix, then precomputing the transition probability matrices for each edge length.}
\item{Determining the best joint history and its log-likelihood.}
\item{Obtaining a conditional likelihood object for stochastic mapping-based sampling.}
\item{Iteratively sampling batches of joint histories using the specified \code{p_method}:
\itemize{
\item{\code{"simmap"}: sample by simulating substitution histories with conditional likelihoods.}
\item{\code{"marginal"}: sample marginal reconstructions and combine them into joint histories.}
\item{\code{"pupko"}: use Pupko-based sampling of joint histories.}}
}
\item{Aggregating sampled histories, removing duplicates, and filtering out (-\infty) ln-likelihoods.}
\item{Optionally removing outlier ln-likelihoods to improve stability.}
\item{Assembling \code{state_df} and \code{tip_state_df} from the remaining histories.} }

}
\examples{ \donttest{
res is a corHMM-style result object produced earlier in your workflow
ci <- compute_joint_ci(res)
str(ci)
}
}
\author{James D. Boyko}
\keyword{models}
