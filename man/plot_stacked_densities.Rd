\name{plot_stacked_densities}
\alias{plot_stacked_densities}
\title{Stacked log-likelihood densities by cluster}
\description{Draws one density curve per cluster of sampled joint reconstructions, stacked vertically on a shared log-likelihood axis. Clusters whose densities overlap are statistically comparable alternative histories; a cluster shifted well to the left is a worse explanation of the data.}
\usage{
plot_stacked_densities(lnliks, clusters, cols, cluster_labels = NULL,
best_lnlik = NULL, xlab = "log-likelihood", ...)
}
\arguments{
\item{lnliks}{Numeric vector of log-likelihoods, one per sampled reconstruction, as returned in \code{lnliks} by \code{\link{compute_joint_ci}}.}
\item{clusters}{Vector of cluster assignments, one per reconstruction. Reconstructions assigned to cluster \code{0} (noise) are dropped.}
\item{cols}{Vector of fill colours, one per cluster.}
\item{cluster_labels}{Optional character vector of names for the clusters, drawn on the left axis. Defaults to the cluster numbers.}
\item{best_lnlik}{Optional log-likelihood of the best joint reconstruction, drawn as a dashed vertical reference line.}
\item{xlab}{Label for the log-likelihood axis.}
\item{...}{Further arguments passed to \code{plot}, e.g. \code{main}.}
}
\value{Called for its side effect. Returns \code{NULL} invisibly (and with a message) if only one cluster is present.}
\seealso{\code{\link{compute_joint_ci}}}
\examples{ \dontrun{
plot_stacked_densities(ci$lnliks, clusters, cols = cluster_colors)
}
}
\author{James D. Boyko}
\keyword{plot}
