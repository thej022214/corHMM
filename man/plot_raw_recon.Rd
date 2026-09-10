\name{plot_raw_recon}
\alias{plot_raw_recon}
\title{Raster plot of a sample of joint reconstructions}
\description{Plots a matrix of sampled joint ancestral state reconstructions as a raster, with one row per reconstruction and one column per internal node. When cluster assignments are supplied, reconstructions are sorted by cluster and a colour bar identifies each cluster, making the "reconstruction axes" recovered by clustering directly visible.}
\usage{
plot_raw_recon(state_matrix, cluster_assignments = NULL, cluster_colors = NULL,
state_colors = NULL, state_labels = NULL, show_legend = TRUE,
xlab = "Node", ylab = "Reconstruction", main = NULL, ...)
}
\arguments{
\item{state_matrix}{A matrix (or data.frame) of sampled reconstructions, rows being reconstructions and columns being internal nodes, with entries coded as integer states \code{1:k}. Typically \code{as.matrix(x$state_df)} from \code{\link{compute_joint_ci}}.}
\item{cluster_assignments}{Optional integer vector of cluster memberships, one per row of \code{state_matrix}. A value of \code{0} marks an unclustered (noise) reconstruction; these are drawn last in grey. If \code{NULL}, all rows are drawn in their original order.}
\item{cluster_colors}{Optional vector of colours for the clusters. Defaults to \code{hcl.colors(..., "Dark 3")}.}
\item{state_colors}{Optional vector of colours for the character states. Defaults to \code{hcl.colors(k, "Zissou 1")}.}
\item{state_labels}{Optional character vector of state names used in the legend. Defaults to \code{1:k}.}
\item{show_legend}{Logical. Draw a state legend. Default \code{TRUE}.}
\item{xlab, ylab, main}{Standard plot labels.}
\item{...}{Further arguments passed to \code{image}.}
}
\value{Invisibly returns a list with \code{order}, the row ordering used, and \code{cluster}, the cluster assignments in that order.}
\seealso{\code{\link{compute_joint_ci}}, \code{\link{plot_clustered_phylo}}}
\examples{ \dontrun{
ci <- compute_joint_ci(res)
plot_raw_recon(as.matrix(ci$state_df))
}
}
\author{James D. Boyko}
\keyword{plot}
