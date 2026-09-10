\name{node_importance_by_cluster}
\alias{node_importance_by_cluster}
\title{Identify nodes that distinguish clusters of reconstructions}
\description{For every internal node and every cluster of sampled joint reconstructions, scores how diagnostic that node is of that cluster. High scores mark the nodes at which alternative evolutionary histories actually disagree, and are the nodes worth plotting and interpreting.}
\usage{
node_importance_by_cluster(state_df, clusters)
}
\arguments{
\item{state_df}{Data frame of sampled reconstructions, rows being reconstructions and columns being internal nodes, as returned in \code{state_df} by \code{\link{compute_joint_ci}}.}
\item{clusters}{Vector of cluster assignments, one per row of \code{state_df}.}
}
\details{For a focal cluster and node, the modal state within the cluster is found. \code{within_match} is the proportion of reconstructions inside the cluster carrying that state and \code{other_match} is the proportion outside the cluster carrying it. Importance is \code{within_match^2 * (1 - other_match)}, which is large only when a node is both consistent within the cluster and different from everything else.}
\value{A data frame with one row per node per cluster and columns \code{node}, \code{cluster}, \code{modal_state}, \code{within_match}, \code{other_match} and \code{importance}.}
\seealso{\code{\link{plot_clustered_phylo}}}
\examples{ \dontrun{
ni <- node_importance_by_cluster(ci$state_df, clusters)
subset(ni, importance > 0.5)
}
}
\author{James D. Boyko}
\keyword{models}
