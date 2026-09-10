\name{plot_clustered_phylo}
\alias{plot_clustered_phylo}
\title{Plot competing ancestral histories on a phylogeny}
\description{Draws the phylogeny with the best joint reconstruction at every internal node, and attaches to the diagnostic nodes of each cluster a pie chart of the states that cluster assigns there. This gives a single figure in which the alternative evolutionary histories implied by different clusters of reconstructions can be read off directly, rather than being blurred together as they would be in a marginal summary.}
\usage{
plot_clustered_phylo(phy, state_df, cluster_assignments, node_states,
important_nodes, cluster_colors = c("#E41A1C", "#377EB8", "#4DAF4A",
"#984EA3", "#FF7F00", "#FFFF33", "#A65628"),
state_colors = c("white", "black"), state_labels = NULL,
outer_cex = 2.2, inner_rad = 3, node_offset = 10, edge_width = 0.5,
show_legend = TRUE, legend_pos = "topleft", legend_title = NULL,
phylogram_params = list(), plot_params = list(),
point_params = list(pch = NULL, cex = NULL, col = NULL),
segment_params = list(lwd = NULL, lty = NULL), ...)
}
\arguments{
\item{phy}{A phylogeny of class \code{"phylo"}, matching the tree used to generate \code{state_df}.}
\item{state_df}{Data frame of sampled reconstructions from \code{\link{compute_joint_ci}}.}
\item{cluster_assignments}{Vector of cluster assignments, one per row of \code{state_df}.}
\item{node_states}{Integer vector of internal node states for the reconstruction drawn on the tree itself, typically \code{res$states} from the best joint reconstruction.}
\item{important_nodes}{A list with one element per cluster, in cluster order, each a data frame with a \code{node} column giving that cluster's diagnostic nodes. Usually built by thresholding \code{\link{node_importance_by_cluster}}.}
\item{cluster_colors}{Vector of colours for the clusters.}
\item{state_colors}{Vector of colours for the character states.}
\item{state_labels}{Optional character vector of state names for the legend.}
\item{outer_cex}{Size of the coloured ring drawn around each cluster pie.}
\item{inner_rad}{Radius of the pie chart drawn inside that ring.}
\item{node_offset}{Horizontal spacing used to fan out multiple clusters that flag the same node.}
\item{edge_width}{Width of the tree edges.}
\item{show_legend, legend_pos, legend_title}{Legend control.}
\item{phylogram_params, plot_params}{Named lists of additional arguments for \code{phylogram.plot} and \code{plot.phylo}.}
\item{point_params}{Named list of \code{pch}, \code{cex} and \code{col} for the base node points.}
\item{segment_params}{Named list of \code{lwd} and \code{lty} for the lines connecting cluster pies to the tree.}
\item{...}{Currently unused.}
}
\details{Because clustering has already separated the sample into distinct global narratives, the pie charts show variation remaining \emph{within} a cluster, not an average over all possible histories. Note that \code{important_nodes[[i]]} must correspond to cluster \code{i}; build the list by looping over cluster numbers rather than over \code{unique(clusters)}.}
\value{Invisibly returns a list with the plotting coordinates \code{XX} and \code{YY}, the per-cluster point positions \code{point_list}, and the pie data \code{pie_dat}.}
\seealso{\code{\link{compute_joint_ci}}, \code{\link{node_importance_by_cluster}}, \code{\link{plot_raw_recon}}}
\examples{ \dontrun{
important_nodes <- lapply(seq_len(k), function(i) {
ni <- node_importance_by_cluster(ci$state_df, clusters)
subset(ni, cluster == i & importance > 0.5)
})
plot_clustered_phylo(phy, ci$state_df, clusters, res$states, important_nodes,
cluster_colors = cols, state_colors = state_cols)
}
}
\author{James D. Boyko}
\keyword{plot}
