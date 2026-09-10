\name{phylo_aware_dist_local}
\alias{phylo_aware_dist_local}
\title{Phylogenetically aware distance between two reconstructions}
\description{Computes a dissimilarity between two joint ancestral state reconstructions that accounts for where on the phylogeny they disagree. Reconstructions that differ at several phylogenetically adjacent nodes are treated as more dissimilar than reconstructions that differ at the same number of scattered nodes, so that the resulting distance matrix separates genuinely different evolutionary narratives rather than isolated node flips.}
\usage{
phylo_aware_dist_local(graph_labels_x, graph_labels_y, phylo_dist)
}
\arguments{
\item{graph_labels_x}{Integer vector of internal node states for the first reconstruction.}
\item{graph_labels_y}{Integer vector of internal node states for the second reconstruction, in the same node order.}
\item{phylo_dist}{A square matrix of pairwise node distances among internal nodes, e.g. \code{ape::dist.nodes} applied to a tree with all branch lengths set to 1, subset to the internal nodes.}
}
\details{Let \eqn{V_{diff}} be the set of internal nodes at which the two reconstructions differ. If \eqn{|V_{diff}| \le 1} the distance is 0. Otherwise the distance is \eqn{\sum \exp(1 / d_P(v_a, v_b))} summed over all unique pairs of nodes in \eqn{V_{diff}}, where \eqn{d_P} is the node distance supplied in \code{phylo_dist}.

The function is intended to be registered with the \pkg{proxy} package and used through \code{proxy::dist}, which supplies the pairwise arguments (see examples).}
\value{A single non-negative numeric distance.}
\seealso{\code{\link{compute_joint_ci}}, \code{\link{node_importance_by_cluster}}}
\examples{ \dontrun{
tmp <- phy; tmp$edge.length <- rep(1, nrow(tmp$edge))
phylo_dist <- ape::dist.nodes(tmp)[-(1:Ntip(tmp)), -(1:Ntip(tmp))]
proxy::pr_DB$set_entry(FUN = phylo_aware_dist_local, names = "PhyloSpreadDistLocal")
D <- proxy::dist(as.matrix(ci$state_df), method = "PhyloSpreadDistLocal",
phylo_dist = phylo_dist)
}
}
\author{James D. Boyko}
\keyword{models}
