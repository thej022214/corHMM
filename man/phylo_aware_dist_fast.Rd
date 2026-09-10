\name{phylo_aware_dist_fast}
\alias{phylo_aware_dist_fast}
\title{Compile a faster version of the reconstruction distance}
\description{Compiles and returns a C++ implementation of \code{\link{phylo_aware_dist_local}}. corHMM ships no compiled code, so the C++ source is only built when you ask for it here. The result is numerically identical to the pure R version and roughly three times faster, which is worth the extra step only when you are clustering many thousands of reconstructions.}
\usage{
phylo_aware_dist_fast(verbose = FALSE)
}
\arguments{
\item{verbose}{Logical. Pass \code{TRUE} to see the compiler output. Default \code{FALSE}.}
}
\details{Requires the \pkg{Rcpp} package and a working C++ toolchain (Rtools on Windows, Xcode command line tools on macOS). Compilation takes a few seconds the first time in a session.

The returned function is only valid in the session that created it: it holds a pointer to code compiled into that session, so do not save it to an \code{.rds} or \code{.RData} file and reload it later. Call \code{phylo_aware_dist_fast()} again in each new session instead.

The C++ source is installed with the package and can be read or adapted from \code{system.file("cpp", "phylo_aware_dist_local.cpp", package = "corHMM")}.}
\value{A function with the same arguments and return value as \code{\link{phylo_aware_dist_local}}, suitable for registering with \code{proxy::pr_DB}.}
\seealso{\code{\link{phylo_aware_dist_local}}, \code{\link{compute_joint_ci}}}
\examples{ \dontrun{
fast_dist <- phylo_aware_dist_fast()
proxy::pr_DB$set_entry(FUN = fast_dist, names = "PhyloSpreadDistFast")
D <- proxy::dist(as.matrix(ci$state_df), method = "PhyloSpreadDistFast",
phylo_dist = phylo_dist)
}
}
\author{James D. Boyko}
\keyword{models}
