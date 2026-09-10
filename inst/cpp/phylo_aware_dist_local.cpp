// Optional compiled version of corHMM::phylo_aware_dist_local().
//
// corHMM ships no compiled code, so this file is NOT built when the package is
// installed. Compile it on demand with corHMM::phylo_aware_dist_fast(), which
// requires the Rcpp package and a working C++ toolchain. The result is
// numerically identical to the pure R implementation, roughly three times
// faster, and worth the trouble only for large samples of reconstructions.

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double phylo_aware_dist_local_cpp(LogicalVector mismatch, NumericMatrix phylo_dist) {
  int n = mismatch.size();
  std::vector<int> indices;
  for (int i = 0; i < n; i++) {
    if (mismatch[i]) indices.push_back(i);
  }
  int m = indices.size();
  if (m <= 1) return 0.0;
  double sum = 0.0;
  for (int i = 0; i < m - 1; i++) {
    for (int j = i + 1; j < m; j++) {
      sum += exp(1.0 / phylo_dist(indices[i], indices[j]));
    }
  }
  return sum;
}
