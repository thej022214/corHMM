######################################################################################################################################
######################################################################################################################################
### Internal helpers for transition probability matrices and tree traversal
######################################################################################################################################
######################################################################################################################################

# Q stays fixed while we walk over every branch of the tree, so decomposing it
# once and reusing the decomposition is far cheaper than calling expm() per
# branch: P(t) = V %*% diag(exp(lambda * t)) %*% solve(V).
#
# Returns NULL when Q is defective or the decomposition is too inaccurate to
# trust, which tells the caller to fall back on expm().
getQdecomp <- function(Q, tol = 1e-8){
  ev <- try(eigen(Q), silent = TRUE)
  if(inherits(ev, "try-error")){
    return(NULL)
  }
  Vi <- try(solve(ev$vectors), silent = TRUE)
  if(inherits(Vi, "try-error")){
    return(NULL)
  }
  # V %*% diag(values) %*% Vi must reproduce Q; a defective Q gives a nearly
  # singular V and a large residual here.
  resid <- max(abs(ev$vectors %*% (ev$values * Vi) - Q))
  scale <- max(1, max(abs(Q)))
  if(!is.finite(resid) || resid > tol * scale){
    return(NULL)
  }
  list(values = ev$values, vectors = ev$vectors, inv = Vi)
}

# P(t) from a decomposition built by getQdecomp().
expmAt <- function(dec, t){
  Re(dec$vectors %*% (exp(dec$values * t) * dec$inv))
}

# P(t) %*% x without ever forming P(t) -- two matrix-vector products instead of
# two matrix-matrix products.
expmAtv <- function(dec, t, x){
  Re(dec$vectors %*% (exp(dec$values * t) * (dec$inv %*% x)))
}

# Builds a P(t) closure for Q, using the eigen decomposition when it is
# trustworthy and expm() otherwise. `reference` is an already computed
# expm(Q * ref.t) used to cross-check the decomposition; pass NULL to skip.
# Set clamp = TRUE when the result is used as a probability vector (e.g. fed to
# sample.int, which rejects negative weights): both this and expm() can return
# entries a hair below zero through rounding.
makeExpmFuns <- function(Q, reference = NULL, ref.t = 1, tol = 1e-8, clamp = FALSE){
  dec <- getQdecomp(Q, tol = tol)
  if(!is.null(dec) && !is.null(reference)){
    if(max(abs(expmAt(dec, ref.t) - reference)) > 1e-6){
      dec <- NULL
    }
  }
  if(is.null(dec)){
    P <- function(t) expm(Q * t, method = c("Ward77"))
    Pv <- function(t, x) expm(Q * t, method = c("Ward77")) %*% x
    exact <- FALSE
  }else{
    P <- function(t) expmAt(dec, t)
    Pv <- function(t, x) expmAtv(dec, t, x)
    exact <- TRUE
  }
  if(clamp){
    inner.P <- P
    inner.Pv <- Pv
    P <- function(t){ out <- inner.P(t); out[out < 0] <- 0; out }
    Pv <- function(t, x){ out <- inner.Pv(t, x); out[out < 0] <- 0; out }
  }
  list(P = P, Pv = Pv, exact = exact)
}

# The pruning loops repeatedly ask "which edges descend from this node?".
# Doing that with which(phy$edge[,1] == focal) inside the loop is O(nodes *
# edges); this builds the whole mapping in one pass instead.
#
# Returns a list parallel to `anc` (the unique ancestors in phy$edge[,1] order),
# indexed positionally so lookups are O(1).
getDesRows <- function(edge1, anc){
  split(seq_along(edge1), factor(match(edge1, anc), levels = seq_along(anc)))
}

# Every node except the root has exactly one parent edge, so which(phy$edge[,2]
# == node) can be a single lookup instead of a scan over all edges. Returns a
# vector indexed by node number; the root maps to 0.
getRowOfChild <- function(edge2, n.total){
  out <- integer(n.total)
  out[edge2] <- seq_along(edge2)
  out
}

# Vectorized version of the per-tip tip fog loops. A tip whose row is all 1s is
# an unknown state and is left alone; every other tip has its observed state(s)
# discounted by the fog probability and its zeros filled in with it.
# `divisor` is rate.cat when there is more than one rate category, else 1.
applyTipFog <- function(liks, nb.tip, tip.fog, divisor = 1){
  idx <- seq_len(nb.tip)
  tipm <- liks[idx, , drop = FALSE]
  is.one <- tipm == 1
  is.zero <- tipm == 0
  active <- rowSums(is.zero) > 0
  if(!any(active)){
    return(liks)
  }
  is.one[!active, ] <- FALSE
  is.zero[!active, ] <- FALSE
  fog.sum <- as.vector((tipm != 1) %*% tip.fog) / divisor
  tipm[is.one] <- matrix(1 - fog.sum, nb.tip, ncol(tipm))[is.one]
  tipm[is.zero] <- matrix(tip.fog, nb.tip, ncol(tipm), byrow = TRUE)[is.zero]
  liks[idx, ] <- tipm
  liks
}

# Column-wise / row-wise positive sums of Q, used for the equilibrium root
# frequencies. Equivalent to the per-column which() loop it replaces.
getEquilRoot <- function(Q){
  pos <- Q * (Q >= 0)
  rowsum <- colSums(pos)
  colsum <- rowSums(pos)
  rowsum / (rowsum + colsum)
}
