fixedObservedLabels <- function(phy, nobs) {
  labels <- phy$node.label
  if(is.null(labels)) return(NULL)
  if(length(labels) != phy$Nnode) stop("Expected one node label per internal node.", call.=FALSE)
  text <- trimws(as.character(labels))
  missing <- is.na(text) | text == ""
  if(all(missing)) return(NULL)
  states <- suppressWarnings(as.numeric(text))
  x <- states[!missing]
  if(any(!is.finite(x) | x != floor(x) | x < 1 | x > nobs))
    stop("Fixed node labels must be observed state indices from 1 to ", nobs, ".", call.=FALSE)
  states[missing] <- NA
  as.integer(states)
}

# Fitted labels identify observed states across all hidden categories.
# Explicit makeSimmap constraints instead identify individual rows of Q.
nodeStateEvidence <- function(phy, k, rate.cat, observed=NULL, nodes=NULL, states=NULL) {
  if(length(nodes) != length(states)) stop("Supply one fixed state per node.", call.=FALSE)
  if(is.null(observed) && !length(nodes)) return(NULL)
  ntip <- length(phy$tip.label)
  mask <- matrix(1, ntip + phy$Nnode, k)
  obs <- rep(seq_len(k/rate.cat), rate.cat)
  for(i in which(!is.na(observed))) mask[ntip+i,] <- obs == observed[i]
  indices <- function(x, n) is.numeric(x) && all(is.finite(x) & x==floor(x) & x>=1 & x<=n)
  if(length(nodes) && (!indices(nodes, nrow(mask)) || !indices(states, k)))
    stop("Invalid fixed node or joint state index.", call.=FALSE)
  for(i in seq_along(nodes)) {
    j <- nodes[i]
    mask[j,] <- mask[j,] * (seq_len(k) == states[i])
    if(!any(mask[j,] != 0)) stop("Conflicting state constraints at node ", j, ".", call.=FALSE)
  }
  mask
}

checkNodeMass <- function(x) {
  if(any(!is.finite(x)) || sum(x) <= 0)
    stop("Node constraints, tip data or root probabilities are incompatible.", call.=FALSE)
}
