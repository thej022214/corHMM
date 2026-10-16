# Endpoint-conditioned CTMC paths by uniformization (Hobolth & Stone 2009),
# essentially stolen from OUwie's hOUwie.bridge.R. R = I + Q / mu is a discrete chain
# embedded in a rate-mu Poisson process. Self-transitions are virtual events.
# Its matrix powers are shared by all branches and maps for a given Q.
makeUniformizationTable <- function(Q, max_branch_length, min.prob,
                                    tol=1e-12, max_events_cap=10000L){
  n <- nrow(Q)
  if(!is.matrix(Q) || !is.numeric(Q) || n != ncol(Q) || n == 0L ||
     any(!is.finite(Q)) || any(Q[row(Q) != col(Q)] < 0) ||
     any(diag(Q) > 0) || any(abs(rowSums(Q)) > 1e-10 * pmax(1, -diag(Q)))){
    stop("Uniformization requires a valid transition rate matrix.", call. = FALSE)
  }
  mu <- max(-diag(Q))
  if(mu == 0){
    return(list(mu=0, R=diag(n),
                powers=array(diag(n), c(n, n, 1L)), max_events=0L))
  }
  lambda <- mu * max_branch_length
  n.needed <- if(is.finite(lambda)){
    stats::qpois(log(tol) + log(min.prob), lambda, lower.tail=FALSE, log.p=TRUE) + 2
  }else{
    Inf
  }
  if(!is.finite(n.needed) || n.needed > max_events_cap ||
     n * n * (n.needed + 1) > 2e7){
    stop("Uniformization table is too large; use method = 'rejection'.", call. = FALSE)
  }
  n.needed <- as.integer(n.needed)
  R <- diag(n) + Q / mu
  R[R < 0] <- 0
  powers <- array(0, c(n, n, n.needed + 1L))
  powers[, , 1L] <- diag(n)
  for(i in seq_len(n.needed)){
    powers[, , i + 1L] <- powers[, , i] %*% R
  }
  list(mu=mu, R=R, powers=powers, max_events=n.needed)
}

sampleBranchBridge <- function(init, final, blen, uni){
  impossible <- function(){
    stop("The branch endpoints have zero probability under the supplied model.", call. = FALSE)
  }
  if(blen == 0 || uni$mu == 0){
    if(init != final) impossible()
    return(stats::setNames(blen, init))
  }
  grid <- 0:uni$max_events
  log.weights <- stats::dpois(grid, uni$mu * blen, log=TRUE) +
    log(uni$powers[init, final, grid + 1L])
  largest <- max(log.weights)
  if(!is.finite(largest)) impossible()
  n.events <- sample.int(length(grid), 1L, prob=exp(log.weights - largest)) - 1L
  if(n.events == 0L){
    return(stats::setNames(blen, init))
  }
  times <- sort(stats::runif(n.events, 0, blen))
  states <- integer(n.events + 1L)
  states[1L] <- init
  states[n.events + 1L] <- final
  if(n.events > 1L){
    for(i in seq_len(n.events - 1L)){
      step <- uni$R[states[i], ] * uni$powers[, final, n.events - i + 1L]
      if(!isTRUE(sum(step) > 0)) impossible()
      states[i + 1L] <- sample.int(nrow(uni$R), 1L, prob=step)
    }
  }
  runs <- rle(states)
  cumulative <- cumsum(diff(c(0, times, blen)))
  durations <- diff(c(0, cumulative[cumsum(runs$lengths)]))
  stats::setNames(durations, runs$values)
}
