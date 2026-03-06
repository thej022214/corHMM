## utility function, copied from lme4
## unnamed elements of ... will be given names from deparse()
namedList <- function (...)  {
    L <- list(...)
    snm <- sapply(substitute(list(...)), deparse)[-1]
    if (is.null(nm <- names(L))) 
        nm <- snm
    if (any(nonames <- nm == "")) 
        nm[nonames] <- snm[nonames]
    setNames(L, nm)
}

mkdev.corhmm_rtmb <- function(p, phy, liks, Q, rate, root.p, rate.cat, order.test, lewis.asc.bias, set.fog=FALSE, fog.vec, pen.type=NULL, lambda=1) {
  ## isolate computations that depend explicitly on p within the RTMB functions
  ## do as much computation as possible (i.e. computation that depends only on other components)
  ##   outside the core function; pass data to core function through tmb_data (not strictly necessary,
  ##   but helps with clarity/modularization)
  ## core function **cannot** contain any hard if() statements (e.g. various bits here that return
  ##  <very-large-number> if some 'bad' condition is met
  
  if (lewis.asc.bias) stop("can't do lewis.asc.bias yet; recursive call to likelihood function ...")
  if (order.test) stop("can't do order.test (non-differentiable ...)")
  if (set.fog) warning("set.fog is untested in RTMB implementation")
  
  nb.node <- Nnode(phy)
  nb.tip <- Ntip(phy)
  TIPS <- seq.int(nb.tip)
  anc <- unique(phy$edge[,1])
  k.rates <- dim(Q)[2] / 2
  
  ## Pre-compute static penalty indices outside prune_fun (no AD tracing needed)
  ## These depend only on rate/rate.cat structure, not on p
  pen_indices <- if (!is.null(pen.type) && rate.cat > 1) {
    rate_class_names <- paste0("R", 1:rate.cat)
    lapply(rate_class_names, function(rc) grep(rc, colnames(rate)))
  } else {
    NULL
  }
  
  ## off-diagonal indices for Q (static structure): used by penalty
  ## We need to know which Q entries are off-diagonal non-zero to compute means
  ## These are positions where rate[] value is <= np (i.e., a real parameter)
  ## We pre-compute this mask outside AD
  nq <- nrow(Q)
  
  ## precompute fog application indices outside prune_fun (static, no AD)
  fog_idx <- if (set.fog) {
    lapply(seq_len(nb.tip), function(tip.index) {
      list(
        zeros = which(liks[tip.index, ] == 0),
        ones  = which(liks[tip.index, ] == 1),
        not1  = which(liks[tip.index, ] != 1)
      )
    })
  } else {
    NULL
  }
  
  tmb_data <- namedList(nb.node, nb.tip, TIPS, anc, k.rates, nq, fog_idx)
  
  prune_fun <- function(pars) {
    "[<-" <- RTMB::ADoverload("[<-")
    "c" <- RTMB::ADoverload("c")
    "diag<-" <- RTMB::ADoverload("diag<-")
    RTMB::getAll(pars, tmb_data)
    p <- exp(p)
    cp_root.p <- root.p
    comp <- numeric(nb.tip + nb.node)
    
    if (set.fog) {
      fog_pars    <- seq.int(length(unique(fog.vec)))
      tip.fog.tmp <- p[fog_pars]
      p           <- p[-fog_pars]
      tip.fog     <- numeric(length(fog.vec))
      tip.fog[]   <- c(tip.fog.tmp, 0)[fog.vec]
      
      for (tip.index in seq(nb.tip)) {
        idx <- fog_idx[[tip.index]]
        if (length(idx$zeros) > 0) {
          if (rate.cat > 1) {
            liks[tip.index, idx$ones] <- 1 - (sum(tip.fog[idx$not1]) / rate.cat)
          } else {
            liks[tip.index, idx$ones] <- 1 - sum(tip.fog[idx$not1])
          }
          liks[tip.index, idx$zeros] <- tip.fog[idx$zeros]
        }
      }
    }
    
    np <- length(p)
    Q <- rate
    Q[Q <= np] <- p[Q[Q <= np]]
    Q[rate == np + 1] <- 0
    
    for (i in 1:nrow(Q)) {
      Q[i,i] <- -1 * sum(Q[i,])
    }
    
    ## ---------------------------------------------------------------
    ## Penalty score - computed on the AD-traced Q, fully differentiable
    ## Mirrors get_penalty_score() logic but without value-dependent branches
    ## ---------------------------------------------------------------
    pen_score <- 0
    if (!is.null(pen.type)) {
      if (rate.cat == 1) {
        ## diagonal of Q holds -rowSums; mean(-diag(Q)) == mean(rowSums of off-diag)
        diag_vals <- numeric(nq)
        for (i in 1:nq) diag_vals[i] <- Q[i, i]
        
        if (pen.type == "l1") {
          pen_score <- mean(-diag_vals)           # mean of row exit rates
        } else if (pen.type == "l2") {
          pen_score <- mean(diag_vals^2)
        } else if (pen.type == "er") {
          ## original: sd(p) if length(p)>1 else 0
          ## sd via mean-of-squared-deviations - differentiable
          if (np > 1) {
            pen_score <- mean((p - mean(p))^2)    # variance (proportional to sd^2)
          }
          ## if np <= 1, pen_score stays 0
        }
      } else {
        ## multi-rate-class: sum penalties per rate class using pre-built indices
        pen_score_vec <- numeric(rate.cat)
        for (rc_i in seq_len(rate.cat)) {
          idx <- pen_indices[[rc_i]]
          ## Extract the off-diagonal sub-block for this rate class
          ## Q[idx, idx] has diagonal already set; we want only off-diagonal entries
          ## off-diag mask was computed outside (static structure)
          rc_offdiag <- numeric(0)
          for (r in idx) {
            for (cc in idx) {
              if (r != cc) rc_offdiag <- c(rc_offdiag, Q[r, cc])
            }
          }
          if (pen.type == "l1") {
            pen_score_vec[rc_i] <- mean(rc_offdiag)
          } else if (pen.type == "l2") {
            pen_score_vec[rc_i] <- mean(rc_offdiag^2)
          } else if (pen.type == "er") {
            n_rc <- length(rc_offdiag)
            if (n_rc > 1) {
              pen_score_vec[rc_i] <- mean((rc_offdiag - mean(rc_offdiag))^2)
            }
          }
        }
        pen_score <- sum(pen_score_vec)
      }
    }
    ## ---------------------------------------------------------------
    
    for (i in seq(from = 1, length.out = nb.node)) {
      focal <- anc[i]
      desRows <- which(phy$edge[,1] == focal)
      desNodes <- phy$edge[desRows, 2]
      v <- 1
      for (desIndex in seq_along(desRows)) {
        v <- drop(as.matrix(v * Matrix::expm(Q * phy$edge.length[desRows[desIndex]]) %*% liks[desNodes[desIndex],]))
      }
      
      if (!is.null(phy$node.label)) {
        if (!is.na(phy$node.label[focal - nb.tip])) {
          fixer.tmp <- numeric(dim(Q)[2] / rate.cat)
          fixer.tmp[phy$node.label[focal - nb.tip]] <- 1
          fixer <- rep(fixer.tmp, rate.cat)
          v <- v * fixer
        }
      }
      
      comp[focal] <- sum(v)
      liks[focal, ] <- v / comp[focal]
    }
    
    root <- nb.tip + 1L
    
    equil.root <- numeric(ncol(Q))
    QQ <- Q
    diag(QQ) <- NA
    for (i in 1:ncol(Q)) {
      rowsum <- sum(Q[, i], na.rm = TRUE)
      colsum <- sum(Q[i, ], na.rm = TRUE)
      equil.root[i] <- rowsum / (rowsum + colsum)
    }
    
    if (is.null(root.p)) {
      flat.root <- equil.root
      k.rates <- 1 / length(which(!is.na(equil.root)))
      flat.root[!is.na(flat.root)] <- k.rates
      flat.root[is.na(flat.root)] <- 0
      loglik <- -(sum(log(comp[-TIPS])) + log(sum(flat.root * liks[root,])))
    }
    
    if (is.character(root.p)) {
      if (root.p == "yang") {
        root.p <- Matrix::expm(10000 * Q)[1,]
        loglik <- -(sum(log(comp[-TIPS])) + log(sum(root.p * liks[root,])))
      } else {
        root.p <- liks[root,] / sum(liks[root,])
        loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p) + log(liks[root,])))))
      }
    } else {
      if (is.numeric(root.p[1])) {
        loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p) + log(liks[root,])))))
      }
    }
    
    if (lewis.asc.bias) {
      p <- log(p)
      dummy.liks.vec <- getLewisLikelihood(p = p, phy = phy, liks = liks, Q = Q, rate = rate, root.p = cp_root.p, rate.cat = rate.cat)
      loglik <- loglik - log(sum(root.p * (1 - exp(dummy.liks.vec))))
    }
    
    ## Apply penalty - additive, zero when pen.type is NULL (pen_score==0)
    loglik <- loglik + (pen_score * lambda)
    
    return(loglik)
  }
  
  RTMB::MakeADFun(prune_fun, list(p = p), silent = TRUE)
}