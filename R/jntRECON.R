############################################ compute uncertainty by sampling
compute_joint_ci <- function(res, batch_size = 100, max_samples = 1000, break_threshold = 0, remove_outliers = TRUE, ncores = 1){
  p <- MatrixToPars(res)
  phy <- res$phy
  phy$node.label <- NULL
  
  corData <- corProcessData(res$data, collapse = res$collapse)
  model.set.final <- rate.cat.set.corHMM.JDB(phy=phy, data=res$data, rate.cat=res$rate.cat, ntraits = NULL, model = "ARD", rate.mat=res$index.mat, collapse = res$collapse)
  index_mat <- model.set.final$index.matrix
  index_mat[is.na(index_mat)] <- 0
  Q <- matrix(0, nrow=nrow(index_mat), ncol=ncol(index_mat))
  Q[index_mat > 0] <- p[index_mat[index_mat > 0]]
  diag(Q) <- -rowSums(Q)
  # pre calculate all necessary matrix expm
  p_mat <- vapply(phy$edge.length, function(x) expm(Q * x, method = "Ward77"), FUN.VALUE = matrix(0, nrow(Q), ncol(Q)))
  
  # you've tried the best
  best_joint <- ancRECON_internal(phy = phy, data = res$data, corData = corData, 
    model.set.final = model.set.final, 
    p_mat = p_mat, p = p, method = "joint", 
    rate.cat = res$rate.cat, ntraits = NULL, 
    rate.mat = res$index.mat, root.p = res$root.p, 
    get.likelihood = FALSE)
  phy$node.label <- best_joint$lik.anc.states
  
  best_lnlik <- ancRECON_internal(phy = phy, data = res$data, corData = corData, 
    model.set.final = model.set.final, 
    p_mat = p_mat, p = p, method = "joint", 
    rate.cat = res$rate.cat, ntraits = NULL, 
    rate.mat = res$index.mat, root.p = res$root.p, 
    get.likelihood = TRUE)
  
  # now try the rest
  unique_joints <- list()
  lnliks <- c()
  total_samples <- 0
  break_ratio <- Inf
  while (total_samples < max_samples) {
    new_joints <- parallel::mclapply(1:batch_size, function(x) sample_joint(res, corData, model.set.final, p_mat, x), mc.cores = ncores)
    total_samples <- total_samples + batch_size
    new_lnliks <- sapply(new_joints, function(x) x$lnlik)
    new_joints_filtered <- new_joints[!duplicated(unlist(lapply(new_joints, "[[", "lnlik")))]
    init_number <- length(new_joints_filtered)
    unique_joints <- c(unique_joints, new_joints_filtered)
    lnliks <- c(lnliks, new_lnliks[!duplicated(new_lnliks)])
    added_number <- length(unique_joints)
    unique_joints <- unique_joints[!duplicated(lnliks)]
    lnliks <- lnliks[!duplicated(lnliks)]
    final_number <- length(unique_joints)
    if(total_samples > batch_size){
      break_ratio <- (added_number - final_number)/init_number
    }
    if (break_ratio <= break_threshold) {
      break
    }
  }
  
  unique_joints <- unique_joints[!lnliks == -Inf]
  lnliks <- lnliks[!lnliks == -Inf]
  
  if(remove_outliers){
    # remove bad recons for increased consistency
    outliers <- boxplot(lnliks, plot = FALSE)$out
    if(length(outliers) > 0){
      iqr_index <- match(outliers, lnliks)
      unique_joints <- unique_joints[-iqr_index]
      lnliks <- lnliks[-iqr_index]
    }
  }
  
  state_df <- as.data.frame(do.call(rbind, lapply(unique_joints, "[[", "anc.states")))
  colnames(state_df) <- (Ntip(phy)+1):(Ntip(phy)+Nnode(phy))
  tip_state_df <- as.data.frame(do.call(rbind, lapply(unique_joints, "[[", "tip.states")))
  colnames(tip_state_df) <- phy$tip.label
  
  # Return results
  return(list(
    best_joint = best_joint,
    best_lnlik = best_lnlik,
    sampled_joints = unique_joints,
    state_df=state_df,
    tip_state_df=tip_state_df,
    lnliks = lnliks
  ))
}

sample_joint <- function(res, corData, model.set.final, p_mat, x){
  p <- MatrixToPars(res)
  res$phy$node.label <- NULL
  states_1 <- ancRECON_internal(phy = res$phy, data = res$data, corData = corData, 
    model.set.final = model.set.final, 
    p_mat = p_mat, p = p, method = "joint_unc", 
    rate.cat = res$rate.cat, ntraits = NULL, 
    rate.mat = res$index.mat, root.p = res$root.p, 
    get.likelihood = FALSE)
  res$phy$node.label <- states_1$lik.anc.states
  lnlik <- ancRECON_internal(phy = res$phy, data = res$data, corData = corData, 
    model.set.final = model.set.final, 
    p_mat = p_mat, p = p, method = "joint", 
    rate.cat = res$rate.cat, ntraits = NULL, 
    rate.mat = res$index.mat, root.p = res$root.p, 
    get.likelihood = TRUE)
  return(list(tip.states = states_1$lik.tip.states, anc.states = states_1$lik.anc.states, lnlik=lnlik))
}

compute_joint_results <- function(sim_result) {
  res <- sim_result$res
  lnliks <- sim_result$lnliks
  tru <- sim_result$data$NodeStates
  marginal_recon <- sim_result$recons$marginal_recon
  
  many_joints <- compute_joint_ci(
    res, 
    batch_size = 100, 
    max_samples = 1000, 
    ncores = 10
  )
  
  perc_correct <- unlist(lapply(
    many_joints$sampled_joints, 
    function(x) sum(x$anc.states == tru) / length(tru)
  )) * 100
  
  return(list(
    lnliks = many_joints$lnliks,
    perc_correct = perc_correct,
    true_lnlik = lnliks[3],
    best_marginal_lnlik = lnliks[2],
    bestjoint_lnlik = lnliks[1]
  ))
}

############################################ phylo distance matrix
Rcpp::cppFunction('
double phylo_aware_dist_local_cpp(LogicalVector mismatch, NumericMatrix phylo_dist) {
  int n = mismatch.size();
  std::vector<int> indices;
  for(int i = 0; i < n; i++) {
    if(mismatch[i]) indices.push_back(i);
  }
  int m = indices.size();
  if(m <= 1) return 0.0;
  double sum = 0.0;
  for(int i = 0; i < m-1; i++) {
    for(int j = i+1; j < m; j++) {
      sum += exp(1.0 / phylo_dist(indices[i], indices[j]));
    }
  }
  return sum;
}')

phylo_aware_dist_local <- function(graph_labels_x, graph_labels_y, phylo_dist) {
  mismatch <- graph_labels_x != graph_labels_y
  return(phylo_aware_dist_local_cpp(mismatch, phylo_dist))
}

############################################ analysis
node_importance_by_cluster <- function(state_df, clusters) {
  unique_clusters <- unique(clusters)
  n_nodes <- ncol(state_df)
  
  # For each cluster
  cluster_stats <- lapply(unique_clusters, function(focal_cluster) {
    # For each node
    node_stats <- lapply(1:n_nodes, function(node) {
      # Get most common state in this cluster for this node
      cluster_states <- state_df[clusters == focal_cluster, node]
      modal_state <- as.numeric(names(which.max(table(cluster_states))))
      
      # Calculate proportion of matches within cluster
      within_match <- mean(cluster_states == modal_state)
      
      # Calculate proportion of matches outside cluster
      other_states <- state_df[clusters != focal_cluster, node]
      other_match <- mean(other_states == modal_state)
      
      # Return stats for this node
      data.frame(
        node = node,
        cluster = focal_cluster,
        modal_state = modal_state,
        within_match = within_match,
        other_match = other_match,
        # importance = within_match * (1 - other_match)
        importance = (within_match^2) * ((1 - other_match)^2)  # High when consistent in cluster but different outside
        # importance = exp((within_match-1)) * exp(-other_match)
      )
    })
    
    do.call(rbind, node_stats)
  })
  
  # Combine all clusters
  do.call(rbind, cluster_stats)
}

############################################ viz
plot_stacked_densities <- function(lnliks, clusters, cols, ...) {
  lnliks <- lnliks[clusters != 0]
  clusters <- clusters[clusters != 0]
  
  x_range <- range(lnliks)
  n_clusters <- length(unique(clusters))
  if(n_clusters == 1){
    print("Only one cluster.")
    return(NULL)
  }
  clusters <- factor(clusters, levels = 1:n_clusters)
  
  plot(NULL, xlim = x_range, ylim = c(0, n_clusters), 
    xlab = "Log Likelihood", ylab = "Cluster", yaxt = "n", ...)
  
  axis(2, at = seq(0.5, n_clusters-0.5, 1), labels = levels(clusters))
  
  for(i in seq_along(levels(clusters))) {
    cluster_data <- lnliks[clusters == levels(clusters)[i]]
    if(length(cluster_data) > 0) {
      d <- density(cluster_data)
      scaled_y <- (d$y/max(d$y)) * 0.9
      polygon(d$x, scaled_y + i - 1, col = cols[i], border = cols[i])
    }
  }
}

plot_cluster_analysis <- function(tsne_result, clusters) {
  cluster_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628")
  par(mfrow = c(1,2))
  plot(tsne_result$Y, pch = 21, bg = cluster_colors[clusters], 
    xlab = "t-SNE Dim 1", ylab = "t-SNE Dim 2")
  return(cluster_colors)
}

get_xx_yy <- function (x, type = "phylogram", use.edge.length = TRUE, node.pos = NULL, 
  show.tip.label = TRUE, show.node.label = FALSE, edge.color = NULL, 
  edge.width = NULL, edge.lty = NULL, node.color = NULL, node.width = NULL, 
  node.lty = NULL, font = 3, cex = par("cex"), adj = NULL, 
  srt = 0, no.margin = FALSE, root.edge = FALSE, label.offset = 0, 
  underscore = FALSE, x.lim = NULL, y.lim = NULL, direction = "rightwards", 
  lab4ut = NULL, tip.color = par("col"), plot = TRUE, rotate.tree = 0, 
  open.angle = 0, node.depth = 1, align.tip.label = FALSE, 
  ...) {
  Ntip <- length(x$tip.label)
  if (Ntip < 2) {
    warning("found fewer than 2 tips in the tree")
    return(NULL)
  }
  .nodeHeight <- function(edge, Nedge, yy) .C(node_height, 
    as.integer(edge[, 1]), as.integer(edge[, 2]), as.integer(Nedge), 
    as.double(yy))[[4]]
  .nodeDepth <- function(Ntip, Nnode, edge, Nedge, node.depth) .C(node_depth, 
    as.integer(Ntip), as.integer(edge[, 1]), as.integer(edge[, 
      2]), as.integer(Nedge), double(Ntip + Nnode), as.integer(node.depth))[[5]]
  .nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge, 
    edge.length) .C(node_depth_edgelength, as.integer(edge[, 
      1]), as.integer(edge[, 2]), as.integer(Nedge), as.double(edge.length), 
      double(Ntip + Nnode))[[5]]
  Nedge <- dim(x$edge)[1]
  Nnode <- x$Nnode
  if (any(x$edge < 1) || any(x$edge > Ntip + Nnode)) 
    stop("tree badly conformed; cannot plot. Check the edge matrix.")
  ROOT <- Ntip + 1
  type <- match.arg(type, c("phylogram", "cladogram", "fan", 
    "unrooted", "radial", "tidy"))
  direction <- match.arg(direction, c("rightwards", "leftwards", 
    "upwards", "downwards"))
  if (is.null(x$edge.length)) {
    use.edge.length <- FALSE
  }
  else {
    if (use.edge.length && type != "radial") {
      tmp <- sum(is.na(x$edge.length))
      if (tmp) {
        warning(paste(tmp, "branch length(s) NA(s): branch lengths ignored in the plot"))
        use.edge.length <- FALSE
      }
    }
  }
  if (is.numeric(align.tip.label)) {
    align.tip.label.lty <- align.tip.label
    align.tip.label <- TRUE
  }
  else {
    if (align.tip.label) 
      align.tip.label.lty <- 3
  }
  if (align.tip.label) {
    if (type %in% c("unrooted", "radial") || !use.edge.length || 
        is.ultrametric(x)) 
      align.tip.label <- FALSE
  }
  if (type %in% c("unrooted", "radial") || !use.edge.length || 
      is.null(x$root.edge) || !x$root.edge) 
    root.edge <- FALSE
  phyloORclado <- type %in% c("phylogram", "cladogram", "tidy")
  horizontal <- direction %in% c("rightwards", "leftwards")
  if (type == "tidy" && any(x$edge.length < 0)) 
    stop("cannot plot in tidy mode with negative branch lengths. Check 'edge.length' vector.")
  xe <- x$edge
  if (phyloORclado) {
    phyOrder <- attr(x, "order")
    if (is.null(phyOrder) || phyOrder != "cladewise") {
      x <- reorder(x)
      if (!identical(x$edge, xe)) {
        ereorder <- match(x$edge[, 2], xe[, 2])
        if (length(edge.color) > 1) {
          edge.color <- rep(edge.color, length.out = Nedge)
          edge.color <- edge.color[ereorder]
        }
        if (length(edge.width) > 1) {
          edge.width <- rep(edge.width, length.out = Nedge)
          edge.width <- edge.width[ereorder]
        }
        if (length(edge.lty) > 1) {
          edge.lty <- rep(edge.lty, length.out = Nedge)
          edge.lty <- edge.lty[ereorder]
        }
      }
    }
    yy <- numeric(Ntip + Nnode)
    TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
    yy[TIPS] <- 1:Ntip
  }
  getStringLengthbyTip <- function(x, lab, sin, cex) {
    s <- strwidth(lab, "inches", cex = cex)
    lim <- getLimit(x, lab, sin, cex)
    alp <- lim/sin
    s * alp
  }
  getLimit <- function(x, lab, sin, cex) {
    s <- strwidth(lab, "inches", cex = cex)
    if (any(s > sin)) 
      return(1.5 * max(x))
    Limit <- 0
    while (any(x > Limit)) {
      i <- which.max(x)
      alp <- x[i]/(sin - s[i])
      Limit <- x[i] + alp * s[i]
      x <- x + alp * s
    }
    Limit
  }
  z <- reorder(x, order = "postorder")
  if (phyloORclado) {
    if (is.null(node.pos)) 
      node.pos <- if (type == "cladogram" && !use.edge.length) 
        2
    else 1
    if (node.pos == 1) {
      yy <- .nodeHeight(z$edge, Nedge, yy)
    }
    else {
      ans <- .C(node_height_clado, as.integer(Ntip), as.integer(z$edge[, 
        1]), as.integer(z$edge[, 2]), as.integer(Nedge), 
        double(Ntip + Nnode), as.double(yy))
      xx <- ans[[5]] - 1
      yy <- ans[[6]]
    }
    if (!use.edge.length) {
      if (node.pos != 2) 
        xx <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, 
          node.depth) - 1
      xx <- max(xx) - xx
    }
    else {
      xx <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge, 
        z$edge.length)
    }
    if (type == "tidy") {
      if (!show.tip.label) {
        yy <- tidy.xy(z$edge, Ntip, Nnode, xx, yy)
      }
      else {
        xx.tips <- xx[1:Ntip]
        pin1 <- par("pin")[1]
        lab.strlength <- getStringLengthbyTip(xx.tips, 
          x$tip.label, pin1, cex)
        xx2 <- xx
        xx2[1:Ntip] <- xx2[1:Ntip] + lab.strlength
        yy <- tidy.xy(z$edge, Ntip, Nnode, xx2, yy)
      }
    }
  }
  else {
    twopi <- 2 * pi
    rotate.tree <- twopi * rotate.tree/360
    if (type != "unrooted") {
      TIPS <- x$edge[which(x$edge[, 2] <= Ntip), 2]
      xx <- seq(0, twopi * (1 - 1/Ntip) - twopi * open.angle/360, 
        length.out = Ntip)
      theta <- double(Ntip)
      theta[TIPS] <- xx
      theta <- c(theta, numeric(Nnode))
    }
    switch(type, fan = {
      theta <- .nodeHeight(z$edge, Nedge, theta)
      if (use.edge.length) {
        r <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, 
          Nedge, z$edge.length)
      } else {
        r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
        max_r <- max(r)
        r <- (max_r - r + 1)/max_r
      }
      theta <- theta + rotate.tree
      if (root.edge) r <- r + x$root.edge
      xx <- r * cos(theta)
      yy <- r * sin(theta)
    }, unrooted = {
      nb.sp <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
      XY <- if (use.edge.length) unrooted.xy(Ntip, Nnode, 
        z$edge, z$edge.length, nb.sp, rotate.tree) else unrooted.xy(Ntip, 
          Nnode, z$edge, rep(1, Nedge), nb.sp, rotate.tree)
      xx <- XY$M[, 1] - min(XY$M[, 1])
      yy <- XY$M[, 2] - min(XY$M[, 2])
    }, radial = {
      r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
      r[r == 1] <- 0
      r <- 1 - r/Ntip
      theta <- .nodeHeight(z$edge, Nedge, theta) + rotate.tree
      xx <- r * cos(theta)
      yy <- r * sin(theta)
    })
  }
  if (phyloORclado) {
    if (!horizontal) {
      tmp <- yy
      yy <- xx
      xx <- tmp - min(tmp) + 1
    }
    if (root.edge) {
      if (direction == "rightwards") 
        xx <- xx + x$root.edge
      if (direction == "upwards") 
        yy <- yy + x$root.edge
    }
  }
  if (no.margin) 
    par(mai = rep(0, 4))
  if (show.tip.label) 
    nchar.tip.label <- nchar(x$tip.label)
  max.yy <- max(yy)
  if (is.null(x.lim)) {
    if (phyloORclado) {
      if (horizontal) {
        xx.tips <- xx[1:Ntip]
        if (show.tip.label) {
          pin1 <- par("pin")[1]
          tmp <- getLimit(xx.tips, x$tip.label, pin1, 
            cex)
          tmp <- tmp + label.offset
        }
        else tmp <- max(xx.tips)
        x.lim <- c(0, tmp)
      }
      else {
        x.lim <- c(1, max(xx[1:Ntip]))
      }
    }
    else switch(type, fan = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
            cex)
        x.lim <- range(xx) + c(-offset, offset)
      } else x.lim <- range(xx)
    }, unrooted = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
            cex)
        x.lim <- c(0 - offset, max(xx) + offset)
      } else x.lim <- c(0, max(xx))
    }, radial = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.03 * cex)
        x.lim <- c(-1 - offset, 1 + offset)
      } else x.lim <- c(-1, 1)
    })
  }
  else if (length(x.lim) == 1) {
    x.lim <- c(0, x.lim)
    if (phyloORclado && !horizontal) 
      x.lim[1] <- 1
    if (type %in% c("fan", "unrooted") && show.tip.label) 
      x.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy * 
          cex)
    if (type == "radial") 
      x.lim[1] <- if (show.tip.label) 
        -1 - max(nchar.tip.label * 0.03 * cex)
    else -1
  }
  if (phyloORclado && direction == "leftwards") 
    xx <- x.lim[2] - xx
  if (is.null(y.lim)) {
    if (phyloORclado) {
      if (horizontal) {
        y.lim <- c(1, max(yy[1:Ntip]))
      }
      else {
        pin2 <- par("pin")[2]
        yy.tips <- yy[1:Ntip]
        if (show.tip.label) {
          tmp <- getLimit(yy.tips, x$tip.label, pin2, 
            cex)
          tmp <- tmp + label.offset
        }
        else tmp <- max(yy.tips)
        y.lim <- c(0, tmp)
      }
    }
    else switch(type, fan = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
            cex)
        y.lim <- c(min(yy) - offset, max.yy + offset)
      } else y.lim <- c(min(yy), max.yy)
    }, unrooted = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
            cex)
        y.lim <- c(0 - offset, max.yy + offset)
      } else y.lim <- c(0, max.yy)
    }, radial = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.03 * cex)
        y.lim <- c(-1 - offset, 1 + offset)
      } else y.lim <- c(-1, 1)
    })
  }
  else if (length(y.lim) == 1) {
    y.lim <- c(0, y.lim)
    if (phyloORclado && horizontal) 
      y.lim[1] <- 1
    if (type %in% c("fan", "unrooted") && show.tip.label) 
      y.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy * 
          cex)
    if (type == "radial") 
      y.lim[1] <- if (show.tip.label) 
        -1 - max(nchar.tip.label * 0.018 * max.yy * cex)
    else -1
  }
  if (phyloORclado && direction == "downwards") 
    yy <- y.lim[2] - yy
  if (phyloORclado && root.edge) {
    if (direction == "leftwards") 
      x.lim[2] <- x.lim[2] + x$root.edge
    if (direction == "downwards") 
      y.lim[2] <- y.lim[2] + x$root.edge
  }
  asp <- if (type %in% c("fan", "radial", "unrooted")) 
    1
  else NA
  return(list(xx = xx, yy = yy))
}


add_node_boxes <- function(XX, YY, box_width = 2, nodes = NULL) {
  plot_dims <- par("usr")
  pin <- par("pin")
  width_ratio <- (plot_dims[2] - plot_dims[1])/((plot_dims[4] - plot_dims[3]) * (pin[1]/pin[2]))
  box_height <- box_width/width_ratio
  
  for(i in seq_along(XX)) {
    if(!is.null(nodes)){
      if(i %in% nodes){
        rect(XX[i] - box_width/2, YY[i] - box_height/2,
          XX[i] + box_width/2, YY[i] + box_height/2,
          border = "black", lwd = 0.5, col = "white")
      }
    }else{
      rect(XX[i] - box_width/2, YY[i] - box_height/2,
        XX[i] + box_width/2, YY[i] + box_height/2,
        border = "black", lwd = 0.5, col = "white")
    }
  }
}

add_state_boxes <- function(XX, YY, n_possible, state_df, cols, box_width = 0.5, nodes = NULL) {
  plot_dims <- par("usr")
  pin <- par("pin")
  width_ratio <- (plot_dims[2] - plot_dims[1])/((plot_dims[4] - plot_dims[3]) * (pin[1]/pin[2]))
  box_height <- box_width/width_ratio
  
  # Calculate grid dimensions to make approximately square boxes
  n_cols <- ceiling(sqrt(n_possible))
  n_rows <- ceiling(n_possible/n_cols)
  
  small_width <- box_width/n_cols
  small_height <- box_height/n_rows
  
  for(i in seq_along(XX)) {
    count <- 1
    if(!is.null(nodes)){
      if(i %in% nodes){
        for(row in n_rows:1) {
          for(col in 1:n_cols) {
            if(count <= n_possible) {
              focal_state <- state_df[count,i]
              rect(XX[i] - box_width/2 + (col-1)*small_width,
                YY[i] - box_height/2 + (row-1)*small_height,
                XX[i] - box_width/2 + col*small_width,
                YY[i] - box_height/2 + row*small_height,
                border = "black", lwd = 0.025, col = cols[focal_state])
              count <- count + 1
            }
          }
        }
      }
    }else{
      for(row in n_rows:1) {
        for(col in 1:n_cols) {
          if(count <= n_possible) {
            focal_state <- state_df[count,i]
            rect(XX[i] - box_width/2 + (col-1)*small_width,
              YY[i] - box_height/2 + (row-1)*small_height,
              XX[i] - box_width/2 + col*small_width,
              YY[i] - box_height/2 + row*small_height,
              border = "black", lwd = 0.025, col = cols[focal_state])
            count <- count + 1
          }
        }
      }
    }
  }
}

plot_clustered_phylo <- function(
    phy,                    # phylogenetic tree object
  state_df,              # data frame of state distributions
  cluster_assignments,    # vector of cluster assignments
  node_states,           # vector of node states
  important_nodes,       # list of important nodes per cluster
  cluster_colors = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628"),  # colors for clusters
  state_colors = c("white", "black"),  # colors for states
  state_labels = NULL,   # labels for states in legend
  outer_cex = 2.2,       # size of outer circle
  inner_rad = 3,         # radius of inner pie chart
  node_offset = 10,           # spacing between offset nodes
  edge_width = 0.5,      # width of tree edges
  show_legend = TRUE,    # whether to show the legend
  legend_pos = "topleft", # position of the legend
  legend_title = NULL,   # title for the legend
  phylogram_params = list(), # additional parameters for phylogram.plot
  plot_params = list(),   # additional parameters for plot.phylo
  point_params = list(pch = NULL, cex = NULL, col = NULL), # style parameters for points (pch, cex, col only)
  segment_params = list(lwd = NULL, lty = NULL), # style parameters for segments (lwd, lty only)
  ...
) {
  # Validate inputs
  if (!inherits(phy, "phylo")) {
    stop("phy must be a phylogenetic tree object")
  }
  if (!is.data.frame(state_df)) {
    stop("state_df must be a data frame")
  }
  
  # Get number of clusters
  n_clusters <- sum(unique(cluster_assignments) != 0)
  
  # Set up state colors
  # Set up state colors and labels
  if (is.null(state_labels)) {
    state_labels <- as.character(1:length(state_colors))
  }
  cols2 <- setNames(state_colors, 1:length(state_colors))
  
  # Get tree coordinates
  xxyy <- get_xx_yy(phy, 
    show.tip.label = FALSE, 
    edge.width = edge_width, 
    no.margin = TRUE, 
    direction = "downwards")
  XX <- xxyy$xx
  YY <- xxyy$yy
  
  # Process pie chart data for each cluster
  examined_nodes <- c()
  pie_dat <- list()
  for(i in 1:n_clusters) {
    focal_df <- state_df[cluster_assignments == i, ]
    focal_nodes <- important_nodes[[i]]$node
    examined_nodes <- c(examined_nodes, focal_nodes)
    pie_dat[[i]] <- t(apply(focal_df, 2, function(x) {
      c(table(factor(x, 1:length(state_colors)))/length(x))
    }))[focal_nodes,]
    rownames(pie_dat[[i]]) <- focal_nodes
  }
  
  # Adjust X coordinates for examined nodes
  for(i in seq_along(examined_nodes)) {
    focal_node <- examined_nodes[i]
    XX[XX > XX[focal_node + Ntip(phy)]] <- XX[XX > XX[focal_node + Ntip(phy)]] + node_offset
  }
  
  # Add extra cushion
  for(i in seq_along(unique(examined_nodes))) {
    focal_node <- unique(examined_nodes)[i]
    XX[XX > XX[focal_node + Ntip(phy)]] <- XX[XX > XX[focal_node + Ntip(phy)]] + node_offset
  }
  
  # Plot base tree
  plot.phylo(plot = FALSE, 
    phy, 
    show.tip.label = FALSE, 
    edge.width = edge_width, 
    no.margin = TRUE, 
    direction = "downwards", 
    x.lim = range(XX))
  
  # Plot phylogram
  # Merge default and user-specified phylogram parameters
  default_phylogram_params <- list(
    edge = phy$edge,
    Ntip = Ntip(phy),
    Nnode = Nnode(phy),
    xx = XX,
    yy = YY,
    horizontal = FALSE,
    edge.color = "black",
    edge.width = edge_width,
    edge.lty = 1
  )
  phylogram_params <- modifyList(default_phylogram_params, phylogram_params)
  do.call(phylogram.plot, phylogram_params)
  
  # Add legend if requested
  if (show_legend) {
    legend_args <- list(
      x = legend_pos,
      pch = 21,
      pt.bg = cols2,
      legend = state_labels
    )
    if (!is.null(legend_title)) {
      legend_args$title <- legend_title
    }
    do.call(legend, legend_args)
  }
  
  # Plot connecting lines between nodes
  nodes_visited <- c()
  point_list <- vector("list", length(pie_dat))
  for(i in seq_along(pie_dat)) {
    tmp <- c()
    for(j in 1:nrow(pie_dat[[i]])) {
      node_number <- as.numeric(rownames(pie_dat[[i]]))[j]
      nodes_visited <- c(nodes_visited, node_number)
      xx <- XX[node_number + Ntip(phy)]
      yy <- YY[node_number + Ntip(phy)]
      current_offset <- node_offset * sum(nodes_visited == node_number)
      tmp <- rbind(tmp, 
        data.frame(node = node_number, 
          xx = xx + current_offset, 
          yy = yy, 
          col = cluster_colors[i]))
    }
    point_list[[i]] <- tmp
  }
  
  # Draw connecting lines
  for(i in seq_along(point_list)) {
    pnt <- point_list[[i]]
    for(j in 1:nrow(pnt)) {
      node_j <- pnt$node[j]
      as_anc <- phy$edge[,1] == node_j + Ntip(phy)
      as_dec <- phy$edge[,2] == node_j + Ntip(phy)
      anc_nodes <- phy$edge[as_dec,1]
      dec_nodes <- phy$edge[as_anc,2]
      main_xy <- unlist(pnt[j, c('xx', 'yy')])
      for(k in seq_along(anc_nodes)) {
        if(anc_nodes[k] %in% (pnt$node + Ntip(phy))) {
          anc_xy <- pnt[(pnt$node + Ntip(phy)) == anc_nodes[k],]
          anc_xy <- c(xx = anc_xy$xx, yy = anc_xy$yy)
        } else {
          anc_xy <- c(xx = XX[anc_nodes[k]],
            yy = YY[anc_nodes[k]])
        }
        # Apply segments with only style parameters customizable
        segments(x0 = main_xy[1],
          y0 = main_xy[2],
          x1 = anc_xy[1],
          y1 = anc_xy[2],
          col = cluster_colors[i],
          lwd = segment_params$lwd %||% edge_width,
          lty = segment_params$lty %||% 2)
      }
    }
  }
  
  # Plot base nodes
  for(i in seq_along(XX)) {
    if(i > Ntip(phy)) {
      # Merge default and user-specified point parameters
      default_point_params <- list(
        x = XX[i],
        y = YY[i],
        pch = 21,
        bg = cols2[node_states[i - Ntip(phy)]]
      )
      final_point_params <- modifyList(default_point_params, point_params)
      do.call(points, final_point_params)
    }
  }
  
  # Plot cluster nodes with pie charts
  nodes_visited <- c()
  point_list <- vector("list", length(pie_dat))
  for(i in seq_along(pie_dat)) {
    tmp <- c()
    for(j in 1:nrow(pie_dat[[i]])) {
      node_number <- as.numeric(rownames(pie_dat[[i]]))[j]
      nodes_visited <- c(nodes_visited, node_number)
      xx <- XX[node_number + Ntip(phy)]
      yy <- YY[node_number + Ntip(phy)]
      current_offset <- node_offset * sum(nodes_visited == node_number)
      points(xx + current_offset, 
        yy, 
        pch = 21, 
        bg = cluster_colors[i], 
        cex = outer_cex)
      floating.pie.asp(xx + current_offset, 
        yy, 
        pie_dat[[i]][j,], 
        radius = inner_rad, 
        col = cols2)
      tmp <- rbind(tmp, 
        data.frame(node = node_number, 
          xx = xx + current_offset, 
          yy = yy, 
          col = cluster_colors[i]))
    }
    point_list[[i]] <- tmp
  }
  
  # Return the coordinates and other plotting information invisibly
  invisible(list(
    XX = XX,
    YY = YY,
    point_list = point_list,
    pie_dat = pie_dat
  ))
}