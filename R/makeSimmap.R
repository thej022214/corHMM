makeSimmap <- function(tree, data, model, rate.cat, root.p="yang", nSim=1, nCores=1, fix.node=NULL, fix.state=NULL, parsimony = FALSE, max.attempt = 100000, collapse=TRUE){
  if(any(tree$edge.length<=.Machine$double.eps)){
    warning(paste0("Branch lengths of 0 detected. Adding ", sqrt(.Machine$double.eps)), immediate. = TRUE)
    tree$edge.length <- tree$edge.length + sqrt(.Machine$double.eps) 
  }
  model[is.na(model)] <- 0
  diag(model) <- 0
  diag(model) <- -rowSums(model)
  conditional.lik <- getConditionalNodeLik(tree, data, model, rate.cat, root.p, parsimony = parsimony, collapse=collapse)
  if(length(fix.node) != length(fix.state)){
    stop("The number of nodes supplied to be fixed does not match the number of states provided.",.call=FALSE)
  }
  if(!is.null(fix.state)){
    if(max(fix.state) > dim(model)[1]){
      stop("One of the states being fixed does not exist in this model.")
    }
  }
  if(!is.null(fix.node) & !is.null(fix.state)){
    # test if we are fixing an external or internal node
    for(i in 1:length(fix.node)){
      focal.fix.node <- fix.node[i]
      focal.fix.state <- fix.state[i]
      if(focal.fix.node <= length(tree$tip.label)){
        conditional.lik$tip.states[focal.fix.node,] <- 0
        conditional.lik$tip.states[focal.fix.node, focal.fix.state] <- 1
      }else{
        fix.internal <- focal.fix.node - length(tree$tip.label)
        conditional.lik$node.states[fix.internal,] <- 0
        conditional.lik$node.states[fix.internal, focal.fix.state] <- 1
      }
    }
  }
  maps <- simSubstHistory(tree, conditional.lik$tip.states, conditional.lik$node.states, model, nSim, nCores, max.attempt)
  mapped.edge <- lapply(maps, function(x) convertSubHistoryToEdge(tree, x))
  obj <- vector("list", nSim)
  CorData <- corProcessData(data, collapse = collapse)
  if(collapse){
    StateNames <- rep(gsub("_", "|", CorData$ObservedTraits), rate.cat)
    RCNames <- rep(paste("R", 1:rate.cat, sep = ""), each = length(CorData$ObservedTraits))
  }else{
    StateNames <- rep(gsub("_", "|", CorData$PossibleTraits), rate.cat)
    RCNames <- rep(paste("R", 1:rate.cat, sep = ""), each = length(CorData$PossibleTraits))
  }
  if(rate.cat > 1){
    StateNames <- paste(RCNames, StateNames)
  }
  for(i in 1:nSim){
    tree.simmap <- tree
    tree.simmap$maps <- maps[[i]]
    tree.simmap$maps <- lapply(maps[[i]], function(x) correctMapName(x, StateNames))
    if(dim(mapped.edge[[i]])[2] == length(StateNames)){
      tree.simmap$mapped.edge <- mapped.edge[[i]]
      colnames(tree.simmap$mapped.edge) <- StateNames
    }else{ # some states are missing from the simulation
      missing_states <- which(!(1:length(StateNames) %in% as.numeric(colnames(mapped.edge[[i]]))))
      mapped.edge[[i]] <- cbind(mapped.edge[[i]], matrix(0, dim(mapped.edge[[i]])[1], length(missing_states), dimnames = list(rownames(mapped.edge[[i]]), missing_states)))
      tree.simmap$mapped.edge <- mapped.edge[[i]]
      colnames(tree.simmap$mapped.edge) <- StateNames[as.numeric(colnames(tree.simmap$mapped.edge))]
    }
    tree.simmap$mapped.edge <- tree.simmap$mapped.edge[,match(StateNames, colnames(tree.simmap$mapped.edge))]
    tree.simmap$Q <- model
    colnames(tree.simmap$Q) <- rownames(tree.simmap$Q) <- StateNames
    attr(tree.simmap, "map.order") <- "right-to-left"
    if (!inherits(tree.simmap, "simmap")) 
      class(tree.simmap) <- c("simmap", setdiff(class(tree.simmap), "simmap"))
    obj[[i]] <- tree.simmap
  }
  if(nSim > 1){
    class(obj) <- c("multiSimmap", "multiPhylo")
  }
  return(obj)
}

correctMapName <- function(map_element, state_names){
  names(map_element) <- state_names[as.numeric(names(map_element))]
  return(map_element)
}

# simulate ancestral states at each internal node
# simSingleCharHistory<- function(phy, model, cladewise.index, state.probability, state.sample){
#   # sample the root state
#   anc.index <- phy$edge[cladewise.index[1],1]
#   state.sample[anc.index,] <- t(rmultinom(1, 1, state.probability[anc.index,]))
#   
#   # sample the nodes
#   for(i in cladewise.index){
#     # A = the probability of transition FROM node i (ancestor)
#     anc.index <- phy$edge[i,1]
#     dec.index <- phy$edge[i,2]
#     A <- state.sample[anc.index,] %*% expm(model * phy$edge.length[i])
#     # B = the conditional probability of node j (descendent)
#     B <- state.probability[dec.index,]
#     # The decendant is sampled from p/sum(p), where p = A * B
#     p <- A * B
#     state.sample[dec.index,] <- t(rmultinom(1, 1, p))
#   }
#   return(state.sample)
# }

# simulate the character history eq 3 of bollback
simCharHistory <- function(phy, Pj, root, model, vector.form = FALSE){
  # organize
  cladewise.index <- reorder.phylo(phy, "cladewise", index.only = TRUE)
  lnliks <- numeric(length(cladewise.index))
  # combine data in order of the edge matrix
  state.sample <- matrix(0, dim(Pj)[3], dim(Pj)[2])
  # sample the root
  root.index <- phy$edge[cladewise.index[1],1]
  state.sample[root.index,] <- c(rmultinom(1, 1, root))
  # single sim function, not needed since we already calculated all possible probablities. now we just need to sample from those probabilities given the ancestral node.
  # state.samples <- simSingleCharHistory(phy, model, cladewise.index, state.probability, state.sample)
  for(i in cladewise.index){
    # A = the probability of transition FROM node i (ancestor)
    anc.index <- phy$edge[i,1]
    dec.index <- phy$edge[i,2]
    p <- Pj[,,i][which(state.sample[anc.index,] == 1),]
    # if all p = 0, then we have an incompatible character history (possible with hidden states and directional models) and we have to try another character history.
    if(all(p == 0)){
      state.sample <- simCharHistory(phy, Pj, root, model, vector.form = FALSE)
      return(state.sample)
    }
    state_i <- c(rmultinom(1, 1, p))
    lnliks[i] <- log(p[state_i==1])
    state.sample[dec.index,] <- state_i
  }
  if(!vector.form){
    state.sample <- apply(state.sample, 1, function(y) which(y == 1))
  }
  return(state.sample)
}

simBranchSubstHistory <- function(init, final, total.bl, model, d.rates, max.attempt=1000){
  current.bl <- c()
  current.state <- init
  restart <- TRUE
  attempt.no <- 0
  while(restart == TRUE){
    if(attempt.no >= max.attempt){
      ViableShortPath <- FloydWalshAlg(model, init, final)
      current.bl <- rep(total.bl/length(ViableShortPath), length(ViableShortPath))
      names(current.bl) <- ViableShortPath
      return(current.bl)
    }
    attempt.no <- attempt.no + 1
    # draw a rate and waiting time based on init
    rate <- d.rates[current.state]
    # if the rate is 0, then we are in a sink state and we will remain in that state for the entire branch length
    if(rate == 0){
      waiting.time <- total.bl
    }else{
      waiting.time <- rexp(1, rate)
    }
    current.bl <- c(current.bl, waiting.time)
    names(current.bl)[length(current.bl)] <- current.state
    # if reach the max attempts, draw a path of viable character states then use a uniform distribution to split the times
    # if the waiting time is smaller than the branch
    if(sum(current.bl) < total.bl){
      # then: a substitution is drawn with probability of leaving that state
      # then: a new wating time is simulated
      p <- model[current.state,]
      p[p < 0] <- 0
      current.state <- which(rmultinom(1, 1, p) == 1)
    }else{
      # if the waiting time is longer than the branch length AND the states are NOT the same
      if(current.state != final){
        current.bl <- c()
        current.state <- init
        # if the waiting time is longer than the branch length AND the states are the same
      }else{
        # if the vector is one current bl is total bl
        if(length(current.bl) == 1){
          current.bl <- total.bl
          names(current.bl) <- current.state
        }else{
          current.bl[length(current.bl)] <- total.bl - sum(current.bl[-length(current.bl)])
        }
        restart <- FALSE
      }
    }
  }
  return(current.bl)
}

simSingleSubstHistory <- function(cladewise.index, CharHistory, phy, model, max.attempt = 1000){
  SubstHistory <- vector("list", length(cladewise.index))
  d.rates <- -diag(model)
  for(i in cladewise.index){
    anc.index <- phy$edge[i,1]
    dec.index <- phy$edge[i,2]
    # init and final
    init <- CharHistory[anc.index]
    final <- CharHistory[dec.index]
    total.bl <- phy$edge.length[i]
    SubstHistory[[i]] <- simBranchSubstHistory(init, final, total.bl, model, d.rates, max.attempt = max.attempt)
  }
  return(SubstHistory)
}

# simulate a substitution history given the simulations of ancestral states
simSubstHistory <- function(phy, tip.states, states, model, nSim, nCores, max.attempt = 1000, vector.form=FALSE, return.char = FALSE){
  # set-up
  cladewise.index <- reorder.phylo(phy, "cladewise", index.only = TRUE)
  # a potential speedup is to calculate all Pij (bollback eq.3) for all branches first
  Pij <- array(0, c(dim(model)[1], dim(model)[2], length(phy$edge.length)))
  # plus one because the root has no edge
  Pj <- array(0, c(dim(model)[1], dim(model)[2], length(phy$edge.length)+1))
  for(i in 1:length(phy$edge.length)){
    Pij[,,i] <- expm(model * phy$edge.length[i])
  }
  # then multiply Pij by l_sigma-1_j
  l_sigma_j <- rbind(tip.states, states)
  dec.index <- phy$edge[,2]
  count <- 1
  for(i in dec.index){
    Pj[,,count] <- sweep(Pij[,,count], MARGIN = 2, l_sigma_j[i,], '*')
    count <- count + 1
  }
  # simulate a character history
  CharHistories <- lapply(1:nSim, function(x) simCharHistory(phy=phy, Pj=Pj, root=states[1,], model=model, vector.form=vector.form))
  if(return.char){
    return(CharHistories)
  }
  obj <- mclapply(CharHistories, function(x) simSingleSubstHistory(cladewise.index, x, phy, model, max.attempt), mc.cores = nCores)
  return(obj)
}

# convert a substitution history into a mapped edge
convertSubHistoryToEdge <- function(phy, map){
  Traits <- as.numeric(unique(names(unlist(map))))
  RowNames <- apply(phy$edge, 1, function(x) paste(x[1], x[2], sep = ","))
  obj <- do.call(rbind, lapply(map, function(x) sapply(Traits, function(y) sum(x[names(x) == y]))))
  rownames(obj) <- RowNames
  colnames(obj) <- Traits
  return(obj)
}

# get the conditional likelihoods of particular nodes
getConditionalNodeLik <- function(tree, data, model, rate.cat, root.p, parsimony=FALSE, collapse=TRUE, model.set.final=NULL){
  phy <- reorder(tree, "pruningwise")
  nb.node <- phy$Nnode
  nb.tip <- length(phy$tip.label)
  # get the liks table
  if(is.null(model.set.final)){
    # process the data to match the liks table
    CorData <- corProcessData(data, collapse=collapse)
    nObs <- length(CorData$ObservedTraits)
    model.set.final <- rate.cat.set.corHMM.JDB(phy=phy,data=data, rate.cat=rate.cat, 
      ntraits = nObs, model = "ER", collapse=collapse, rate.mat=model)
  }
  if(parsimony==TRUE){
    model <- model/1000
  }
  liks <- model.set.final$liks
  anc <- unique(phy$edge[,1])
  for (i in seq(from = 1, length.out = nb.node)) {
    #the ancestral node at row i is called focal
    focal <- anc[i]
    #Get descendant information of focal
    desRows <- which(phy$edge[,1]==focal)
    desNodes <- phy$edge[desRows,2]
    v <- 1
    #Loops through all descendants of focal (how we deal with polytomies):
    for (desIndex in sequence(length(desRows))){
      v <- v*expm(model * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks[desNodes[desIndex],]
    }
    #Divide each likelihood by the sum to obtain probabilities:
    liks[focal, ] <- v/sum(v)
  }
  input.root.p <- root.p
  if(is.numeric(input.root.p)){
    root.p <- input.root.p/sum(input.root.p)
  }
  if(is.character(input.root.p)){
      if(input.root.p == "yang"){
        root.p <- Null(model)
        root.p <- c(root.p/sum(root.p))
        }
      if(input.root.p == "maddfitz"){
        root.p <- liks[focal, ]
        }
      if(input.root.p == "flat"){
        root.p <- rep(1/dim(model)[1], dim(model))
        }
  }
  if(is.null(input.root.p)){
    root.p <- rep(1/dim(model)[1], dim(model))
  }
  liks[focal, ] <-  liks[focal, ] * root.p
  return(list(tip.states = liks[1:nb.tip,],
              node.states = liks[(nb.tip+1):(nb.node+nb.tip),]))
}

# simulate a Q along a phylogeny
simMarkov <- function(phy, Q, root.freqs){

  #Randomly choose starting state at root using the root.values as the probability:
  root.value <- sample.int(dim(Q)[2], 1, FALSE, prob=root.freqs/sum(root.freqs))
  #Reorder the phy:
  phy <- reorder(phy, "postorder")
  ntips <- length(phy$tip.label)
  N <- dim(phy$edge)[1]
  ROOT <- ntips + 1 #perhaps use an accessor to get the root node id

  #Generate vector that contains the simulated states:
  CharacterHistory <- integer(ntips + phy$Nnode)
  CharacterHistory[ROOT] <- as.integer(root.value)
  anc <- phy$edge[, 1]
  des <- phy$edge[, 2]
  edge.length <- phy$edge.length
  diag(Q) = 0
  diag(Q) = -rowSums(Q)

  # setting up the alternative Q matrix at the node of interest
  # if(!any(is.na(Q2))){
  #   diag(Q2) = 0
  #   diag(Q2) = -rowSums(Q2)
  # }
  # if(!is.na(NoI)){
  #   NewQDesc <- getDescendants(phy, NoI)
  # }

  #standard simulation protocol
  # if(any(is.na(Q2)) | is.na(NoI)){
  for (i in N:1) {
    p <- expm(Q * edge.length[i], method="Ward77")[CharacterHistory[anc[i]], ]
    CharacterHistory[des[i]] <- sample.int(dim(Q)[2], size = 1, FALSE, prob = p)
  }
  # }

  # simulating a clade under a different (Q2) evolutionary model
  # if(!any(is.na(Q2)) & !is.na(NoI)){
  #   for (i in N:1) {
  #     if(anc[i] %in% NewQDesc){
  #       p <- expm(Q2 * edge.length[i], method="Ward77")[CharacterHistory[anc[i]], ]
  #       CharacterHistory[des[i]] <- sample.int(dim(Q2)[2], size = 1, FALSE, prob = p)
  #     }else{
  #       p <- expm(Q * edge.length[i], method="Ward77")[CharacterHistory[anc[i]], ]
  #       CharacterHistory[des[i]] <- sample.int(dim(Q)[2], size = 1, FALSE, prob = p)
  #     }
  #   }
  # }

  TipStates <-  CharacterHistory[1:ntips]
  names(TipStates) <- phy$tip.label
  NodeStates <- CharacterHistory[ROOT:(N+1)]
  names(NodeStates) <- ROOT:(N+1)

  res <- list(TipStates = TipStates, NodeStates = NodeStates)
  return(res)
}

# get the probability of a particular painting of a branch
getSimmapBranchProb <- function(branch, Q){
  P.Branch <- numeric(length(branch))
  count <- 1
  # 1. start rootward and determine if transition
  while(length(branch) > 1){
    # 2. if transition calculate probabality
    from <- as.numeric(names(branch)[1])
    to <- as.numeric(names(branch)[2])
    qii <- -diag(Q)[from]
    qij <- Q[from, to]
    t <- branch[1]
    P.waiting.time <- qii * exp(-qii * t)
    P.trans <- qij/qii
    P.Branch[count] <- P.waiting.time * P.trans
  # 3. start at new starting place and determine if transition
    count <- count + 1
    branch <- branch[-1]
  }
  # 4. if transition repeat step 1.
  # 5. else no transition calculate probabliity
  from <- as.numeric(names(branch)[1])
  qii <- -diag(Q)[from]
  t <- branch[1]
  P.Branch[count] <- 1 - (qii * exp(-qii * t))
  return(prod((P.Branch)))
}

# get the likelihood/ of a particular simmap
getSimmapLik <- function(simmap, Q){
  maps <- simmap$maps
  lik <- prod(unlist(lapply(maps, function(x) getSimmapBranchProb(x, Q))))
  return(lik)
}

FloydWalshAlg <- function(model, init, final){
  nStates <- dim(model)[1]
  Dist <- matrix(Inf, nStates, nStates)
  Next <- matrix(NA, nStates, nStates)
  for(V_index in sequence(nStates)){
    Dist[V_index, V_index] <- 0
    Next[V_index, V_index] <- V_index
  }
  for(V_index in sequence(nStates)){
    To <- which(model[V_index, ] > 0)
    Dist[V_index, To] <- 1
    Next[V_index, To] <- To
  }
  for(k in sequence(nStates)){
    for(i in sequence(nStates)){
      for(j in sequence(nStates)){
        if(Dist[i,j] > Dist[i,k] + Dist[k,j]){
          Dist[i,j] = Dist[i,k] + Dist[k,j]
          Next[i,j] = Next[i,k]
        }
      }
    }
  }
  if(is.na(Next[init, final])){
    stop("Impossible transition on branch detected...")
  }else{
    path = init
    u = init
    v = final
    while(u != v){
      u <- Next[u, v]
      path <- c(path, u)
    }
  }
  return(path)
}

summarize_single_simmap <- function(simmap){
  transition_index <- which(unlist(lapply(simmap$maps, length)) > 1)
  bt <- branching.times(simmap)
  transition_df <- data.frame(branch_id = NA, transition = NA, age = NA)
  for(i in seq_along(transition_index)){
    single_map <- simmap$maps[[transition_index[i]]]
    for(j in 2:length(single_map)){
      transition <- paste(names(single_map)[j-1], names(single_map)[j], sep = ",")
      branch_index <- simmap$edge[transition_index[i],1] - Ntip(simmap)
      timing <- bt[branch_index] - sum(single_map[1:j-1])
      tmp <- data.frame(branch_id = transition_index[i], transition = transition, age = timing)
      transition_df <- rbind(transition_df, tmp)
    }
  }
  transition_df <- transition_df[-1,]
  rownames(transition_df) <- NULL
  return(transition_df)
}

summarize_transition_stats <- function(simmap_summaries) {
  num_sims <- length(simmap_summaries)
  all_transitions_df <- do.call(rbind, simmap_summaries)
  unique_transitions <- unique(all_transitions_df$transition)
  counts_per_sim <- lapply(simmap_summaries, function(df) {
    table(factor(df$transition, levels = unique_transitions))
  })
  counts_matrix <- do.call(rbind, counts_per_sim)
  avg_counts <- colMeans(counts_matrix)
  sd_counts <- apply(counts_matrix, 2, sd)
  avg_ages <- aggregate(age ~ transition, data = all_transitions_df, FUN = mean)
  sd_ages <- aggregate(age ~ transition, data = all_transitions_df, FUN = sd)
  summary_df <- data.frame(
    transition = names(avg_counts),
    avg_count = avg_counts,
    sd_count = sd_counts
  )
  summary_df <- merge(summary_df, avg_ages, by = "transition")
  summary_df <- merge(summary_df, sd_ages, by = "transition", suffixes = c("_avg", "_sd"))
  names(summary_df)[names(summary_df) == "age_avg"] <- "avg_age"
  names(summary_df)[names(summary_df) == "age_sd"] <- "sd_age"
  
  first_transition_times <- do.call(
    rbind,
    lapply(seq_along(simmap_summaries), function(i) {
      df <- simmap_summaries[[i]]
      first_df <- aggregate(age ~ transition, data = df, FUN = max)
      first_df$sim <- i
      first_df
    })
  )
  
  avg_first_age <- aggregate(age ~ transition, data = first_transition_times, mean)
  sd_first_age  <- aggregate(age ~ transition, data = first_transition_times, sd)
  
  names(avg_first_age)[2] <- "avg_first_age"
  names(sd_first_age)[2]  <- "sd_first_age"
  
  summary_df <- merge(summary_df, avg_first_age, by = "transition", all.x = TRUE)
  summary_df <- merge(summary_df, sd_first_age,  by = "transition", all.x = TRUE)
  
  rownames(summary_df) <- NULL
  return(summary_df)
}

plot_transition_summary <- function(simmap_summaries, cols = NULL) {
  all_transitions_df <- do.call(rbind, simmap_summaries)
  if (nrow(all_transitions_df) == 0) {
    warning("The provided simmap summaries are empty. No plots will be generated.")
    return(invisible(NULL))
  }
  
  unique_transitions <- sort(unique(all_transitions_df$transition))
  counts_per_sim <- lapply(simmap_summaries, function(df) {
    table(factor(df$transition, levels = unique_transitions))
  })
  counts_matrix <- do.call(rbind, counts_per_sim)
  counts_df <- as.data.frame(counts_matrix)
  
  if (is.null(cols)) {
    cols <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
      "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC")
  }
  plot_cols <- rep(cols, length.out = length(unique_transitions))
  transparent_cols <- adjustcolor(plot_cols, alpha.f = 0.3)
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(1, 2), oma = c(0, 0, 3, 0))
  
  par(mar = c(5, 4, 4, 1))
  boxplot(counts_df,
    las = 2,
    col = plot_cols,
    main = "Distribution of Transition Counts",
    ylab = "Count per Simulation",
    xlab = "",
    border = "gray40",
    outline = FALSE)
  
  stripchart(counts_df,
    vertical = TRUE,
    method = "jitter",
    add = TRUE,
    pch = 19,
    col = transparent_cols)
  
  par(mar = c(5, 5, 4, 2))
  n_transitions <- length(unique_transitions)
  overlap_factor <- 0.8
  xlim_age <- range(all_transitions_df$age, na.rm = TRUE)
  ylim_age <- c(0.5, n_transitions + overlap_factor)
  
  plot(0, type = 'n',
    xlim = xlim_age,
    ylim = ylim_age,
    main = "Distribution of Transition Ages",
    xlab = "Age (Time from present)",
    ylab = "",
    yaxt = "n")
  
  axis(2, at = (1:n_transitions) + (overlap_factor / 3), labels = unique_transitions, las = 1)
  
  for (i in 1:n_transitions) {
    current_transition <- unique_transitions[i]
    ages_current <- all_transitions_df$age[all_transitions_df$transition == current_transition]
    
    if (length(ages_current) > 1) {
      d <- density(ages_current, from = xlim_age[1], to = xlim_age[2])
      scaled_y <- d$y / max(d$y, na.rm = TRUE) * overlap_factor
      polygon(c(d$x, rev(d$x)), c(scaled_y + i, rep(i, length(d$y))),
        col = plot_cols[i],
        border = "gray40")
    }
  }
  
  mtext("Summary of Transition Counts and Ages", outer = TRUE, cex = 1.5, line = 0.5)
}
# simmap <- makeSimmap(tree=phy, tip.states=tip.states, states=states, model=model, 
#                      nSim=10, nCores=1)
# 
# liks <- unlist(lapply(simmap, function(x) getSimmapLik(x, Q)))
# log(sum(exp(liks)))

# 
# library(castor)
# library(corHMM)
# 
# data(primates)
# phy <- primates[[1]]
# phy <- multi2di(phy)
# data <- primates[[2]][,-3]
# 
# ##run corhmm
# MK <- corHMM(phy, data, 1)
# 
# ##get simmap from corhmm solution
# model <- MK$solution
# simmap <- makeSimmap(tree=phy, data=data, model=model, rate.cat=1, nSim=1, nCores=1)[[1]]
# 
# ltt_dat <- count_lineages_through_time(phy, 100)
# 
# get_edge_span <- function(simmap){
#   bt <- branching.times(simmap)
#   time_frame <- simmap$edge
#   for(i in 1:length(bt)){
#     index <- simmap$edge[,1] == names(bt)[i]
#     time_frame[index, 1] <- bt[i]
#     time_frame[index, 2] <- bt[i] - simmap$edge.length[index]
#   }
#   return(time_frame)
# }
# 
# get_char_at_time <- function(time, simmap){
#   edge_span <- round(get_edge_span(simmap), 4)
#   edge_index <- which((time <= edge_span[,1]) & (time >= edge_span[,2]))
#   bt <- branching.times(simmap)
#   char_at_time <- c()
#   for(i in edge_index){
#     focal_map <- simmap$maps[[i]]
#     if(length(focal_map) == 1){
#       char <- names(focal_map)
#     }else{
#       starting_age <- bt[names(bt) == simmap$edge[i,1]]
#       time_span <- round(starting_age - cumsum(focal_map), 5)
#       possible_index <- which(time >= time_span)
#       char <- names(focal_map)[possible_index[1]]
#     }
#     char_at_time <- c(char_at_time, char)
#   }
#   char_at_time <- c(char_at_time, colnames(simmap$mapped.edge))
#   out <- table(char_at_time)
#   out <- out - 1
#   return(out)
# }
# 
# ages <- round(ltt_dat$ages, 4)
# test <- sapply(ages, function(x) get_char_at_time(x, simmap))
# plot_table <- data.frame(ages = ages, t(test))
# 
# library(ggplot2)
# 
# ggplot(plot_table, aes(x = rev(ages))) +
#   geom_line(aes(y = log(X0), color = "X0"), size = 1) +
#   geom_line(aes(y = log(X1), color = "X1"), size = 1) +
#   scale_color_manual(values = c("X0" = "blue", "X1" = "red")) +
#   labs(x = "Ages", y = "Values") +
#   theme_minimal()

