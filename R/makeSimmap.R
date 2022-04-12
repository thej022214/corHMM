# exported function for use
makeSimmap <- function(tree, data, model, rate.cat, root.p="yang", nSim=1, nCores=1, fix.node=NULL, fix.state=NULL, parsimony = FALSE, max.attempt = 1000, collapse=TRUE){
  model[is.na(model)] <- 0
  diag(model) <- 0
  diag(model) <- -rowSums(model)
  conditional.lik <- getConditionalNodeLik(tree, data, model, rate.cat, root.p, parsimony = parsimony, collapse=collapse)
  # if(!is.null(tip.probs)){
  #   for(i in 1:dim(conditional.lik$tip.states)[1]){
  #     nCol <- dim(tip.probs)[2]
  #     tmp <- vector("numeric", nCol)
  #     tmp[sample(1:nCol, 1, prob = tip.probs[i,])] <- 1
  #     conditional.lik$tip.states[i, ] <- tmp
  #   }
  # }
  if(length(fix.node) != length(fix.state)){
    stop("The number of nodes supplied to be fixed does not match the number of states provided.",.call=FALSE)
  }
  if(max(fix.state) > dim(model)[1]){
    stop("One of the states being fixed does not exist in this model.")
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
  legend <- getStateMat4Dat(data, collapse = collapse)$legend
  if(rate.cat > 1){
    StateNames <- paste("(", rep(legend, rate.cat), ",", rep(paste("R", 1:rate.cat, sep = ""), each = length(legend)), ")", sep = "")
    names(StateNames) <- 1:length(StateNames)
  }else{
    StateNames <- legend
  }
  for(i in 1:nSim){
    tree.simmap <- tree
    tree.simmap$maps <- maps[[i]]
    tree.simmap$maps <- lapply(maps[[i]], function(x) correctMapName(x, StateNames))
    tree.simmap$mapped.edge <- mapped.edge[[i]]
    colnames(tree.simmap$mapped.edge) <- StateNames
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
  names(map_element) <- state_names[match(as.numeric(names(map_element)), as.numeric(names(state_names)))]
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
  # combine data in order of the edge matrix
  state.sample <- matrix(0, dim(Pj)[3], dim(Pj)[2])
  # sample the root
  anc.index <- phy$edge[cladewise.index[1],1]
  state.sample[anc.index,] <- c(rmultinom(1, 1, root))
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
    state.sample[dec.index,] <- c(rmultinom(1, 1, p))
  }
  if(vector.form == FALSE){
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
simSubstHistory <- function(phy, tip.states, states, model, nSim, nCores, max.attempt = 1000){
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
  CharHistories <- lapply(1:nSim, function(x) simCharHistory(phy=phy, Pj=Pj, root=states[1,], model=model))
  obj <- mclapply(CharHistories, function(x) simSingleSubstHistory(cladewise.index, x, phy, model, max.attempt), mc.cores = nCores)
  return(obj)
}

# convert a substitution history into a mapped edge
convertSubHistoryToEdge <- function(phy, map){
  Traits <- levels(as.factor(unique(names(unlist(map)))))
  RowNames <- apply(phy$edge, 1, function(x) paste(x[1], x[2], sep = ","))
  obj <- do.call(rbind, lapply(map, function(x) sapply(Traits, function(y) sum(x[names(x) == y]))))
  rownames(obj) <- RowNames
  return(obj)
}

# get the conditional likelihoods of particular nodes
getConditionalNodeLik <- function(tree, data, model, rate.cat, root.p, parsimony=FALSE, collapse=TRUE){
  phy <- reorder(tree, "pruningwise")
  nb.node <- phy$Nnode
  nb.tip <- length(phy$tip.label)
  # process the data to match the liks table
  CorData <- corProcessData(data, collapse=collapse)
  nObs <- length(CorData$ObservedTraits)
  # get the liks table
  model.set.final <- rate.cat.set.corHMM.JDB(phy=phy,data=data, rate.cat=rate.cat, ntraits = nObs, model = "ER", collapse=collapse, rate.mat=model)
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
  if(!is.null(input.root.p)){
    if(!is.character(input.root.p)){
      root.p <- input.root.p/sum(input.root.p)
    }
  }
  if(is.null(input.root.p) | input.root.p == "flat"){
    root.p <- rep(1/dim(model)[1], dim(model))
  }
  if(input.root.p == "yang"){
    root.p <- Null(model)
    root.p <- c(root.p/sum(root.p))
  }
  if(input.root.p == "maddfitz"){
    root.p <- liks[focal, ]
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

# The Floyd–Warshall algorithm (also known as Floyd's algorithm, the Roy–Warshall algorithm, the Roy–Floyd algorithm, or the WFI algorithm) is an algorithm for finding shortest paths in a directed weighted graph with positive or negative edge weights (but with no negative cycles).
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

# simmap <- makeSimmap(tree=phy, tip.states=tip.states, states=states, model=model, 
#                      nSim=10, nCores=1)
# 
# liks <- unlist(lapply(simmap, function(x) getSimmapLik(x, Q)))
# log(sum(exp(liks)))
