# simulate ancestral states at each internal node
simSingleCharHistory<- function(phy, model, cladewise.index, state.probability, state.sample){
  # sample the root state
  anc.index <- phy$edge[cladewise.index[1],1]
  state.sample[anc.index,] <- t(rmultinom(1, 1, state.probability[anc.index,]))
  
  # sample the nodes
  for(i in cladewise.index){
    # A = the probability of transition FROM node i (ancestor)
    anc.index <- phy$edge[i,1]
    dec.index <- phy$edge[i,2]
    A <- state.sample[anc.index,] %*% expm(model * phy$edge.length[i])
    # B = the conditional probability of node j (descendent)
    B <- state.probability[dec.index,]
    # The decendant is sampled from p/sum(p), where p = A * B
    p <- A * B
    state.sample[dec.index,] <- t(rmultinom(1, 1, p))
  }
  return(state.sample)
}


simCharHistory <- function(phy, tip.states, states, model, vector.form = FALSE){
  # warnings
  if(!all(round(rowSums(tip.states)) == 1)){
    return(cat("\nWarning: Tip states must be reconstructed to run Simmap.\n"))
  }
  
  # organize
  nTip <- Ntip(phy)
  cladewise.index <- reorder.phylo(phy, "cladewise", index.only = TRUE)
  # combine data in order of the edge matrix
  state.probability <- rbind(tip.states, states)
  state.sample <- matrix(0, dim(state.probability)[1], dim(state.probability)[2])
  
  # single sim function
  state.samples <- simSingleCharHistory(phy, model, cladewise.index, state.probability, state.sample)
  if(vector.form == FALSE){
    state.samples <- apply(state.samples, 1, function(y) which(y == 1))
  }
  return(state.samples)
}

simBranchSubstHistory <- function(init, final, total.bl, model){
  current.bl <- c()
  current.state <- init
  restart <- TRUE
  while(restart == TRUE){
    # draw a rate and waiting time based on init
    rate <- -diag(model)[current.state]
    waiting.time <- rexp(1, rate)
    current.bl <- c(current.bl, waiting.time)
    names(current.bl)[length(current.bl)] <- current.state
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

simSingleSubstHistory <- function(cladewise.index, CharHistory, phy, model){
  SubstHistory <- vector("list", length(cladewise.index))
  for(i in cladewise.index){
    anc.index <- phy$edge[i,1]
    dec.index <- phy$edge[i,2]
    # init and final
    init <- CharHistory[anc.index]
    final <- CharHistory[dec.index]
    total.bl <- phy$edge.length[i]
    SubstHistory[[i]] <- simBranchSubstHistory(init, final, total.bl, model)
  }
  return(SubstHistory)
}

# simulate a substitution history given the simulations of ancestral states
simSubstHistory <- function(phy, tip.states, states, model){
  # set-up
  cladewise.index <- reorder.phylo(phy, "cladewise", index.only = TRUE)
  # simulate a character history
  CharHistories <- simCharHistory(phy=phy, tip.states=tip.states, states=states, model=model)
  obj <- simSingleSubstHistory(cladewise.index, CharHistories, phy, model)
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

# exported function for use
makeSimmap <- function(tree, tip.states, states, model, nSim=1, nCores=1){
  maps <- mclapply(1:nSim, function(x) simSubstHistory(tree, tip.states, states, model), mc.cores = nCores)
  mapped.edge <- lapply(maps, function(x) convertSubHistoryToEdge(tree, x))
  obj <- vector("list", nSim)
  for(i in 1:nSim){
    tree.simmap <- tree
    tree.simmap$maps <- maps[[i]]
    tree.simmap$mapped.edge <- mapped.edge[[i]]
    tree.simmap$Q <- model
    attr(tree.simmap, "map.order") <- "right-to-left"
    if (!inherits(tree.simmap, "simmap")) 
      class(tree.simmap) <- c("simmap", setdiff(class(tree.simmap), "simmap"))
    obj[[i]] <- tree.simmap
  }
  return(obj)
}






##### a set of secret functions being used to test why number of substitutions is underestimated
makeSimmapTipsy <- function(tree, tip.states, states, model, nSim=1, nCores=1){
  # the substitution history needs to allow any tip history
  maps <- mclapply(1:nSim, function(x) simSubstHistoryTipsy(tree, tip.states, states, model), mc.cores = nCores)
  mapped.edge <- lapply(maps, function(x) convertSubHistoryToEdge(tree, x))
  obj <- vector("list", nSim)
  for(i in 1:nSim){
    tree.simmap <- tree
    tree.simmap$maps <- maps[[i]]
    tree.simmap$mapped.edge <- mapped.edge[[i]]
    tree.simmap$Q <- model
    attr(tree.simmap, "map.order") <- "right-to-left"
    if (!inherits(tree.simmap, "simmap")) 
      class(tree.simmap) <- c("simmap", setdiff(class(tree.simmap), "simmap"))
    obj[[i]] <- tree.simmap
  }
  return(obj)
}

simSubstHistoryTipsy <- function(phy, tip.states, states, model){
  # set-up
  cladewise.index <- reorder.phylo(phy, "cladewise", index.only = TRUE)
  # simulate a character history
  StateMax <- dim(tip.states)[2]
  tipsy.states <- matrix(0, dim(tip.states)[1], dim(tip.states)[2])
  new.tips <- round(runif(length(phy$tip.label), 1, StateMax))
  for(i in 1:length(new.tips)){
    tipsy.states[i,new.tips[i]] <- 1
  }
  CharHistories <- simCharHistory(phy=phy, tip.states=tipsy.states, states=states, model=model)
  obj <- simSingleSubstHistory(cladewise.index, CharHistories, phy, model)
  return(obj)
}

simMarkov <- function(phy, Q, root.freqs, Q2 = NA, NoI = NA){

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
  if(!any(is.na(Q2))){
    diag(Q2) = 0
    diag(Q2) = -rowSums(Q2)
  }
  if(!is.na(NoI)){
    NewQDesc <- getDescendants(phy, NoI)
  }

  #standard simulation protocol
  if(any(is.na(Q2)) | is.na(NoI)){
    for (i in N:1) {
      p <- expm(Q * edge.length[i], method="Ward77")[CharacterHistory[anc[i]], ]
      CharacterHistory[des[i]] <- sample.int(dim(Q)[2], size = 1, FALSE, prob = p)
    }
  }

  # simulating a clade under a different (Q2) evolutionary model
  if(!any(is.na(Q2)) & !is.na(NoI)){
    for (i in N:1) {
      if(anc[i] %in% NewQDesc){
        p <- expm(Q2 * edge.length[i], method="Ward77")[CharacterHistory[anc[i]], ]
        CharacterHistory[des[i]] <- sample.int(dim(Q2)[2], size = 1, FALSE, prob = p)
      }else{
        p <- expm(Q * edge.length[i], method="Ward77")[CharacterHistory[anc[i]], ]
        CharacterHistory[des[i]] <- sample.int(dim(Q)[2], size = 1, FALSE, prob = p)
      }
    }
  }

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
  return(sum(log(P.Branch)))
}

# get the likelihood/ of a particular simmap
getSimmapLik <- function(simmap, Q){
  maps <- simmap$maps
  lik <- sum(unlist(lapply(maps, function(x) getSimmapBranchProb(x, Q))))
  return(lik)
}

# simmap <- makeSimmap(tree=phy, tip.states=tip.states, states=states, model=model, 
#                      nSim=10, nCores=1)
# 
# liks <- unlist(lapply(simmap, function(x) getSimmapLik(x, Q)))
# log(sum(exp(liks)))
