simCharacterHistory <- function(phy, Q, root.freqs, Q2 = NA, NoI = NA){
  
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
  #return(CharacterHistory)
}

simProbalisticHistory <- function(phy, Q, root.freqs, Q2 = NA, NoI = NA){
  #Randomly choose starting state at root using the root.values as the probability:
  root.value <- root.freqs
  #Reorder the phy:
  phy <- reorder(phy, "postorder")
  ntips <- length(phy$tip.label)
  N <- dim(phy$edge)[1]
  ROOT <- ntips + 1 #perhaps use an accessor to get the root node id
  #Generate vector that contains the simulated states:
  CharacterHistory <- matrix(0, ntips + phy$Nnode, dim(Q)[1])
  CharacterHistory[ROOT,] <- as.integer(root.value)
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
      p <- CharacterHistory[anc[i],] %*% expm(Q * edge.length[i], method="Ward77")
      CharacterHistory[des[i],] <- c(p)/sum(p)
    }
  }
  
  # simulating a clade under a different (Q2) evolutionary model
  if(!any(is.na(Q2)) & !is.na(NoI)){
    for (i in N:1) {
      if(anc[i] %in% NewQDesc){
        p <- CharacterHistory[anc[i],] %*% expm(Q2 * edge.length[i], method="Ward77")
        CharacterHistory[des[i],][sample.int(dim(Q2)[2], size = 1, FALSE, prob = c(p)/sum(p))] <- 1
      }else{
        p <- CharacterHistory[anc[i],] %*% expm(Q * edge.length[i], method="Ward77")
        CharacterHistory[des[i],][sample.int(dim(Q)[2], size = 1, FALSE, prob = c(p)/sum(p))] <- 1
      }
    }
  }
  TipStates <-  CharacterHistory[1:ntips]
  names(TipStates) <- 1:ntips
  NodeStates <- CharacterHistory[ROOT:(N+1)]
  names(NodeStates) <- ROOT:(N+1)
  
  res <- list(TipStates = TipStates, NodeStates = NodeStates)
  return(res)
  #return(CharacterHistory)
}

