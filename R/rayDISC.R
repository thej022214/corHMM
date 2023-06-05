#############################
#EVOLUTION OF DISCRETE TRAITS, ALLOWING POLYMORPHIC AND MISSING STATES
#library(expm)
#library(phangorn)
#############################

#written by Jeremy M. Beaulieu & Jeffrey C. Oliver

rayDISC<-function(phy, data, ...){

  return(corHMM(phy = phy,
                data = data,
                rate.cat = 1, 
                ...))  
}

##Keeping this because other functions in other packages use this (i.e., selac):
dev.raydisc <- function(p, phy, liks, Q, rate, root.p, lewis.asc.bias){
  
  p.new <- exp(p)
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  TIPS <- 1:nb.tip
  comp <- numeric(nb.tip + nb.node)
  
  #Obtain an object of all the unique ancestors
  anc <- unique(phy$edge[,1])
  #This bit is to allow packages like "selac" the ability to deal with this function directly:
  if(is.null(rate)){
    Q=Q
  }else{
    if (any(is.nan(p.new)) || any(is.infinite(p.new))) return(1000000)
    Q[] <- c(p.new, 0)[rate]
    diag(Q) <- -rowSums(Q)
  }
  
  for (i  in seq(from = 1, length.out = nb.node)) {
    #the ancestral node at row i is called focal
    focal <- anc[i]
    #Get descendant information of focal
    desRows<-which(phy$edge[,1]==focal)
    desNodes<-phy$edge[desRows,2]
    v <- 1
    for (desIndex in sequence(length(desRows))){
      v <- v * expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks[desNodes[desIndex],]
    }
    comp[focal] <- sum(v)
    liks[focal, ] <- v/comp[focal]
  }
  
  #Specifies the root:
  root <- nb.tip + 1L
  #If any of the logs have NAs restart search:
  if (is.na(sum(log(comp[-TIPS])))){return(1000000)}
  else{
    equil.root <- NULL
    for(i in 1:ncol(Q)){
      posrows <- which(Q[,i] >= 0)
      rowsum <- sum(Q[posrows,i])
      poscols <- which(Q[i,] >= 0)
      colsum <- sum(Q[i,poscols])
      equil.root <- c(equil.root,rowsum/(rowsum+colsum))
    }
    if (is.null(root.p)){
      flat.root = equil.root
      k.rates <- 1/length(which(!is.na(equil.root)))
      flat.root[!is.na(flat.root)] = k.rates
      flat.root[is.na(flat.root)] = 0
      root.p <- flat.root
      loglik <- sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,]))))
    }else{
      if(is.character(root.p)){
        # root.p==yang will fix root probabilities based on the inferred rates: q10/(q01+q10)
        if(root.p == "yang"){
          root.p <- Null(Q)
          root.p <- c(root.p/sum(root.p))
          loglik <- sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,]))))
          if(is.infinite(loglik)){
            return(1000000)
          }
        }else{
          # root.p==maddfitz will fix root probabilities according to FitzJohn et al 2009 Eq. 10:
          root.p = liks[root,] / sum(liks[root,])
          loglik <- sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,]))))
        }
      }
      # root.p!==NULL will fix root probabilities based on user supplied vector:
      else{
        loglik <- sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,]))))
        if(is.infinite(loglik)){
          return(1000000)
        }
      }
    }
  }
  
  if(lewis.asc.bias == TRUE){
    dummy.liks.vec <- numeric(dim(Q)[1])
    for(state.index in 1:dim(Q)[1]){
      dummy.liks.vec[state.index] <- CalculateLewisLikelihood(p=p.new, phy=phy, liks=liks, Q=Q, rate=rate, root.p=root.p, state.num=state.index)
    }
    loglik <- loglik - log(sum(root.p * (1 - exp(dummy.liks.vec))))
  }
  -loglik
}



CalculateLewisLikelihood <- function(p, phy, liks, Q, rate, root.p, state.num=1){
  
  p.new <- p
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  TIPS <- 1:nb.tip
  comp <- numeric(nb.tip + nb.node)
  
  #Obtain an object of all the unique ancestors
  anc <- unique(phy$edge[,1])
  #This bit is to allow packages like "selac" the ability to deal with this function directly:
  if(is.null(rate)){
    Q=Q
  }else{
    if (any(is.nan(p.new)) || any(is.infinite(p.new))) return(1000000)
    Q[] <- c(p.new, 0)[rate]
    diag(Q) <- -rowSums(Q)
  }
  
  liks.dummy <- liks
  liks.dummy[TIPS,] = 0
  liks.dummy[TIPS,state.num] = 1
  comp.dummy <- comp
  for (i  in seq(from = 1, length.out = nb.node)) {
    #the ancestral node at row i is called focal
    focal <- anc[i]
    #Get descendant information of focal
    desRows <- which(phy$edge[,1]==focal)
    desNodes <- phy$edge[desRows,2]
    v.dummy <- 1
    for(desIndex in sequence(length(desRows))){
      v.dummy <- v.dummy * expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks.dummy[desNodes[desIndex],]
    }
    comp.dummy[focal] <- sum(v.dummy)
    liks.dummy[focal, ] <- v.dummy/comp.dummy[focal]
  }
  #Specifies the root:
  root <- nb.tip + 1L
  #If any of the logs have NAs restart search:
  if(is.na(sum(log(comp[-TIPS])))){return(1000000)}
  else{
    equil.root <- NULL
    for(i in 1:ncol(Q)){
      posrows <- which(Q[,i] >= 0)
      rowsum <- sum(Q[posrows,i])
      poscols <- which(Q[i,] >= 0)
      colsum <- sum(Q[i,poscols])
      equil.root <- c(equil.root,rowsum/(rowsum+colsum))
    }
    if (is.null(root.p)){
      flat.root = equil.root
      k.rates <- 1/length(which(!is.na(equil.root)))
      flat.root[!is.na(flat.root)] = k.rates
      flat.root[is.na(flat.root)] = 0
      loglik <- (sum(log(comp.dummy[-TIPS])) + log(sum(exp(log(flat.root)+log(liks.dummy[root,])))))
    }else{
      if(is.character(root.p)){
        # root.p==yang will fix root probabilities based on the inferred rates: q10/(q01+q10)
        if(root.p == "yang"){
          root.p <- Null(Q)
          root.p <- c(root.p/sum(root.p))
          loglik <- (sum(log(comp.dummy[-TIPS])) + log(sum(exp(log(root.p)+log(liks.dummy[root,])))))
          if(is.infinite(loglik)){
            return(1000000)
          }
        }else{
          # root.p==maddfitz will fix root probabilities according to FitzJohn et al 2009 Eq. 10:
          root.p = liks.dummy[root,] / sum(liks.dummy[root,])
          loglik <- -(sum(log(comp.dummy[-TIPS])) + log(sum(exp(log(root.p)+log(liks.dummy[root,])))))
        }
      }else{# root.p!==NULL will fix root probabilities based on user supplied vector:
        loglik <- (sum(log(comp.dummy[-TIPS])) + log(sum(exp(log(root.p)+log(liks.dummy[root,])))))
        if(is.infinite(loglik)){
          return(1000000)
        }
      }
    }
  }
  loglik
  
}


#########################
#    match.tree.data    #
#########################
# Compares a tree and data to make sure they include the same taxa
# Taxa which are in the tree, but not the data matrix, are added to the matrix and coded as missing data.
# Any taxa in the data matrix which are not in the tree are removed from the matrix
# The function returns an object with three parts:
#	$phy: the tree
#	$data: the matrix, omitting taxa not in tree and taxa that were present in the tree but not in the matrix
#	$message.data: a brief message explaining modifications (if any) to the data
#	$message.tree: a brief message explaining modificatoins (if any) to the tree
match.tree.data <- function(phy, data){
  matchobj <- NULL
  matchobj$phy <- phy
  matchobj$data <- data
  matchobj$message.data <- NULL
  matchobj$message.tree <- NULL
  # First look at data matrix to see if each taxon in matrix is also in tree
  missing.fromtree <- NULL
  for(datarow in 1:length(data[,1])){
    if(is.na(match(data[datarow,1],phy$tip.label))){
      missing.fromtree <- c(missing.fromtree,datarow)
    }
  }
  if(length(missing.fromtree) > 0){ # At least one taxa is listed in the matrix, but is not in the tree
    # Make message so user knows taxa have been removed
    matchobj$message.data <- "The following taxa in the data matrix were not in the tree and were excluded from analysis: "
    first <- TRUE
    for(toRemove in 1:length(missing.fromtree)){
      if(first){
        matchobj$message.data <- paste(matchobj$message.data,as.character(data[missing.fromtree[toRemove],1]),sep="")
        first <- FALSE
      } else { #not the first one, so add leading comma
        matchobj$message.data <- paste(matchobj$message.data,", ",as.character(data[missing.fromtree[toRemove],1]),sep="")
      }
    }
    matchobj$data <- data[-missing.fromtree,] # omits those data rows which have no match in the tree
    for(datacol in 2:length(matchobj$data[1,])){
      matchobj$data[,datacol] <- factor(matchobj$data[,datacol]) # have to use factor to remove any factors not present in the final dataset
    }
  }
  
  missing.taxa <- NULL
  for(tip in 1:length(phy$tip.label)){
    if(is.na(match(phy$tip.label[tip],matchobj$data[,1]))){
      if(is.null(matchobj$message.tree)){ # The first missing taxon
        missing.taxa <- as.character(phy$tip.label[tip])
        matchobj$message.tree <- "The following taxa were in the tree but did not have corresponding data in the data matrix.  They are coded as missing data for subsequent analyses: "
      } else { # not the first missing taxon, add with leading comma
        missing.taxa <- paste(missing.taxa,", ",as.character(phy$tip.label[tip]),sep="")
      }
      # missing taxa will be coded as having missing data "?"
      addtaxon <- as.character(phy$tip.label[tip])
      numcols <- length(matchobj$data[1,])
      newrow <- matrix(as.character("\x3F"),1,numcols) # absurd, but it works
      newrow[1,1] <- addtaxon
      newrowdf <- data.frame(newrow)
      colnames(newrowdf) <- colnames(matchobj$data)
      matchobj$data <- rbind(matchobj$data,newrowdf)
    }
  }
  rownames(matchobj$data) <- matchobj$data[,1] # Use first column (taxon names) as row names
  matchobj$data <- matchobj$data[matchobj$phy$tip.label,] # Sort by order in tree
  rownames(matchobj$data) <- NULL # remove row names after sorting
  if(!is.null(missing.taxa)){
    matchobj$message.tree <- paste(matchobj$message.tree,missing.taxa,sep="")
  }
  return(matchobj)
}
