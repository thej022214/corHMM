
######################################################################################################################################
######################################################################################################################################
### tipwise homology probability functions
######################################################################################################################################
######################################################################################################################################

## get descendents
dev.get_desc <- function(phy, node){
  desc <- phy$edge[phy$edge[,1] == node, 2]
  return(desc)
}

get_desc <- function(phy, node){
  to_check <- node
  n_tip <- Ntip(phy)
  desc <- c()
  while(length(to_check) > 0){
    curr_desc <- dev.get_desc(phy, to_check[1])
    desc <- c(desc, curr_desc)
    curr_desc < curr_desc[curr_desc > n_tip]
    to_check <- c(to_check, curr_desc)
    to_check <- to_check[-1]
  }
  return(desc)
}

#written by James Boyko
tipHomology <- function(corhmm_obj, type="strict", node=NULL, return.likelihoods = FALSE){
  
  phy=corhmm_obj$phy
  data=corhmm_obj$data
  # some prereqs
  p = sapply(1:max(corhmm_obj$index.mat, na.rm = TRUE), function(x) 
    na.omit(c(corhmm_obj$solution))[na.omit(c(corhmm_obj$index.mat) == x)][1])
  Q <- corhmm_obj$solution
  Q[is.na(Q)] <- 0
  diag(Q) <- -rowSums(Q)
  p_mat_by_edge <- vapply(corhmm_obj$phy$edge.length, function(x) expm(Q * x, method = "Ward77"), 
    FUN.VALUE = matrix(0, nrow(Q), ncol(Q)))
  corData <- corProcessData(data, collapse = corhmm_obj$collapse)
  model.set.final <- rate.cat.set.corHMM.JDB(
    phy=corhmm_obj$phy, 
    data=corhmm_obj$data, 
    rate.cat=corhmm_obj$rate.cat, 
    ntraits = length(corData$ObservedTraits), 
    model = "ARD", 
    rate.mat = corhmm_obj$rate.mat, 
    collapse = corhmm_obj$collapse)
  
  if(is.null(node)){
    tips <- 1:Ntip(phy)
  }else{
    tips <- get_desc(phy, node)
    tips <- tips[tips <= Ntip(phy)]
  }
  tip_liks <- sapply(tips, function(x) 
    dev.tip_homology(phy = corhmm_obj$phy,
      corData = corData, 
      p = p, 
      rate.cat = corhmm_obj$rate.cat, 
      tip = x, 
      node = node, 
      type = type, 
      rate.mat = corhmm_obj$index.mat , 
      root.p = corhmm_obj$root.p, 
      collapse = corhmm_obj$collapse, 
      p_mat_by_edge = p_mat_by_edge,
      model.set.final = model.set.final))
  tip_liks <- setNames(tip_liks, phy$tip.label[tips])
  if(return.likelihoods){
    return(tip_liks)
  }else{
    tip_liks <- exp(tip_liks - corhmm_obj$loglik)
    return(tip_liks)
  }
}


dev.tip_homology <- function(phy, corData, p, rate.cat, tip, node, type, rate.mat=NULL, root.p=NULL, collapse = TRUE, p_mat_by_edge, model.set.final){
	
	data <- corData
	
	#Ensures that weird root state probabilities that do not sum to 1 are input:
	if(!is.null(root.p)){
		if(!is.character(root.p)){
			root.p <- root.p/sum(root.p)
		}
	}
	# because we overwrite root.p below and i need it later
	root.p_input <- root.p
	if (is.null(rate.cat)){
		rate.cat <- 1
	}
	model="ARD"
	
	#Note: Does not like zero branches at the tips. Here I extend these branches by just a bit:
	phy$edge.length[phy$edge.length<=1e-5]=1e-5
	data.sort <- data.frame(data[,2], data[,2],row.names=data[,1])
	data.sort <- data.sort[phy$tip.label,]
	levels <- levels(as.factor(data.sort[,1]))
	
	#Some initial values for use later
	obj <- NULL
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	
	k <- ntraits <- length(corData$ObservedTraits)
	drop.states <- NULL
	if(is.null(rate.mat)){
		rate.mat <- model.set.final$index.matrix
		rate <- model.set.final$rate
		num.dropped.states <- NULL
	}else{
		rate <- rate.mat
		col.sums <- which(colSums(rate.mat, na.rm=TRUE) == 0)
		row.sums <- which(rowSums(rate.mat, na.rm=TRUE) == 0)
		drop.states <- col.sums[which(col.sums == row.sums)]
		if(length(drop.states > 0)){
			model.set.final$liks[,drop.states] <- 0
			num.dropped.states <- length(drop.states)
		}else{
			num.dropped.states <- NULL
		}
		rate[is.na(rate)] <- max(rate,na.rm=TRUE)+1
	}
	
	x <- data.sort[,1]
	TIPS <- 1:nb.tip
	tranQ <- Q <- model.set.final$Q
	liks <- model.set.final$liks
	
	if(length(drop.states > 0)){
		liks[,drop.states] <- 0
	}
	
	p[p==0] = exp(-21)
	Q[] <- c(p, 0)[rate]
	diag(Q) <- -rowSums(Q)
	phy <- reorder(phy, "pruningwise")
	TIPS <- 1:nb.tip
	anc <- unique(phy$edge[,1])
	
	#A temporary likelihood matrix so that the original does not get written over:
	liks.down <- liks
	focal_state <- liks[tip,]
	tip_index <- which(liks[tip,] == 1)
	#A transpose of Q for assessing probability of j to i, rather than i to j:
	tranQ <- t(Q)
	comp <- numeric(nb.tip + nb.node)
	#The first down-pass: The same algorithm as in the main function to calculate the conditional likelihood at each node:
	for (i in seq(from = 1, length.out = nb.node)) {
		#the ancestral node at row i is called focal
		focal <- anc[i]
		#Get descendant information of focal
		desRows<-which(phy$edge[,1]==focal)
		desNodes<-phy$edge[desRows,2]
		v <- 1
		for (desIndex in sequence(length(desRows))){
		  if(tip %in% desNodes){
		    if(tip %in% desNodes[desIndex]){
		      if(type == "strict"){
		        v <- v * focal_state * exp(diag(Q)[tip_index] * phy$edge.length[desRows[desIndex]])
		      }
		      if(type == "wagner"){
		        v <- v * focal_state * p_mat_by_edge[,,desRows[desIndex]] %*% liks.down[desNodes[desIndex],]
		      }
		    }else{
		      v <- v * p_mat_by_edge[,,desRows[desIndex]] %*% liks.down[desNodes[desIndex],]
		    }
		    v <- v * focal_state
		  }else{
		    v <- v * p_mat_by_edge[,,desRows[desIndex]] %*% liks.down[desNodes[desIndex],]
		  }
		}
		tip <- ifelse(tip %in% desNodes, focal, tip)
		if(!is.null(node)){
		  tip <- ifelse(tip == node, NA, tip)
		}
		comp[focal] <- sum(v)
		liks.down[focal, ] <- v/comp[focal]
	}
	root <- nb.tip + 1L
	# liks.down[root,] <- focal_state
	#Enter the root defined root probabilities if they are supplied by the user:
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
		liks.down[root,  ] <- flat.root * liks.down[root,  ]
		liks.down[root,  ] <- liks.down[root, ] / sum(liks.down[root,  ])
		root.p = flat.root
	}else{
		if(is.character(root.p)){
			# root.p==yang will fix root probabilities based on the inferred rates: q10/(q01+q10), q01/(q01+q10), etc.
			if(root.p == "yang"){
				root.p <- Null(Q)
				root.p <- c(root.p/sum(root.p))
			}else{
				# root.p==maddfitz will fix root probabilities according to FitzJohn et al 2009 Eq. 10:
				root.p = liks.down[root, ] / sum(liks.down[root, ])
			}
		}
	}
	loglik <- (sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks.down[root,])))))
  return(loglik)
}


