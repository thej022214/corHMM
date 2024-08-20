
######################################################################################################################################
######################################################################################################################################
### Standalone ancestral state reconstruction
######################################################################################################################################
######################################################################################################################################

#written by Jeremy M. Beaulieu and Jeffrey C. Oliver

ancRECON <- function(phy, data, p, method=c("joint", "marginal", "scaled"), rate.cat, ntraits=NULL, rate.mat=NULL, model="ARD", root.p=NULL, get.likelihood=FALSE, get.tip.states = FALSE, tip.fog=NULL, collapse = TRUE){
	
  # if(hasArg(corHMM_fit)){
  #   corHMM_fit$phy$node.label <- NULL
  #   phy = corHMM_fit$phy
  #   data = corHMM_fit$data
  #   p = sapply(na.omit(c(corHMM_fit$index.mat)), function(x) na.omit(c(corHMM_fit$solution))[x])
  #   rate.cat = corHMM_fit$rate.cat
  #   rate.mat = corHMM_fit$index.mat
  #   root.p = corHMM_fit$root.p
  # }else{
  #   if(!hasArg(phy) & hasArg(data) & hasArg(p) & hasArg(rate.cat)){
  #     return(cat("Warning: please provide either a corHMM results object or a phylogeny, dataset, parameter vector, and rate category."))
  #   }
  # }
	
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
	
	#data consistency stuff
	input.data <- data
	corData <- corProcessData(data, collapse = collapse)
	data <- corData$corData

	matching <- match.tree.data(phy,data)
	data <- matching$data
	phy <- matching$phy
	
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
		model.set.final <- rate.cat.set.corHMM.JDB(phy=phy,data=input.data, rate.cat=rate.cat, ntraits = ntraits, model = model, rate.mat = rate.mat, collapse = collapse)
		rate.mat <- model.set.final$index.matrix
		rate <- model.set.final$rate
		num.dropped.states <- NULL
	}else{
		model.set.final <- rate.cat.set.corHMM.JDB(phy=phy,data=input.data, rate.cat=rate.cat, ntraits = ntraits, model = model, rate.mat=rate.mat, collapse = collapse)
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
	
	#if((max(na.omit(as.vector(rate)))-1) < length(p)){
	#    return(cat("You have given a vector of transition rates greater than the number of parameters in the model."))
	#}
	#if((max(na.omit(as.vector(rate)))-1) > length(p)){
	#    return(cat("You have given a vector of transition rates less than the number of parameters in the model."))
	#}
	
	#Makes a matrix of tip states and empty cells corresponding
	#to ancestral nodes during the optimization process.
	if(is.null(tip.fog)){
		tip.fog <- numeric(length(levels))
	}
	
	if(sum(tip.fog) != 0){
		if(length(tip.fog) == 1){
			#Default option, but need to replicate these values across the observed states
			tip.fog <- rep(tip.fog, dim(model.set.final$Q)[2])
		}
		if(rate.cat > 1){
			#Error only applies to observed states, but need to replicate across the rate categories:
			tip.fog <- rep(tip.fog, rate.cat)
			for(tip.index in 1:Ntip(phy)){
				#Why is this here? What happens if someone does not know the state. We would code all states as 1. So here, we just alter if there are zeros for a tip:
				num.zeros <- length(model.set.final$liks[tip.index,which(model.set.final$liks[tip.index,]==0)])
				if(num.zeros > 0){
					model.set.final$liks[tip.index,which(model.set.final$liks[tip.index,]==1)] <- 1 - (sum(tip.fog[which(model.set.final$liks[tip.index,]!=1)])/rate.cat)
					model.set.final$liks[tip.index,which(model.set.final$liks[tip.index,]==0)] <- tip.fog[which(model.set.final$liks[tip.index,]==0)]
				}
			}
		}else{
			for(tip.index in 1:Ntip(phy)){
				#Why is this here? What happens if someone does not know the state. We would code all states as 1. So here, we just alter if there are zeros for a tip:
				num.zeros <- length(model.set.final$liks[tip.index,which(model.set.final$liks[tip.index,]==0)])
				if(num.zeros > 0){
					model.set.final$liks[tip.index,which(model.set.final$liks[tip.index,]==1)] <- 1 - sum(tip.fog[which(model.set.final$liks[tip.index,]!=1)])
					model.set.final$liks[tip.index,which(model.set.final$liks[tip.index,]==0)] <- tip.fog[which(model.set.final$liks[tip.index,]==0)]
				}
			}
		}
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
	
	# # if the q matrix has columns not estimated, remove them
	# row2rm <- apply(rate, 1, function(x) all(x == max(rate)))
	# col2rm <- apply(rate, 2, function(x) all(x == max(rate)))
	# Q.root <- Q[!row2rm | !col2rm, !row2rm | !col2rm]
	
	if(method=="joint"){
		if(!is.null(phy$node.label)){
			tip.state.vector <- rep(NA, Ntip(phy))
			#We remove the first, because the root state probability comes in through root.p:
			known.state.vector <- phy$node.label
			known.state.vector <- c(tip.state.vector, known.state.vector)
		}else{
			tip.state.vector <- rep(NA, Ntip(phy))
			known.state.vector <- rep(NA, Nnode(phy))
			known.state.vector <- c(tip.state.vector, known.state.vector)
		}
		lik.states<-numeric(nb.tip + nb.node)
		pupko.L <- matrix(NA,nrow=nb.tip + nb.node,ncol(liks))
		pupko.C <- matrix(NA,nrow=nb.tip + nb.node,ncol(liks))
		for (i  in seq(from = 1, length.out = nb.node)) {
			#The ancestral node at row i is called focal:
			focal <- anc[i]
			#Get descendant information of focal:
			desRows<-which(phy$edge[,1]==focal)
			#Get node information for each descendant:
			desNodes<-phy$edge[desRows,2]
			
			#Initiates a loop to check if any nodes are tips:
			for (desIndex in sequence(length(desRows))){
				#If a tip calculate C_y(i) for the tips and stores in liks matrix:
				if(any(desNodes[desIndex]==phy$edge[,1])==FALSE){
					v <- c(rep(1, k*rate.cat))
					Pij <- expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77"))
					#Pij <- matrix(c(0.7, 0.45, 0.3, 0.55), 2, 2)
					v <- v * liks[desNodes[desIndex],]
					L <- Pij %*% v
					#liks: rows are taxa + internal nodes, cols are # states
					if(is.na(known.state.vector[focal])){
						pupko.L[desNodes[desIndex],] <- L
						pupko.C[desNodes[desIndex],] <- which.is.max(L==max(L))
					}else{
						pupko.L[desNodes[desIndex],] <- L[known.state.vector[focal],]
						pupko.C[desNodes[desIndex],] <- known.state.vector[focal]
					}
				}
			}
			#Collects t_z, or the branch subtending focal:
			tz <- phy$edge.length[which(phy$edge[,2] == focal)]
			if(length(tz)==0){
				#The focal node is the root, calculate P_k:
				root.state=1
				for (desIndex in sequence(length(desRows))){
					#This is the basic marginal calculation:
					root.state <- root.state * pupko.L[desNodes[desIndex],]
				}
				if(is.na(known.state.vector[focal])){
					equil.root <- NULL
					for(i in 1:ncol(Q)){
						posrows <- which(Q[,i] >= 0)
						rowsum <- sum(Q[posrows,i])
						poscols <- which(Q[i,] >= 0)
						colsum <- sum(Q[i,poscols])
						equil.root <- c(equil.root,rowsum/(rowsum+colsum))
					}
					if (is.null(root.p)){
						if(is.na(known.state.vector[focal])){
							flat.root = equil.root
							k.rates <- 1/length(which(!is.na(equil.root)))
							flat.root[!is.na(flat.root)] = k.rates
							flat.root[is.na(flat.root)] = 0
							#root.p <- c(.6,.4)
							root.p <- flat.root
							pupko.L[focal, ] <- root.state
						}else{
							root.p <- rep(0, dim(Q)[2])
							root.p[known.state.vector[focal]] <- 1
							pupko.L[focal, ] <- root.state
						}
					}
					else{
						if(is.character(root.p)){
							# root.p==yang will fix root probabilities based on the inferred rates: q10/(q01+q10), q01/(q01+q10), etc.
							if(root.p == "yang"){
								root.p <- Null(Q)
								root.p <- c(root.p/sum(root.p))
								pupko.L[focal, ] <- root.state
							}else{
								# root.p==maddfitz will fix root probabilities according to FitzJohn et al 2009 Eq. 10:
								root.p <- (root.state / sum(root.state))[ ]
								pupko.L[focal,] <- root.state
							}
						}
						# root.p!==NULL will fix root probabilities based on user supplied vector:
						else{
							root.p <- root.p
							pupko.L[focal, ] <- root.state
						}
					}
				}else{
					root.p = rep(0, dim(Q)[1])
					root.p[known.state.vector[focal]] <- 1
					pupko.L[focal, ] <- root.state
				}
			}
			#All other internal nodes, except the root:
			else{
				#Calculates P_ij(t_z):
				Pij <- expm(Q * tz, method=c("Ward77"))
				#Pij <- matrix(c(0.7, 0.45, 0.3, 0.55), 2, 2)
				#Calculates L_z(i):
				v <- c(rep(1, k*rate.cat))
				if(is.na(known.state.vector[focal])){
					for (desIndex in sequence(length(desRows))){
						v <- v * pupko.L[desNodes[desIndex],]
					}
					focalRow <- which(phy$edge[,2]==focal)
					motherRow <- which(phy$edge[,1]==phy$edge[focalRow,1])
					motherNode <- phy$edge[focalRow,1]
					if(is.na(known.state.vector[motherNode])){
						for(row.index in 1:dim(Pij)[1]){
							L <- Pij[row.index,] * v
							pupko.L[focal, row.index] <- max(L)
							pupko.C[focal, row.index] <- which.is.max(L)
						}
					}else{
						L <- Pij[known.state.vector[motherNode],] * v
						pupko.L[focal,] <- L
						pupko.C[focal,] <- which.is.max(L)
					}
				}else{
					for (desIndex in sequence(length(desRows))){
						v <- v * pupko.L[desNodes[desIndex],]
					}
					focalRow <- which(phy$edge[,2] == focal)
					motherRow <- which(phy$edge[,1] == phy$edge[focalRow,1])
					motherNode <- phy$edge[focalRow,1]
					if(is.na(known.state.vector[motherNode])){
						for(row.index in 1:dim(Pij)[1]){
							L <- Pij[row.index,] * v
							pupko.L[focal, row.index] <- L[known.state.vector[focal]]
							pupko.C[focal, row.index] <- known.state.vector[focal]
						}
					}else{
						L <- Pij[known.state.vector[motherNode],] * v
						pupko.L[focal,] <- L[known.state.vector[focal]]
						pupko.C[focal,] <- known.state.vector[focal]
					}
				}
				if(sum(pupko.L[focal,])<1e-200){
					cat("Kicking in arbitrary precision package Rmpfr due to very low probabilities.\n")
					#Kicks in arbitrary precision calculations:
					pupko.L <- mpfr(pupko.L, 15)
				}
			}
		}
		root <- nb.tip + 1L
		if(get.likelihood == TRUE){
			loglik <- log(sum(exp(log(root.p)+log(pupko.L[root, ]))))
			return(as.numeric(loglik))
		}else{
			root <- nb.tip + 1L
			if(is.na(known.state.vector[root])){
				pupko.L[root, ] <- log(root.p)+log(pupko.L[root, ])
				lik.states[root] <- which(pupko.L[root,] == max(pupko.L[root,]))[1]
			}else{
				lik.states[root] <- known.state.vector[root]
			}
			N <- dim(phy$edge)[1]
			for(i in N:1){
				anc <- phy$edge[i,1]
				des <- phy$edge[i,2]
				lik.states[des] <- pupko.C[des,lik.states[anc]]
			}
			#Outputs likeliest tip states
			obj$lik.tip.states <- lik.states[TIPS]
			#Outputs likeliest node states
			obj$lik.anc.states <- lik.states[-TIPS]
			#Outputs the information gained (in bits) per node
			obj$info.anc.states <- NULL
			return(obj)
		}
	}
	if(method=="marginal"){
		#if(is.character(root.p)){
		#Fairly certain this is right. See Maddison (1995) and what to do at the root.
		#   root.p = NULL
		#}
		#A temporary likelihood matrix so that the original does not get written over:
		liks.down <- liks
		#A transpose of Q for assessing probability of j to i, rather than i to j:
		tranQ <- t(Q)
		comp <- matrix(0,nb.tip + nb.node,ncol(liks))
		#The first down-pass: The same algorithm as in the main function to calculate the conditional likelihood at each node:
		for (i in seq(from = 1, length.out = nb.node)) {
			#the ancestral node at row i is called focal
			focal <- anc[i]
			#Get descendant information of focal
			desRows<-which(phy$edge[,1]==focal)
			desNodes<-phy$edge[desRows,2]
			v <- 1
			for (desIndex in sequence(length(desRows))){
				v <- v*expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks.down[desNodes[desIndex],]
			}
			
			##Allows for fixed nodes based on user input tree.
			if(!is.null(phy$node.label)){
				if(!is.na(phy$node.label[focal - nb.tip])){
					fixer.tmp = numeric(dim(Q)[2]/rate.cat)
					fixer.tmp[phy$node.label[focal - nb.tip]] = 1
					fixer = rep(fixer.tmp, rate.cat)
					v <- v * fixer
				}
			}
			
			comp[focal] <- sum(v)
			liks.down[focal, ] <- v/comp[focal]
		}
		root <- nb.tip + 1L
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
					liks.down[root,  ] <-  liks.down[root,  ] * root.p
					liks.down[root,  ] <- liks.down[root, ] / sum(liks.down[root, ])
				}else{
					# root.p==maddfitz will fix root probabilities according to FitzJohn et al 2009 Eq. 10:
					root.p = liks.down[root, ] / sum(liks.down[root, ])
					liks.down[root,  ] <- root.p * liks.down[root,  ]
					liks.down[root,  ] <- liks.down[root, ] / sum(liks.down[root,  ])
				}
			}else{
				liks.down[root,  ] <- root.p * liks.down[root,  ]
				liks.down[root,  ] <- liks.down[root, ] / sum(liks.down[root,  ])
			}
		}
		#The up-pass
		obsRoot <- root.p
		root.p <- vector("numeric", dim(liks)[2])
		root.p[ ] <- obsRoot
		liks.up <- liks
		states<-apply(liks,1,which.max)
		N <- dim(phy$edge)[1]
		comp <- numeric(nb.tip + nb.node)
		for(i in length(anc):1){
			focal <- anc[i]
			if(!focal==root){
				#Gets mother and sister information of focal:
				focalRow <- which(phy$edge[,2]==focal)
				motherRow <- which(phy$edge[,1]==phy$edge[focalRow,1])
				motherNode <- phy$edge[focalRow,1]
				desNodes <- phy$edge[motherRow,2]
				sisterNodes <- desNodes[(which(!desNodes==focal))]
				sisterRows <- which(phy$edge[,2]%in%sisterNodes==TRUE)
				#If the mother is not the root then you are calculating the probability of being in either state.
				#But note we are assessing the reverse transition, j to i, rather than i to j, so we transpose Q to carry out this calculation:
				if(motherNode != root){
					v <- expm(tranQ * phy$edge.length[which(phy$edge[,2]==motherNode)], method=c("Ward77")) %*% liks.up[motherNode,]
					#Allows for fixed nodes based on user input tree.
					if(!is.null(phy$node.label)){
						if(!is.na(phy$node.label[motherNode - nb.tip])){
							fixer.tmp <- numeric(dim(Q)[2]/rate.cat)
							fixer.tmp[phy$node.label[motherNode - nb.tip]] <- 1
							fixer <- rep(fixer.tmp, rate.cat)
							v <- v * fixer
						}
					}
				}else{
					#If the mother is the root then just use the marginal. This can also be the prior, which I think is the equilibrium frequency.
					#But for now we are just going to use the marginal at the root -- it is unclear what Mesquite does.
					v <- root.p
				}
				#Now calculate the probability that each sister is in either state. Sister can be more than 1 when the node is a polytomy.
				#This is essentially calculating the product of the mothers probability and the sisters probability:
				for (sisterIndex in sequence(length(sisterRows))){
					v <- v * expm(Q * phy$edge.length[sisterRows[sisterIndex]], method=c("Ward77")) %*% liks.down[sisterNodes[sisterIndex],]
				}
				
				comp[focal] <- sum(v)
				liks.up[focal,] <- v/comp[focal]
			}
		}
		#The final pass
		liks.final <- liks
		comp <- numeric(nb.tip + nb.node)
		#In this final pass, root is never encountered. But its OK, because root likelihoods are set after the loop:
		for (i in seq(from = 1, length.out = nb.node-1)) {
			#the ancestral node at row i is called focal
			focal <- anc[i]
			focalRows <- which(phy$edge[,2]==focal)
			#Now you are assessing the change along the branch subtending the focal by multiplying the probability of
			#everything at and above focal by the probability of the mother and all the sisters given time t:
			v <- liks.down[focal,] * (expm(tranQ * phy$edge.length[focalRows], method=c("Ward77")) %*% liks.up[focal,])
			comp[focal] <- sum(v)
			liks.final[focal, ] <- v/comp[focal]
		}
		if(get.tip.states == TRUE){
			#Now get the states for the tips (will do, not available for general use):
			liks.final[TIPS,] <- GetTipStateBruteForce(p=p, phy=phy, data=input.data, rate.mat=rate.mat, rate.cat=rate.cat, ntraits=ntraits, model=model, root.p=root.p_input, collapse = collapse, tip.fog = tip.fog)
		}else{
			liks.final[TIPS,] <- liks.down[TIPS,]
		}
		#Just add in the marginal at the root calculated on the original downpass or if supplied by the user:
		liks.final[root,] <- liks.down[root,]
		#root.final <- liks.down[root,] * root.p
		#comproot <- sum(root.final)
		#liks.final[root,] <- root.final/comproot
		
		if(get.likelihood == TRUE){
			#NEED TO FIGURE OUT LOG COMPENSATION ISSUE --- see line 397.
			loglik <- as.numeric(log(liks[root,lik.states[root]]))
			return(loglik)
		}else{
			#Outputs likeliest tip states
			obj$lik.tip.states <- liks.final[TIPS,]
			#Outputs likeliest node states
			obj$lik.anc.states <- liks.final[-TIPS,]
			#Outputs the information gained (in bits) per node
			obj$info.anc.states <- getInfoPerNode(obj$lik.anc.states, Q)
			return(obj)
		}
	}
	
	if(method=="scaled"){
		comp<-matrix(0,nb.tip + nb.node,ncol(liks))
		root <- nb.tip + 1L
		#The same algorithm as in the main function. See comments in either corHMM.R, corDISC.R, or rayDISC.R for details:
		for (i  in seq(from = 1, length.out = nb.node)) {
			#the ancestral node at row i is called focal
			focal <- anc[i]
			#Get descendant information of focal
			desRows<-which(phy$edge[,1]==focal)
			desNodes<-phy$edge[desRows,2]
			v <- 1
			for (desIndex in sequence(length(desRows))){
				v <- v*expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks[desNodes[desIndex],]
			}
			comp[focal] <- sum(v)
			liks[focal, ] <- v/comp[focal]
		}
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
			liks[root, ] <- flat.root * liks[root, ]
			liks[root, ] <- liks[root, ] / sum(liks[root, ])
		}else{
			if(is.character(root.p)){
				# root.p==yang will fix root probabilities based on the inferred rates: q10/(q01+q10), q01/(q01+q10), etc.
				if(root.p == "yang"){
					root.p <- Null(Q)
					root.p <- c(root.p/sum(root.p))
					liks[root, ] <- root.p * liks[root,  ]
					liks[root, ] <- liks[root, ] / sum(liks[root, ])
				}else{
					# root.p==maddfitz will fix root probabilities according to FitzJohn et al 2009 Eq. 10:
					root.p = liks[root, ] / sum(liks[root, ])
					liks[root, ] <- root.p * liks[root, ]
					liks[root, ] <- liks[root, ] / sum(liks[root, ])
				}
			}
			# root.p!==NULL will fix root probabilities based on user supplied vector:
			else{
				liks[root, ] <- root.p * liks[root, ]
				liks[root, ] <- liks[root, ] / sum(liks[root, ])
			}
		}
		#Reports the probabilities for all internal nodes as well as tips:
		obj$lik.tip.states <- liks[TIPS,]
		#Outputs likeliest node states
		obj$lik.anc.states <- liks[-TIPS,]
		#Outputs the information gained (in bits) per node
		obj$info.anc.states <- NULL
		return(obj)
	}
}


######################################################################################################################################
######################################################################################################################################
### Utility functions
######################################################################################################################################
######################################################################################################################################

getInfoPerNode <- function(lik.anc.states, Q){
  Eq <- c(Null(Q))/sum(Null(Q))
  #nStates <- dim(lik.anc.states)[2]
  #H_uncond <- sum(1/nStates * -log2(rep(1/nStates, nStates)))
  H_uncond <- sum(Eq * -log2(Eq))
  H_cond <- rowSums(lik.anc.states * -log2(lik.anc.states))
  Info <- H_uncond - H_cond
  return(Info)
}


GetTipStateBruteForce <- function(p, phy, data, rate.mat, rate.cat, ntraits, model, root.p, collapse = TRUE, tip.fog = NULL){
	
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	
	data.for.likelihood.function <- corHMM:::rate.cat.set.corHMM.JDB(phy=phy, data=data, rate.cat=rate.cat, ntraits = ntraits, model = model, collapse = collapse)
	
	if(!is.null(rate.mat)){
		rate <- rate.mat
		data.for.likelihood.function$np <- max(rate, na.rm=TRUE)
		rate[is.na(rate)]=max(rate, na.rm=TRUE)+1
		data.for.likelihood.function$rate <- rate
		data.for.likelihood.function$index.matrix <- rate.mat
		## for precursor type models ##
		col.sums <- which(colSums(rate.mat, na.rm=TRUE) == 0)
		row.sums <- which(rowSums(rate.mat, na.rm=TRUE) == 0)
		drop.states <- col.sums[which(col.sums == row.sums)]
		if(length(drop.states > 0)){
			data.for.likelihood.function$liks[,drop.states] <- 0
			num.dropped.states <- length(drop.states)
		}else{
			num.dropped.states <- NULL
		}
		###############################
	}else{
		num.dropped.states <- NULL
	}
	
	if(!is.null(tip.fog)){
		if(length(tip.fog) == 1){
			#Default option, but need to replicate these values across the observed states
			tip.fog <- rep(tip.fog, dim(data.for.likelihood.function$Q)[2])
		}
		if(rate.cat > 1){
			#Error only applies to observed states, but need to replicate across the rate categories:
			tip.fog <- rep(tip.fog, rate.cat)
			for(tip.index in 1:Ntip(phy)){
				#Why is this here? What happens if someone does not know the state. We would code all states as 1. So here, we just alter if there are zeros for a tip:
				num.zeros <- length(data.for.likelihood.function$liks[tip.index,which(data.for.likelihood.function$liks[tip.index,]==0)])
				if(num.zeros > 0){
					data.for.likelihood.function$liks[tip.index,which(data.for.likelihood.function$liks[tip.index,]==1)] <- 1 - (sum(tip.fog[which(data.for.likelihood.function$liks[tip.index,]!=1)])/rate.cat)
					data.for.likelihood.function$liks[tip.index,which(data.for.likelihood.function$liks[tip.index,]==0)] <- tip.fog[which(data.for.likelihood.function$liks[tip.index,]==0)]
				}
			}
		}else{
			for(tip.index in 1:Ntip(phy)){
				#Why is this here? What happens if someone does not know the state. We would code all states as 1. So here, we just alter if there are zeros for a tip:
				num.zeros <- length(data.for.likelihood.function$liks[tip.index,which(data.for.likelihood.function$liks[tip.index,]==0)])
				if(num.zeros > 0){
					data.for.likelihood.function$liks[tip.index,which(data.for.likelihood.function$liks[tip.index,]==1)] <- 1 - sum(tip.fog[which(data.for.likelihood.function$liks[tip.index,]!=1)])
					data.for.likelihood.function$liks[tip.index,which(data.for.likelihood.function$liks[tip.index,]==0)] <- tip.fog[which(data.for.likelihood.function$liks[tip.index,]==0)]
				}
			}
		}
	}

	nodes <- unique(phy$edge[,1])
	marginal.probs <- matrix(0, nb.tip, dim(data.for.likelihood.function$Q)[2])
	total.state.space <- 1:dim(data.for.likelihood.function$Q)[2]
	unique.state.space <- dim(data.for.likelihood.function$Q)[2] / rate.cat

	for(taxon.index in 1:Ntip(phy)){
		if(!is.null(tip.fog)){
			marginal.probs.tmp <- numeric(dim(data.for.likelihood.function$Q)[2])
			states.keep <- data.for.likelihood.function$liks[taxon.index,]
			nstate.index <- total.state.space[total.state.space%%unique.state.space == 1]
			count <- 1
			nstates <- c()
			for(index in nstate.index){
				tmp.location <- c(index:(unique.state.space * count))
				tmp.state <- which.max(data.for.likelihood.function$liks[taxon.index,tmp.location])
				nstates <- c(nstates, tmp.location[tmp.state])
				count <- count + 1
			}
			count <- 1
			for(state.index in setdiff(nstates, drop.states)){
				data.for.likelihood.function$liks[taxon.index,] <- 0
				data.for.likelihood.function$liks[taxon.index,c(index:(unique.state.space * count))] <- states.keep[1:unique.state.space]
				marginal.probs.tmp[state.index] <- -corHMM:::dev.corhmm(p=log(p), phy=phy, liks=data.for.likelihood.function$liks, Q=data.for.likelihood.function$Q, rate=data.for.likelihood.function$rate, root.p=root.p, rate.cat = rate.cat, order.test = FALSE, lewis.asc.bias = FALSE)
				count <- count + 1
			}
			data.for.likelihood.function$liks[taxon.index,] <- states.keep
			best.probs = max(marginal.probs.tmp[nstates])
			marginal.probs.rescaled = marginal.probs.tmp[nstates] - best.probs
			marginal.probs[taxon.index,nstates] = exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
			
		}else{
			marginal.probs.tmp <- numeric(dim(data.for.likelihood.function$Q)[2])
			nstates <- which(!data.for.likelihood.function$liks[taxon.index,] == 0)
			states.keep <- data.for.likelihood.function$liks[taxon.index,]
			for(state.index in setdiff(nstates, drop.states)){
				data.for.likelihood.function$liks[taxon.index,] <- 0
				data.for.likelihood.function$liks[taxon.index,state.index] = 1
				marginal.probs.tmp[state.index] <- -dev.corhmm(p=log(p), phy=phy, liks=data.for.likelihood.function$liks, Q=data.for.likelihood.function$Q, rate=data.for.likelihood.function$rate, root.p=root.p, rate.cat = rate.cat, order.test = FALSE, lewis.asc.bias = FALSE)
			}
			data.for.likelihood.function$liks[taxon.index,] <- states.keep
			best.probs = max(marginal.probs.tmp[nstates])
			marginal.probs.rescaled = marginal.probs.tmp[nstates] - best.probs
			marginal.probs[taxon.index,nstates] = exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
		}
	}
	tip.states <- marginal.probs[1:nb.tip,]
	rownames(tip.states) <- phy$tip.label
	return(tip.states)
}


# dev.ancRECON.marginal <- function(p, phy, liks, Q, rate, root.p, rate.cat, get.tip.states){
#
#   obj <- NULL
#   nb.tip <- length(phy$tip.label)
#   nb.node <- phy$Nnode
#   phy <- reorder(phy, "pruningwise")
#   TIPS <- 1:nb.tip
#   anc <- unique(phy$edge[,1])
#   p[p==0] = exp(-21)
#   Q[] <- c(p, 0)[rate]
#   diag(Q) <- -rowSums(Q)
#   # if the q matrix has columns not estimated, remove them
#   row2rm <- apply(rate, 1, function(x) all(x == max(rate)))
#   col2rm <- apply(rate, 2, function(x) all(x == max(rate)))
#   Q.root <- Q[!row2rm | !col2rm, !row2rm | !col2rm]
#
#   #A temporary likelihood matrix so that the original does not get written over:
#   liks.down <- liks
#   #A transpose of Q for assessing probability of j to i, rather than i to j:
#   tranQ <- t(Q)
#   comp <- matrix(0,nb.tip + nb.node,ncol(liks))
#   #The first down-pass: The same algorithm as in the main function to calculate the conditional likelihood at each node:
#   for (i in seq(from = 1, length.out = nb.node)) {
#     #the ancestral node at row i is called focal
#     focal <- anc[i]
#     #Get descendant information of focal
#     desRows<-which(phy$edge[,1]==focal)
#     desNodes<-phy$edge[desRows,2]
#     v <- 1
#     for (desIndex in sequence(length(desRows))){
#       v <- v*expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks.down[desNodes[desIndex],]
#     }
#
#     ##Allows for fixed nodes based on user input tree.
#     if(!is.null(phy$node.label)){
#       if(!is.na(phy$node.label[focal - nb.tip])){
#         fixer.tmp = numeric(dim(Q)[2]/rate.cat)
#         fixer.tmp[phy$node.label[focal - nb.tip]] = 1
#         fixer = rep(fixer.tmp, rate.cat)
#         v <- v * fixer
#       }
#     }
#
#     comp[focal] <- sum(v)
#     liks.down[focal, ] <- v/comp[focal]
#   }
#   root <- nb.tip + 1L
#   #Enter the root defined root probabilities if they are supplied by the user:
#   equil.root <- NULL
#   for(i in 1:ncol(Q.root)){
#     posrows <- which(Q.root[,i] >= 0)
#     rowsum <- sum(Q.root[posrows,i])
#     poscols <- which(Q.root[i,] >= 0)
#     colsum <- sum(Q.root[i,poscols])
#     equil.root <- c(equil.root,rowsum/(rowsum+colsum))
#   }
#   if (is.null(root.p)){
#     flat.root = equil.root
#     k.rates <- 1/length(which(!is.na(equil.root)))
#     flat.root[!is.na(flat.root)] = k.rates
#     flat.root[is.na(flat.root)] = 0
#     liks.down[root, !col2rm] <- flat.root * liks.down[root, !col2rm]
#     liks.down[root, !col2rm] <- liks.down[root,!col2rm] / sum(liks.down[root, !col2rm])
#     root.p = flat.root
#   }else{
#     if(is.character(root.p)){
#       # root.p==yang will fix root probabilities based on the inferred rates: q10/(q01+q10), q01/(q01+q10), etc.
#       if(root.p == "yang"){
#         root.p <- Null(Q.root)
#         root.p <- c(root.p/sum(root.p))
#         liks.down[root, !col2rm] <-  liks.down[root, !col2rm] * root.p
#         liks.down[root, !col2rm] <- liks.down[root,!col2rm] / sum(liks.down[root,!col2rm])
#       }else{
#         # root.p==maddfitz will fix root probabilities according to FitzJohn et al 2009 Eq. 10:
#         root.p = liks.down[root,!col2rm] / sum(liks.down[root,!col2rm])
#         liks.down[root, !col2rm] <- root.p * liks.down[root, !col2rm]
#         liks.down[root, !col2rm] <- liks.down[root,!col2rm] / sum(liks.down[root, !col2rm])
#       }
#     }else{
#       liks.down[root, !col2rm] <- root.p * liks.down[root, !col2rm]
#       liks.down[root, !col2rm] <- liks.down[root,!col2rm] / sum(liks.down[root, !col2rm])
#     }
#   }
#   #The up-pass
#   obsRoot <- root.p
#   root.p <- vector("numeric", dim(liks)[2])
#   root.p[!col2rm] <- obsRoot
#   liks.up <- liks
#   states<-apply(liks,1,which.max)
#   N <- dim(phy$edge)[1]
#   comp <- numeric(nb.tip + nb.node)
#   for(i in length(anc):1){
#     focal <- anc[i]
#     if(!focal==root){
#       #Gets mother and sister information of focal:
#       focalRow <- which(phy$edge[,2]==focal)
#       motherRow <- which(phy$edge[,1]==phy$edge[focalRow,1])
#       motherNode <- phy$edge[focalRow,1]
#       desNodes <- phy$edge[motherRow,2]
#       sisterNodes <- desNodes[(which(!desNodes==focal))]
#       sisterRows <- which(phy$edge[,2]%in%sisterNodes==TRUE)
#       #If the mother is not the root then you are calculating the probability of being in either state.
#       #But note we are assessing the reverse transition, j to i, rather than i to j, so we transpose Q to carry out this calculation:
#       if(motherNode != root){
#         v <- expm(tranQ * phy$edge.length[which(phy$edge[,2]==motherNode)], method=c("Ward77")) %*% liks.up[motherNode,]
#         #Allows for fixed nodes based on user input tree.
#         if(!is.null(phy$node.label)){
#           if(!is.na(phy$node.label[motherNode - nb.tip])){
#             fixer.tmp <- numeric(dim(Q)[2]/rate.cat)
#             fixer.tmp[phy$node.label[motherNode - nb.tip]] <- 1
#             fixer <- rep(fixer.tmp, rate.cat)
#             v <- v * fixer
#           }
#         }
#       }else{
#         #If the mother is the root then just use the marginal. This can also be the prior, which I think is the equilibrium frequency.
#         #But for now we are just going to use the marginal at the root -- it is unclear what Mesquite does.
#         v <- root.p
#       }
#       #Now calculate the probability that each sister is in either state. Sister can be more than 1 when the node is a polytomy.
#       #This is essentially calculating the product of the mothers probability and the sisters probability:
#       for (sisterIndex in sequence(length(sisterRows))){
#         v <- v * expm(Q * phy$edge.length[sisterRows[sisterIndex]], method=c("Ward77")) %*% liks.down[sisterNodes[sisterIndex],]
#       }
#
#       comp[focal] <- sum(v)
#       liks.up[focal,] <- v/comp[focal]
#     }
#   }
#   #The final pass
#   liks.final <- liks
#   comp <- numeric(nb.tip + nb.node)
#   #In this final pass, root is never encountered. But its OK, because root likelihoods are set after the loop:
#   for (i in seq(from = 1, length.out = nb.node-1)) {
#     #the ancestral node at row i is called focal
#     focal <- anc[i]
#     focalRows <- which(phy$edge[,2]==focal)
#     #Now you are assessing the change along the branch subtending the focal by multiplying the probability of
#     #everything at and above focal by the probability of the mother and all the sisters given time t:
#     v <- liks.down[focal,] * (expm(tranQ * phy$edge.length[focalRows], method=c("Ward77")) %*% liks.up[focal,])
#     comp[focal] <- sum(v)
#     liks.final[focal, ] <- v/comp[focal]
#   }
#   if(get.tip.states == TRUE){
#     #Now get the states for the tips (will do, not available for general use):
#     liks.final[TIPS,] <- GetTipStateBruteForce(p=p, phy=phy, data=input.data, rate.mat=rate.mat, rate.cat=rate.cat, ntraits=ntraits, model=model, root.p=root.p_input)
#   }else{
#     liks.final[TIPS,] <- liks.down[TIPS,]
#   }
#   #Just add in the marginal at the root calculated on the original downpass or if supplied by the user:
#   liks.final[root,] <- liks.down[root,]
#   #root.final <- liks.down[root,] * root.p
#   #comproot <- sum(root.final)
#   #liks.final[root,] <- root.final/comproot
#
#   #Outputs likeliest tip states
#   obj$lik.tip.states <- liks.final[TIPS,]
#   #Outputs likeliest node states
#   obj$lik.anc.states <- liks.final[-TIPS,]
#   return(obj)
# }

######################################################################################################################################
######################################################################################################################################
### time_slice ancestral state reconstruction
######################################################################################################################################
######################################################################################################################################

#Marginal reconstruction function for our custom corHMM
GetEdgeMarginal <- function(p, phy, data, rate.mat, rate.cat, ntraits, model, root.p, collapse = TRUE, tip.fog = NULL){
  
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  data.for.likelihood.function <- rate.cat.set.corHMM.JDB(phy=phy, data=data, rate.cat=rate.cat, ntraits = ntraits, model = model, collapse = collapse)
  
  if(!is.null(rate.mat)){
	rate <- rate.mat
	data.for.likelihood.function$np <- max(rate, na.rm=TRUE)
	rate[is.na(rate)]=max(rate, na.rm=TRUE)+1
	data.for.likelihood.function$rate <- rate
	data.for.likelihood.function$index.matrix <- rate.mat
	## for precursor type models ##
	col.sums <- which(colSums(rate.mat, na.rm=TRUE) == 0)
	row.sums <- which(rowSums(rate.mat, na.rm=TRUE) == 0)
	drop.states <- col.sums[which(col.sums == row.sums)]
	if(length(drop.states > 0)){
		data.for.likelihood.function$liks[,drop.states] <- 0
		num.dropped.states <- length(drop.states)
	}else{
		drop.states <- NULL
		num.dropped.states <- NULL
	}
	###############################
  }else{
	drop.states <- NULL
	num.dropped.states <- NULL
  }
  
  if(!is.null(tip.fog) | tip.fog == 0){
	  if(length(tip.fog) == 1){
		  #Default option, but need to replicate these values across the observed states
		  tip.fog <- rep(tip.fog, dim(data.for.likelihood.function$Q)[2])
	  }
	  if(rate.cat > 1){
		  #Error only applies to observed states, but need to replicate across the rate categories:
		  tip.fog <- rep(tip.fog, rate.cat)
		  for(tip.index in 1:Ntip(phy)){
			  #Why is this here? What happens if someone does not know the state. We would code all states as 1. So here, we just alter if there are zeros for a tip:
			  num.zeros <- length(data.for.likelihood.function$liks[tip.index,which(data.for.likelihood.function$liks[tip.index,]==0)])
			  if(num.zeros > 0){
				  data.for.likelihood.function$liks[tip.index,which(data.for.likelihood.function$liks[tip.index,]==1)] <- 1 - (sum(tip.fog[which(data.for.likelihood.function$liks[tip.index,]!=1)])/rate.cat)
				  data.for.likelihood.function[tip.index,which(data.for.likelihood.function$liks[tip.index,]==0)] <- tip.fog[which(data.for.likelihood.function$liks[tip.index,]==0)]
			  }
		  }
	  }else{
		  for(tip.index in 1:Ntip(phy)){
			  #Why is this here? What happens if someone does not know the state. We would code all states as 1. So here, we just alter if there are zeros for a tip:
			  num.zeros <- length(data.for.likelihood.function$liks[tip.index,which(data.for.likelihood.function$liks[tip.index,]==0)])
			  if(num.zeros > 0){
				  data.for.likelihood.function$liks[tip.index,which(data.for.likelihood.function$liks[tip.index,]==1)] <- 1 - sum(tip.fog[which(data.for.likelihood.function$liks[tip.index,]!=1)])
				  data.for.likelihood.function$liks[tip.index,which(data.for.likelihood.function$liks[tip.index,]==0)] <- tip.fog[which(data.for.likelihood.function$liks[tip.index,]==0)]
			  }
		  }
	  }
  }

  phy <- reorder(phy, "pruningwise")
  nodes <- unique(phy$edge[,1])
  taxon.index <- grep("FakeyMcFakerson", x=phy$tip.label)
  
  marginal.probs.tmp <- numeric(dim(data.for.likelihood.function$Q)[2])
  nstates = which(!data.for.likelihood.function$liks[taxon.index,] == 0)
  states.keep = data.for.likelihood.function$liks[taxon.index,]
  for(state.index in setdiff(1:dim(data.for.likelihood.function$Q)[2], drop.states)){
	data.for.likelihood.function$liks[taxon.index,] = 0
	data.for.likelihood.function$liks[taxon.index,state.index] = 1
	# print(data.for.likelihood.function$liks)
	marginal.probs.tmp[state.index] <- -dev.corhmm(p=log(p), phy=phy, liks=data.for.likelihood.function$liks, Q=data.for.likelihood.function$Q, rate=data.for.likelihood.function$rate, root.p=root.p, rate.cat = rate.cat, order.test = FALSE, lewis.asc.bias = FALSE)
  }
  data.for.likelihood.function$liks[taxon.index,] <- states.keep
  best.probs <- max(marginal.probs.tmp[nstates])
  marginal.probs.rescaled <- marginal.probs.tmp[nstates] - best.probs
  marginal.probs.final <- exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
  
  return(marginal.probs.final)
}

# creates a table for locating where to add nodes. first column is the node index (tip or internal). second column is the time period from the node to add
getNodePlacementsForSlice <- function(phy, time_slice){
  to_add <- data.frame(node = 1:length(phy$tip.label), position = time_slice)
  nodes_younger_than_slice <- branching.times(phy)[branching.times(phy) < time_slice]
  node_to_add <- data.frame(node = as.numeric(names(nodes_younger_than_slice)), position = time_slice - nodes_younger_than_slice)
  to_add <- rbind(to_add, node_to_add)
  edge_lengths <- phy$edge.length[match(to_add[,1], phy$edge[,2])]
  to_add <- to_add[edge_lengths - to_add[,2] > 0,] # exludes nodes too far away
  return(to_add)
}

# the main function for doing a marginal time slice reconstruction
ancRECON_slice <- function(corhmm.obj, time_slice, collapse=TRUE, ncores = 1){
  if(max(time_slice) >= max(branching.times(corhmm.obj$phy))){
	stop("time_slice must be less than the maximum branching time of the tree")
  }
  # get the fake node locations for the reconstruction based on the time slice
  to_recon <- sapply(time_slice, function(x) getNodePlacementsForSlice(corhmm.obj$phy, x))
  to_recon <- (apply(to_recon, 2, function(x) do.call(cbind, x)))
  for(i in seq_len(length(time_slice))){
	to_recon[[i]] <- cbind(time_slice = time_slice[i], to_recon[[i]])
  }
  to_recon <- do.call(rbind, to_recon)
  # create a dummy tip and dummy data
  tip.name <- "FakeyMcFakerson"
  tip <- list(edge=matrix(c(2,1),1,2), tip.label=tip.name, edge.length=0, Nnode=1)
  class(tip) <- "phylo"
  new.data <- corhmm.obj$data.legend
  new.data[,1] <- as.character(new.data[,1])
  new.data <- rbind(new.data, "?")
  new.data[nrow(new.data), 1] <- "FakeyMcFakerson"
  ntraits <- max(as.numeric(corhmm.obj$data.legend$d))
  # get all the new phys
  print(paste0("Reconstructing ", nrow(to_recon), " nodes for ", length(time_slice), " time slices..."))
  corhmm.obj$phy$node.label <- NULL
  all_trees <- apply(to_recon, 1, function(x) bind.tree(corhmm.obj$phy, tip, x[2], x[3]))
  #Now lets calculate some marginals!
  p <- sapply(1:max(corhmm.obj$index.mat, na.rm = TRUE), function(x) na.omit(c(corhmm.obj$solution))[na.omit(c(corhmm.obj$index.mat) == x)][1])
  all_recon <- mclapply(all_trees, function(x) GetEdgeMarginal(p=p, phy=x, data=new.data, rate.mat=corhmm.obj$index.mat, rate.cat=corhmm.obj$rate.cat, ntraits=ntraits, root.p=corhmm.obj$root.p, model = "ARD", collapse = collapse), mc.cores = ncores)
  all_recon <-  do.call(rbind, all_recon)
  colnames(all_recon) <- colnames(corhmm.obj$solution)
  out <- cbind(to_recon, all_recon)
  return(out)
}

plot_slice_recon <- function(phy, slice_df, col=NULL){
  plot(phy, show.tip.label = FALSE)
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  node <- (lastPP$Ntip + 1):length(lastPP$xx)
  XX <- lastPP$xx[node]
  YY <- lastPP$yy[node]
  for(i in 1:nrow(slice_df)){
	focal <- slice_df[i,]
	xx <- lastPP$xx[focal[2]] - focal[3]
	yy <- lastPP$yy[focal[2]]
	pp <- focal[-c(1:3)]
	floating.pie.asp(xx, yy, pp, radius = 0.5, col = col)
  }
}

# # testing
# library(corHMM)
# library(viridis)
# data(primates)
# phy <- multi2di(primates[[1]])
# data <- primates[[2]]
# corhmm.obj <- corHMM(phy = phy, data = data, rate.cat = 1)
# test <- ancRECON_slice(corhmm.obj, time_slice = c(1, 10, 20, 30, 40, 49), collapse = TRUE, ncores = 4)
# corHMM:::plot_slice_recon(phy, test, col = viridis(3))
