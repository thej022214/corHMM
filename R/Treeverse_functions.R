# Includes all of the functions necessary for Treeverse
# all functions are new with the exception of print.corhmm which is an updated version compatible with Treeverse

mass.root <- function(phy_list, phy, output= c("newick", "phylo")) {
  if (!ape::is.rooted(phy)) stop("The species tree must be rooted prior to rooting nni trees.")

  Ntip <- length(phy$tip.label)
  root <- Ntip + 1
  root.row <- which(phy$edge[,1] == root)
  root.desc <- phy$edge[root.row,2]

  l <- if (root.desc[1] < root) 1 else length(ape::extract.clade(phy, root.desc[1])$tip.label)
  r <- if (root.desc[2] < root) 1 else length(ape::extract.clade(phy, root.desc[2])$tip.label)
	
  out.node <- if (r > l) root.desc[1] else root.desc[2]

  get_tips <- function(edge, node) {
	descendants <- list()
	current_nodes <- node
	while(length(current_nodes) > 0) {
	  new_nodes <- c()
	  for(n in current_nodes) {
		children <- edge[edge[,1] == n, 2]
		new_nodes <- c(new_nodes, children)
	  }
	  descendants <- c(descendants, list(new_nodes))
	  current_nodes <- new_nodes
	}
	all_descendants <- unlist(descendants)
	tips <- all_descendants[!all_descendants %in% edge[,1]]
	return(sort(tips))
  }

  if (out.node < root) {
	tips <- out.node
  } else {
	tips <- get_tips(phy$edge, out.node)
  }
  
  tip.labels <- phy$tip.label[tips]
  rooted_trees <- list()
  
  edge_finder <- function(tree, outgroup_tips) {
  tip_indices <- which(tree$tip.label %in% outgroup_tips)
  if (length(tip_indices) == 0) return(integer(0))
  
  is_outgroup <- logical(max(tree$edge))
  is_outgroup[tip_indices] <- TRUE
  
  repeat {
    new_nodes <- tree$edge[is_outgroup[tree$edge[,2]], 1]
    new_nodes <- new_nodes[!is_outgroup[new_nodes]]
    if (length(new_nodes) == 0) break
    is_outgroup[new_nodes] <- TRUE
  }
  
  out_edges <- which(is_outgroup[tree$edge[,1]] != is_outgroup[tree$edge[,2]])
  
  if (length(out_edges) > 0) {
    return(out_edges[1])  
  } else {
    return(integer(0))
  }
}
  for (i in 1:length(phy_list)) {
    if (inherits(phy_list, "multiPhylo")) {
      t.tree <- phy_list[[i]]
    } else {
      t.tree <- ape::read.tree(text = phy_list[i])
    }

    tp <- match(tip.labels, t.tree$tip.label)
    
    if (all(is.na(tp))) {
	  message(paste("No outgroup tip labels found in tree", i, "- skipping this tree."))
	  next
	} else if (any(is.na(tp))) {
	  tipward.node <- ape::getMRCA(t.tree, tp[!is.na(tp)])
      root.edge <- which(t.tree$edge[,2] == tipward.node)
    } else if (length(tp) == 1) {
      root.edge <- which(t.tree$edge[,2] == tp)
    } else {
      tipward.node <- ape::getMRCA(t.tree, tips)
      root.edge <- which(t.tree$edge[,2] == tipward.node)
    }
	
	if (length(root.edge) == 0) {
	  root.edge<-edge_finder(t.tree, tip.labels)
	  t.tree <- castor::root_in_edge(t.tree, root_edge = root.edge, location = 0.5)
	} else {
	  t.tree <- castor::root_in_edge(t.tree, root_edge = root.edge, location = 0.5)
	}

    if (output == "phylo") {
      rooted_trees[[length(rooted_trees) + 1]] <- t.tree
    } else if (output == "newick") {
      rooted_trees[i] <- ape::write.tree(t.tree)
    }
  }
    if (output == "phylo") {
      class(rooted_trees) <- "multiPhylo"
    } else if (output == "newick") {
      class(rooted_trees) <- "character"
    }

  return(rooted_trees)
}

mass.date <- function(phy_list, phy, method=c("chronos","treePL"), output=c("newick", "phylo"), ...) {
if (!ape::is.rooted(phy)) {
              stop("******The reference tree (phy) must be rooted prior using this function.******")
          }
if (!ape::is.ultrametric(phy)) {
              stop("******The reference tree (phy) must be ultrametric prior using this function.******")
          }

output <- match.arg(output)
method <- match.arg(method)

if (inherits(phy_list, "character")){
    trees <- lapply(phy_list, function(x) ape::read.tree(text = x))
  } else if (inherits(phy_list, "multiPhylo")) {
    trees <- phy_list
  }
root.age <- max(ape::branching.times(phy))
dated_trees <- list()

for (i in seq_along(trees)) {
  tree <- trees[[i]]
  message("Dating tree ", i, " of ", length(phy_list), " ... ")
  # ======================
  # ====== CHRONOS ======
  # ======================
  if (method == "chronos") {
    temp <- ape::chronos(tree, ...)
    temp$edge.length <- temp$edge.length * root.age
    if (is.ultrametric(temp) == FALSE) {
      temp <- phangorn::nnls.tree(cophenetic(temp), temp, method = "ultrametric", rooted = is.rooted(temp), trace = 0)
    }
    message("Tree ", i, " successfully dated. Writing output ...")
      if (output == "newick") {
        dated_trees[[i]] <- ape::write.tree(temp)
      } else if (output == "phylo") {
        dated_trees[[i]] <- temp
      }
    
  # ======================
  # ======= treePL =======
  # ======================
  } else if (method == "treePL") {
    scale <- "treePL"
    temp <- geiger::congruify.phylo(reference = phy,target = tree,tol = 1e-8,scale = scale, ...)
    if (!is.null(temp$phy)) {
      temp$phy <- phangorn::nnls.tree(cophenetic(temp$phy), temp$phy, method = "ultrametric", rooted = is.rooted(temp$phy), trace = 0)
      message("Tree ", i, " successfully dated. Writing output ...")
      if (output == "newick") {
        dated_trees[[i]] <- ape::write.tree(temp$phy)
      } else if (output == "phylo") {
        dated_trees[[i]] <- temp$phy
      }

      } else {
        message("Tree ", i, " could not be dated. Returning null result.")
         if (output == "newick") {
        dated_trees[[i]] <- ape::write.tree(tree)
      } else if (output == "phylo") {
        dated_trees[[i]] <- tree
      }
    }
  }
  
    
  }

  if (output == "phylo") {
      class(dated_trees) <- "multiPhylo"
    } else if (output == "newick") {
      class(dated_trees) <- "character"
    }

  return(dated_trees)
}

node.match<-function(ref, tar){
  anc <- unique(ref$edge[,1])
  anc <- anc[order(anc)]
  coals <- unique(tar$edge[,1])

  node.frame <- data.frame(reference = anc, coalnode = numeric(length(anc)))
  
  for (i in seq_along(anc)){
  node<-anc[i]
  ref.tips <- ape::extract.clade(ref, node)$tip.label
  ref.tips <- ref.tips[ref.tips %in% tar$tip.label]
  if (length(ref.tips) < 2){
    node.frame$coalnode[i] <- 0 
    next
  }
  tar.node <- ape::getMRCA(tar, ref.tips)
  node.frame$coalnode[i] <- tar.node
  }

return(node.frame)
# Example output 
# node.match(trees[[1]], trees[[2]])
#     reference coalnode
# 1         144      132
# 2         145      133
# 3         146      133
# 4         147      155
# 5         148      155 
}

print.corhmm<-function(x,...){
    if (inherits(x$phy, "multiPhylo") || is.list(x$phy)){
    ntips=Ntip(x$phy[[1]])
    output<-data.frame(x$loglik,x$AIC,x$AICc,x$rate.cat,ntips,x$ntrees, row.names="")
    names(output)<-c("lnL","AIC","AICc","Rate.cat","ntax","ntrees")
    cat("\nFit\n")
    print(output)
    cat("\n")
    }else {
    ntips=Ntip(x$phy)
    output<-data.frame(x$loglik,x$AIC,x$AICc,x$rate.cat,ntips, row.names="")
    names(output)<-c("lnL","AIC","AICc","Rate.cat","ntax")
    cat("\nFit\n")
    print(output)
    cat("\n")
    }

    UserStates <- gsub("_", "|", corProcessData(x$data)$PossibleTraits)
    ColNames <- paste0(colnames(x$data)[-1], collapse = "|")

    cat("Legend\n")
    print(ColNames)
    print(UserStates)
    cat("\n")
    
    param.est<- x$solution
    cat("Rates\n")
    print(param.est)
    cat("\n")

	if(!is.null(x$tip.fog.probs)){
		cat("Tip fog\n")
		print(x$tip.fog.probs)
		cat("\n")
	}
	
    if(any(x$eigval<0)){
        index.matrix <- x$index.mat
        #If any eigenvalue is less than 0 then the solution is not the maximum likelihood solution
        if (any(x$eigval<0)) {
            cat("The objective function may be at a saddle point", "\n")
        }
    }
    else{
        cat("Arrived at a reliable solution","\n")
    }

}

corHMM_treeverse_dev <- function(p, phy_list, data, weights_vector, rate.cat, rate.mat, model, root.p, tip.fog, collapse, node.states="none"){
	p <- exp(p) # Taking the exponential of p (rate transitions vector)
	ntrees <- length(phy_list) # length of nni trees
	tmp_logliks <- numeric(ntrees) # Initialize a vector to store log likelihoods -- initializing first actually saves time.
	for(phy_index in seq_len(ntrees)){
		tmp_logliks[phy_index] <- corHMM(phy_list[[phy_index]], data, rate.cat=rate.cat, rate.mat=rate.mat, model=model, p=p, root.p=root.p, tip.fog=tip.fog, collapse=collapse)$loglik
	}
	# Convert weights to log space
	log_weights <- log(weights_vector)
	# Adding is the same as multiplying in log space:
	log_probs <- log_weights + tmp_logliks
	# Use a log-sum-exp trick Jon Chang showed me:
	max_log_prob <- max(log_probs)
	log_weighted_avg_prob <- max_log_prob + log(sum(exp(log_probs - max_log_prob)))
	# Output the final weighted log likelihood
	return(-log_weighted_avg_prob)
}

corHMM_treeverse <- function(phy_list, data, weights_vector=NULL, acsr_weights=NULL, rate.cat=1, rate.mat=NULL, model = "ER", node.states = "none", fixed.nodes=FALSE, root.p="yang", tip.fog=NULL, ip=NULL, fog.ip = 0.01, get.tip.states = FALSE, collapse=TRUE, lower.bound = 1e-9, upper.bound = 100, opts=NULL){
	
	# Checks to make sure node.states is not NULL.  If it is, just returns a diagnostic message asking for value.
	# does node.states need to exist? currently ancTreeverse only
  # takes marginal
  # if we implement C* then this may be useful
  if(is.null(node.states)){
		obj <- NULL
		obj$loglik <- NULL
		obj$diagnostic <- paste("No model for ancestral states selected.  Please pass one of the following to corHMM command for parameter \'node.states\': joint, marginal, scaled, or none.")
		return(obj)
	} else { # even if node.states is not NULL, need to make sure its one of the three valid options
		valid.models <- c("joint", "marginal", "scaled", "none")
		if(!any(valid.models == node.states)){
			obj <- NULL
			obj$loglik <- NULL
			obj$diagnostic <- paste("\'",node.states, "\' is not valid for ancestral state reconstruction method.  Please pass one of the following to corHMM command for parameter \'node.states\': joint, marginal, scaled, or none.",sep="")
			return(obj)
		}
		if(length(node.states) > 1){ # User did not enter a value, so just pick marginal.
			node.states <- "marginal"
			cat("No model selected for \'node.states\'. Will perform marginal ancestral state estimation.\n")
		}
	}

	if(!inherits(phy_list, "multiPhylo")){
		obj <- NULL
		obj$loglik <- NULL
		obj$diagnostic <- paste("Only a single tree provided. Use corHMM() instead.")
		return(obj)
	}
	
	if(!inherits(data, "data.frame")){
	  data <- as.data.frame(data)
	  cat("Input data is required to be of class 'data.frame' - converting now.\n")
	}
  
	if(is.null(weights_vector)){
		weights_vector <- rep(1/length(phy_list), length(phy_list))
	}
	# What should I do if the user a weight vector with support values?
	# The format will be a list of numeric vectors, which will break the rate optimization

	#Ensures that weird root state probabilities that do not sum to 1 are input:
	if(!is.null(root.p)){
		if(!is.character(root.p)){
			root.p <- root.p/sum(root.p)
		}
	}

	input.data <- data
	nCol <- dim(data)[2]
	CorData <- corProcessData(data, collapse = collapse)
	data.legend <- data <- CorData$corData
	nObs <- length(CorData$ObservedTraits)

	for (i in seq_along(phy_list)) {
	phy <- phy_list[[i]]
	if(any(phy$edge.length<=.Machine$double.eps)){
	  warning(paste0("Branch lengths of 0 detected. Adding ", sqrt(.Machine$double.eps)), immediate. = TRUE)
	  phy$edge.length <- phy$edge.length + sqrt(.Machine$double.eps)
	}
	phy_list[[i]] <- phy
	}

	#Creates the data structure and orders the rows to match the tree.
	data.sort <- data.frame(data[,2], data[,2],row.names=data[,1])
	data.sort <- data.sort[phy_list[[1]]$tip.label,]
  levels <- levels(as.factor(data.sort[, 1]))

	if(upper.bound < lower.bound){
	  cat("Your upper bound is smaller than your lower bound.\n")
	}
	lb <- log(lower.bound)
	ub <- log(upper.bound)
	order.test <- TRUE
	set.fog <- FALSE

	obj <- NULL
  	nb.tip <- length(phy_list[[1]]$tip.label)
	nb.node <- phy_list[[1]]$Nnode
	rate.cat <- rate.cat
	root.p <- root.p

	model.set.final <- corHMM:::rate.cat.set.corHMM.JDB(phy=phy_list[[1]], data=input.data, rate.cat=rate.cat, ntraits=nObs, model=model, rate.mat=rate.mat, collapse=collapse)
  
	# this allows for custom rate matricies!
	if(!is.null(rate.mat)){
		order.test <- FALSE
		rate.mat[rate.mat == 0] <- NA
		rate <- rate.mat
		model.set.final$np <- max(rate, na.rm=TRUE)
		rate[is.na(rate)]=max(rate, na.rm=TRUE)+1
		model.set.final$rate <- rate
		model.set.final$index.matrix <- rate.mat
		model.set.final$Q <- matrix(0, dim(rate.mat)[1], dim(rate.mat)[2])
		## for precursor type models ##
		col.sums <- which(colSums(rate.mat, na.rm=TRUE) == 0)
		row.sums <- which(rowSums(rate.mat, na.rm=TRUE) == 0)
		drop.states <- col.sums[which(col.sums == row.sums)]
		if(length(drop.states > 0)){
			model.set.final$liks[,drop.states] <- 0
		}
		###############################
	}
  
	lower = rep(lb, model.set.final$np)
	upper = rep(ub, model.set.final$np)

	if(!is.null(tip.fog)){
		if(rate.cat > 1){
			tip.fog <- rep(tip.fog, rate.cat)
			model.set.final$fog.vec <- tip.fog
		}else{
			model.set.final$fog.vec <- tip.fog
		}
		set.fog <- TRUE
	}

	if(is.null(ip)){
		taxa.missing.data.drop <- which(is.na(data.sort[,1]))
		if(length(taxa.missing.data.drop) != 0){
			tip.labs <- names(taxa.missing.data.drop)
			dat <- as.matrix(data.sort)
			dat.red <- dat[-taxa.missing.data.drop,]
			phy.red <- drop.tip(phy_list[[1]], taxa.missing.data.drop)
			dat.red <- phyDat(dat.red,type="USER", levels=levels)
			phy.tmp <- multi2di(phy.red)
			par.score <- parsimony(phy.tmp, dat.red, method="fitch")/2
		}else{
			dat <- as.matrix(data.sort)
			dat <- phyDat(dat,type="USER", levels=levels)
			phy.tmp <- multi2di(phy_list[[1]])
			par.score <- parsimony(phy.tmp, dat, method="fitch")/2
		}
		tl <- sum(phy_list[[1]]$edge.length)
		mean.change <- par.score/tl
		if(mean.change ==0){
			ip<- rep(0.01 + exp(lb), model.set.final$np)
		}else{
			ip <- rep(mean.change, model.set.final$np)
		}
	}
	
	if(set.fog == TRUE){
		ip <- c(rep(0.01, length(unique(model.set.final$fog.vec))), ip)
		lower <- c(rep(lb, length(unique(model.set.final$fog.vec))), lower)
		upper <- c(rep(log(0.50), length(unique(model.set.final$fog.vec))), upper)
	}

	# the optimization options:
	if(is.null(opts)){
		opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
	}

	# calls the corHMM function:
	out <- nloptr(x0=log(ip), eval_f=corHMM_treeverse_dev, lb=lower, ub=upper, opts=opts, phy_list=phy_list, data=data, weights_vector=weights_vector, rate.cat=rate.cat, rate.mat=rate.mat, model=model, root.p=root.p, tip.fog=tip.fog, collapse=collapse, node.states=node.states)
	
	loglik <- -out$objective
	est.pars <- exp(out$solution)

	if(set.fog == TRUE){
		fog.est <- est.pars[1:length(unique(model.set.final$fog.vec))]
		est.pars <- est.pars[-c(1:length(unique(model.set.final$fog.vec)))]
	}else{
		if(is.numeric(tip.fog)){
			fog.est <- tip.fog
		}else{
			fog.est <- NULL
		}
	}
	
  if (acsr_weights == "same"){
	acsr_weights <- weights_vector
  }
  species.phy<-phy_list[[1]]
  liks.anc<-ancTreeverse(phy_list=phy_list, data=input.data, p=est.pars, weights_vector=acsr_weights, rate.cat=rate.cat, 
            ntraits = NULL, rate.mat = rate.mat, root.p = root.p, model = model, 
            get.tip.states = get.tip.states, collapse = collapse, tip.fog = fog.est)
	pr <- apply(liks.anc$comb.anc.states, 1, which.max)
  species.phy$node.label <- pr
  
	solution <- matrix(est.pars[model.set.final$index.matrix], dim(model.set.final$index.matrix))
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
	rownames(solution) <- colnames(solution) <- StateNames
	if(set.fog == TRUE){
		tip.fog.probs <- numeric(length(model.set.final$fog.vec))
		tip.fog.probs[] <- c(fog.est, 0)[model.set.final$fog.vec]
		names(tip.fog.probs) <- StateNames
		model.set.final$np <- model.set.final$np + length(unique(model.set.final$fog.vec))
		AIC <- -2*loglik+2*model.set.final$np
		AICc <- -2*loglik+(2*model.set.final$np*(nb.tip/(nb.tip-model.set.final$np-1)))
	}else{
		tip.fog.probs <- NULL
		AIC <- -2*loglik+2*model.set.final$np
		AICc <- -2*loglik+(2*model.set.final$np*(nb.tip/(nb.tip-model.set.final$np-1)))
	}
	if (is.character(node.states)) {
		if (node.states == "marginal" || node.states == "scaled"){
			colnames(lik.anc$lik.anc.states) <- StateNames
			colnames(tip.states) <- StateNames
		}
	}
	if(loglik == -1e+06){
	  warning("corHMM may have failed to optimize correctly, consider checking inputs and running again.", immediate. = TRUE)
	}
	
	obj = list(loglik = loglik,
	AIC = AIC,
	AICc = AICc,
	rate.cat=rate.cat,
	solution=solution,
	tip.fog.probs=tip.fog.probs,
	index.mat=model.set.final$index.matrix,
	fog.vec=model.set.final$fog.vec,
	data=input.data,
	data.legend = data.legend,
	phy=phy_list,
  	states=liks.anc$comb.anc.states,
  	ntrees=length(phy_list),
  	counts=liks.anc$counts,
	iterations=out$iterations,
	collapse=collapse,
	root.p=root.p)
	class(obj)<-"corhmm"
	return(obj)
}

ancTreeverse <- function(phy_list, data, p, weights = NULL, rate.cat, ntraits = NULL, rate.mat = NULL, model = "ER", root.p = NULL, get.likelihood = FALSE, get.tip.states = FALSE, tip.fog = NULL, get.info = FALSE, collapse = TRUE) {
  ntrees <- length(phy_list)
# ===== Node Matching =====
for (i in seq_along(phy_list)){
    phy_list[[i]]$node.frame<-node.match(phy_list[[1]], phy_list[[i]])
    }
# ===== Ancestral States Reconstruction =====
for (i in seq_along(phy_list)){
    phy_list[[i]]$probs <- ancRECON(phy_list[[i]], data, p, method="marginal", rate.cat, ntraits, rate.mat,
model, root.p, get.likelihood, get.tip.states,tip.fog, get.info, collapse)
print(paste("Finished reconstructing on tree", i, "of", length(phy_list), "trees."))
}
# ===== Weights Vector Construction =====
nref <- length(phy_list[[1]]$node.frame$reference)
  for (i in seq_len(ntrees)) {
    phy_list[[i]]$node.frame$weights <- numeric(nref)
  }

if (is.null(weights)){
for (i in seq_len(ntrees)){
    phy_list[[i]]$node.frame$weights <- rep(1/ntrees, nref)
} 
} else if (is.numeric(weights)){
    if (length(weights) != ntrees){
      stop("The length of weight_vector must match the number of trees provided.")
    }
    weights <- weights / sum(weights)
    for (i in 1:ntrees){
      phy_list[[i]]$node.frame$weights <- rep(weights[i], nref)
    }

} else if (is.list(weights)) {
  if (length(weights) != ntrees){
    stop("The length of weight_vector list must the match number of trees provided.")
    }
  for (i in seq_len(ntrees)) {
    phy_list[[i]]$node.frame$weights <- numeric(nref)
    ntip <- Ntip(phy_list[[i]])
    for (j in seq_along(phy_list[[i]]$node.frame$reference)){
      node <- phy_list[[i]]$node.frame$coalnode[j]
      if (node == 0){
        phy_list[[i]]$node.frame$weights[j] <- 0
      } else {
        index <- node - ntip
        phy_list[[i]]$node.frame$weights[j] <- weights[[i]][index]
      }
    }
  }

  weight.matrix <- matrix(0, nrow = nref, ncol = ntrees)
  for (i in seq_len(ntrees)) {
    weight.matrix[, i] <- phy_list[[i]]$node.frame$weights
  }
  w_sums <- rowSums(weight.matrix)
  w_sums[w_sums == 0] <- 1
  for (i in seq_len(ntrees)) {
    phy_list[[i]]$node.frame$weights <- phy_list[[i]]$node.frame$weights / w_sums
  }
}
# ===== Combine Ancestral States =====
comb.probs <- matrix(0, nrow = nrow(phy_list[[1]]$probs$lik.anc.states), ncol = ncol(phy_list[[1]]$probs$lik.anc.states))
for (i in seq_len(ntrees)){
    foc.phy<-phy_list[[i]]
    for (j in seq_along(foc.phy$node.frame$reference)){
      n2<- foc.phy$node.frame$coalnode[j]
      if (n2 == 0) next
      n2row<- n2-Ntip(foc.phy)
      y <- foc.phy$probs$lik.anc.states[n2row,]
      weight <- foc.phy$node.frame$weights[j]
      y <- y * weight
      comb.probs[j, ] <- comb.probs[j, ] + y
    }
}
final<-list()
final$trees <- phy_list
final$comb.anc.states <- comb.probs
return(final)
}
