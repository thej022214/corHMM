### corHMM -- Generalized hidden Markov Models
# automatic fitting 
corHMMDredge <- function(phy, data, max.rate.cat, root.p="maddfitz", pen.type = "l1", lambda = 1, drop.par = TRUE, drop.threshold = 1e-7,  info.threshold=2, criterion="AIC", merge.params=TRUE, merge.threshold=0, rate.mat=NULL, node.states = "marginal", fixed.nodes=FALSE, ip=NULL, nstarts=0, n.cores=1, get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = FALSE, lower.bound = 1e-10, upper.bound = 100, opts=NULL, verbose=TRUE, p=NULL, rate.cat=NULL, grad=FALSE){
  
  if((is.null(p) & !is.null(rate.cat))){
    print("A rate category was given without specifying a parameter vector (p)")
    return(NULL)
  }
  if((!is.null(p) & is.null(rate.cat))){
    print("A parameter vector (p) was given without specifying a rate category")
    return(NULL)
  }
  
  # FIXED FIT
  if(!is.null(p) & !is.null(rate.cat)){
    if(verbose){
      cat("Evaluating fixed parameters p =", p, "\n")
    }
    fixd_fit <- corHMMDredgeBase(phy=phy, 
                                 data=data, 
                                 rate.cat=rate.cat, 
                                 root.p=root.p,
                                 pen.type = pen.type, 
                                 lambda = lambda, 
                                 rate.mat=rate.mat, 
                                 node.states = node.states, 
                                 fixed.nodes=fixed.nodes, 
                                 ip=ip, 
                                 nstarts=nstarts, 
                                 n.cores=n.cores, 
                                 get.tip.states = get.tip.states, 
                                 lewis.asc.bias = lewis.asc.bias, 
                                 collapse = collapse, 
                                 lower.bound = lower.bound, 
                                 upper.bound = upper.bound, 
                                 opts=opts, 
                                 p=p,
                                 grad=FALSE)
    return(fixd_fit)
  }
  
  # automatic search starting at rate cat 1
  curr_index_mat <- rate.mat
  fit_set <- list()
  model_improved <- TRUE
  count <- 1
  current_rate_category <- 1
  if(verbose){
    cat("Begining dredge...\n")
  }
  while(model_improved){
    curr_fit <- try(corHMMDredgeBase(phy=phy, 
                                 data=data, 
                                 rate.cat=current_rate_category, 
                                 root.p=root.p,
                                 pen.type = pen.type, 
                                 lambda = lambda, 
                                 rate.mat=curr_index_mat, 
                                 node.states = node.states, 
                                 fixed.nodes=fixed.nodes, 
                                 ip=ip, 
                                 nstarts=nstarts, 
                                 n.cores=n.cores, 
                                 get.tip.states = get.tip.states, 
                                 lewis.asc.bias = lewis.asc.bias, 
                                 collapse = collapse, 
                                 lower.bound = lower.bound, 
                                 upper.bound = upper.bound, 
                                 opts=opts, 
                                 p=NULL,
                                 grad=grad))
    if(inherits(curr_fit, "try-error")){
      warning("Model fitting failed. Stopping dredge.")
      model_improved <- FALSE
      fit_set[[count]] <- curr_fit
      next
    }
    curr_info_criterion <- curr_fit[[criterion]]
    fit_set[[count]] <- curr_fit
    if(verbose){
      cat("\n")
      cat("AIC:", curr_info_criterion)
      cat("\nMapping matrix:\n")
      print(curr_fit$index.mat)
    }
    best_info_criterion <- get_best_info_criterion(fit_set, criterion, current_rate_category)
    model_improved <- diff(c(best_info_criterion, curr_info_criterion)) < info.threshold
    count <- count + 1
    if(model_improved | current_rate_category < max.rate.cat){
      # try dropping pars
      curr_index_mat <- drop_pars(curr_fit, drop.threshold)
      if(is.null(curr_index_mat)){
        # try merging pars
        curr_index_mat <- merge_pars(curr_fit, merge.threshold)
        if(is.null(curr_index_mat)){
          current_rate_category <- current_rate_category + 1
          model_improved <- TRUE
          if(current_rate_category > max.rate.cat){
            model_improved <- FALSE
          }else{
            if(verbose){
              cat("\n", rep("*", count-1), "Continuing dredge", rep("*", count-1), "\n")
            }
          }
        }else{
          if(verbose){
            cat("\n", rep("*", count-1), "Continuing dredge", rep("*", count-1), "\n")
          }
        }
      }else{
        if(verbose){
          cat("\n", rep("*", count-1), "Continuing dredge", rep("*", count-1), "\n")
        }
      }
    }
  }
  if(verbose){
    cat("\nDone.\n")
  }
  class(fit_set) <- "corhmm.dredge"
  return(fit_set)
}

get_best_info_criterion <- function(corhmm.obj.list, criterion, rate.cat){
  if(length(corhmm.obj.list) == 1){
    return(Inf)
  }
  info_criteria <- unlist(lapply(corhmm.obj.list, "[[", criterion))
  rate_cats <- unlist(lapply(corhmm.obj.list, "[[", "rate.cat"))
  if(any(rate_cats == rate.cat)){
    best_info_criterion <- min(info_criteria[rate_cats == rate.cat])
  }else{
    return(Inf)
  }
  return(best_info_criterion)
}

drop_pars <- function(corhmm.obj, drop.threshold){
  Q <- corhmm.obj$solution
  Q[is.na(Q)] <- Inf
  to_drop <- unique(corhmm.obj$index.mat[which(Q < drop.threshold)])
  if(length(to_drop) == 0){
    return(NULL)
  }
  index_mat <- dropStateMatPars(corhmm.obj$index.mat, to_drop)
  return(index_mat)
}

merge_pars <- function(corhmm.obj, merge.threshold){
  index_mat <- corhmm.obj$index.mat
  if(corhmm.obj$rate.cat > 1){
    current_pars <- MatrixToPars(corhmm.obj)
    rate_classes <- paste("R", 1:corhmm.obj$rate.cat, sep = "")
    par_list <- vector("list", corhmm.obj$rate.cat+1)
    index_list <- vector("list", corhmm.obj$rate.cat+1)
    for(i in seq(corhmm.obj$rate.cat)){
      search_string <- paste0(rate_classes[i], " .* -> ", rate_classes[i])
      index_list[[i]] <- grep(search_string, names(current_pars))
      par_list[[i]] <- current_pars[index_list[[i]]]
    }
    index_list[[corhmm.obj$rate.cat+1]] <- (1:length(current_pars))[-unlist(index_list)]
    par_list[[corhmm.obj$rate.cat+1]] <- current_pars[index_list[[corhmm.obj$rate.cat+1]]]
    er_test <- all(unlist(lapply(index_list, function(x) length(x) <= 1)))
    if(er_test){
      return(NULL)
    }
    potential_mergers <- lapply(par_list, 
      function(x) merge_current_pars(x, merge.threshold))
    diffs <- mapply(function(x,y){abs(diff(y[x]))}, x=potential_mergers, y=par_list)
    focal_merger <- index_list[[which.min(diffs)]][potential_mergers[[which.min(diffs)]]]
    index_mat_merged <- equateStateMatPars(index_mat, focal_merger)
  }else{
    current_pars <- MatrixToPars(corhmm.obj)
    focal_merger <- merge_current_pars(current_pars, merge.threshold)
    if(is.null(focal_merger)){
      return(NULL)
    }
    index_mat_merged <- equateStateMatPars(index_mat, focal_merger)
  }
  return(index_mat_merged)
}

merge_current_pars <- function(current_pars, merge.threshold){
  if(length(current_pars) <= 1){
    return(NULL)
  }
  dist_mat <- as.matrix(dist(current_pars))
  dist_mat[upper.tri(dist_mat)] <- 0
  focal_merger <- which(dist_mat == min(dist(current_pars)), arr.ind = TRUE)
  avg_par <- mean(current_pars[focal_merger])
  additional_mergers <- which(abs(current_pars - avg_par) < merge.threshold)
  if(length(additional_mergers) > length(focal_merger)){
    focal_merger <- additional_mergers
  }
  return(focal_merger)
}

# this is the function that does most of the heavy lifting
corHMMDredgeBase <- function(phy, data, rate.cat, root.p="maddfitz", pen.type = "l1", lambda = 1, rate.mat=NULL, node.states = "marginal", fixed.nodes=FALSE, ip=NULL, nstarts=0, n.cores=1, get.tip.states = FALSE,lewis.asc.bias = FALSE, collapse = FALSE, lower.bound = 1e-10, upper.bound = 100, opts=NULL, p=NULL, grad=FALSE){
  
  # Checks to make sure node.states is not NULL.  If it is, just returns a diagnostic message asking for value.
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
  
  if(fixed.nodes == FALSE){
    if(!is.null(phy$node.label)){
      phy$node.label <- NULL
      cat("You specified \'fixed.nodes=FALSE\' but included a phy object with node labels. These node labels have been removed.\n")
    }
  }
  
  #Ensures that weird root state probabilities that do not sum to 1 are input:
  if(!is.null(root.p)){
    if(!is.character(root.p)){
      root.p <- root.p/sum(root.p)
    }
  }
  
  if(pen.type == "unreg"){
    lambda <- 0
  }
  
  # rescale phy to height of one
  phy_original <- phy
  H <- max(node.depth.edgelength(phy))
  phy$edge.length <- phy$edge.length/H
  upper.bound <- upper.bound * H
  lower.bound <- lower.bound / H
  
  input.data <- data
  
  nCol <- dim(data)[2]
  
  CorData <- corProcessData(data, collapse = collapse)
  data.legend <- data <- CorData$corData
  nObs <- length(CorData$ObservedTraits)

    # Checks to make sure phy & data have same taxa. Fixes conflicts (see match.tree.data function).
  matching <- match.tree.data(phy,data)
  data <- matching$data
  phy <- matching$phy
  
  # Will not perform reconstructions on invariant characters (unless rate params have been given!)
  if(nlevels(as.factor(data[,1])) <= 1 & !is.null(p)){
    obj <- NULL
    obj$loglik <- NULL
    obj$diagnostic <- paste("Character is invariant. Analysis stopped.",sep="")
    return(obj)
  } else {
    # Still need to make sure second level isnt just an ambiguity
    lvls <- as.factor(data[,1])
    if(nlevels(as.factor(data[,1])) == 2 && length(which(lvls == "?"))){
      obj <- NULL
      obj$loglik <- NULL
      obj$diagnostic <- paste("Character is invariant. Analysis stopped.",sep="")
      return(obj)
    }
  }
  
  if(any(phy$edge.length<=.Machine$double.eps)){
    warning(paste0("Branch lengths of 0 detected. Adding ", sqrt(.Machine$double.eps)), immediate. = TRUE)
    #   phy$edge.length[phy$edge.length<=1e-5] <- 1e-5
    phy$edge.length <- phy$edge.length + sqrt(.Machine$double.eps) 
  }
  #Creates the data structure and orders the rows to match the tree.
  data.sort <- data.frame(data[,2], data[,2],row.names=data[,1])
  data.sort <- data.sort[phy$tip.label,]
  
  counts <- table(data.sort[,1])
  levels <- levels(as.factor(data.sort[,1]))
  cols <- as.factor(data.sort[,1])

  #Some initial values for use later
  k=2
  if(upper.bound < lower.bound){
    cat("Your upper bound is smaller than your lower bound.\n")
  }
  lb <- log(lower.bound)
  ub <- log(upper.bound)
  order.test <- FALSE
  
  obj <- NULL
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  root.p <- root.p
  nstarts <- nstarts
  ip <- ip
  model = "ARD"
  
    model.set.final <- rate.cat.set.corHMM.JDB(phy=phy,data=input.data,rate.cat=rate.cat,ntraits=nObs,model=model,rate.mat=rate.mat, collapse=collapse)
    # adjusting the matrix for corhmm dredge which allows for independent rate classes
    if(is.null(rate.mat)){
      model.set.final$index.matrix[!is.na(model.set.final$index.matrix)] <- 1:length(model.set.final$index.matrix[!is.na(model.set.final$index.matrix)])
      model.set.final$rate <- model.set.final$index.matrix
      model.set.final$rate[is.na(model.set.final$rate)] <- max(model.set.final$rate, na.rm = TRUE) + 1
      model.set.final$np <- max(model.set.final$index.matrix, na.rm = TRUE)
    }
  phy <- reorder(phy, "pruningwise")
  
  if(collapse){
    StateNames <- gsub("_", "|", CorData$ObservedTraits)
  }else{
    StateNames <- gsub("_", "|", CorData$PossibleTraits)
  }  
  print_counts <- rep(0, length(StateNames))
  if(length(grep("&", names(counts))) > 0){
    counts <- counts[-grep("&", names(counts))]
  }
  print_counts[as.numeric(names(counts))] <- counts
  # cat("State distribution in data:\n")
  # cat("States:",StateNames,"\n",sep="\t")
  # cat("Counts:",print_counts,"\n",sep="\t")
  
  lower = rep(lb, model.set.final$np)
  upper = rep(ub, model.set.final$np)
  
  if(is.null(opts)){
    if(grad){
      opts <- list("algorithm"="NLOPT_LD_MMA", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
    }else{
      opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
    }
  }
  if(!is.null(p)){
    cat("Calculating likelihood from a set of fixed parameters", "\n")
    out<-NULL
    p <- p*H
    est.pars<-log(p)
    out$objective<-dev.corhmm.dredge(est.pars,phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p, rate.cat = rate.cat, order.test = order.test, lewis.asc.bias = lewis.asc.bias, pen.type = pen.type, lambda = lambda)
    loglik <- -out$objective
    est.pars <- exp(est.pars)
  }else{
    if(is.null(ip)){
      #If a user-specified starting value(s) is not supplied this begins loop through a set of randomly chosen starting values:
      #Sets parameter settings for random restarts by taking the parsimony score and dividing
      #by the total length of the tree
      # cat("Beginning thorough optimization search -- performing", nstarts, "random restarts", "\n")
      taxa.missing.data.drop <- which(is.na(data.sort[,1]))
      if(length(taxa.missing.data.drop) != 0){
        tip.labs <- names(taxa.missing.data.drop)
        dat <- as.matrix(data.sort)
        dat.red <- dat[-taxa.missing.data.drop,]
        phy.red <- drop.tip(phy, taxa.missing.data.drop)
        dat.red <- phyDat(dat.red,type="USER", levels=levels)
        phy.tmp <- multi2di(phy.red)
        par.score <- parsimony(phy.tmp, dat.red, method="fitch")/2
      }else{
        dat <- as.matrix(data.sort)
        dat <- phyDat(dat,type="USER", levels=levels)
        phy.tmp <- multi2di(phy)
        par.score <- parsimony(phy.tmp, dat, method="fitch")/2
      }
      tl <- sum(phy$edge.length)
      mean.change = par.score/tl
      
      random.restart<-function(nstarts){
        tmp = matrix(,1,ncol=(1+model.set.final$np))
        if(mean.change==0){
          starts=rep(0.01+exp(lb), model.set.final$np)
        }else{
          starts<-sort(rexp(model.set.final$np, 1/mean.change), decreasing = TRUE)
        }
        starts[starts < exp(lb)] = exp(lb)
        starts[starts > exp(ub)] = exp(lb)
        out = nloptr(x0=log(starts), eval_f=dev.corhmm.dredge, lb=lower, ub=upper, opts=opts, phy=phy, liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p, rate.cat = rate.cat, order.test = order.test, lewis.asc.bias = lewis.asc.bias, pen.type = pen.type, lambda = lambda, grad=grad)
        tmp[,1] = out$objective
        tmp[,2:(model.set.final$np+1)] = out$solution
        tmp
      }
      # this is the first model pass
      if(n.cores > 1){
        restart.set<-mclapply(1:nstarts, random.restart, mc.cores=n.cores)
      }else{
        restart.set<-lapply(1:nstarts, random.restart)
      }
      #Finds the best fit within the restart.set list
      best.fit<-which.min(unlist(lapply(restart.set, function(x) x[1])))
      #Generates an object to store results from restart algorithm:
      out<-NULL
      out$objective=unlist(restart.set[[best.fit]][,1])
      out$solution=unlist(restart.set[[best.fit]][,2:(model.set.final$np+1)])
      loglik <- -out$objective
      est.pars <- exp(out$solution)
    }else{
      # the user has specified initial params
      # cat("Beginning subplex optimization routine -- Starting value(s):", ip, "\n")
      ip=ip
      out = nloptr(x0=rep(log(ip), length.out = model.set.final$np), eval_f=dev.corhmm.dredge, lb=lower, ub=upper, opts=opts, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p, rate.cat = rate.cat, order.test = order.test, lewis.asc.bias = lewis.asc.bias, pen.type = pen.type, lambda = lambda, grad=grad)
      loglik <- -out$objective
      est.pars <- exp(out$solution)
    }
  }
  
  #Starts the ancestral state reconstructions:
  if(node.states != "none") {
    # cat("Finished. Inferring ancestral states using", node.states, "reconstruction.","\n")
  }
  TIPS <- 1:nb.tip
  if (node.states == "marginal" || node.states == "scaled"){
    lik.anc <- ancRECON(phy, input.data, est.pars, rate.cat, rate.mat=model.set.final$rate, method=node.states, ntraits=NULL, root.p=root.p, model = model, get.tip.states = get.tip.states, collapse = collapse)
    pr<-apply(lik.anc$lik.anc.states,1,which.max)
    phy$node.label <- pr
    tip.states <- lik.anc$lik.tip.states
    row.names(tip.states) <- phy$tip.label
  }
  if (node.states == "joint"){
    lik.anc <- ancRECON(phy, input.data, est.pars, rate.cat, rate.mat=model.set.final$rate, method=node.states, ntraits=NULL, root.p=root.p, model = model, get.tip.states = get.tip.states, collapse = collapse)
    phy$node.label <- lik.anc$lik.anc.states
    tip.states <- lik.anc$lik.tip.states
  }
  if (node.states == "none") {
    lik.anc <- list(lik.tip.states=NA, lik.anc.states=NA, info.anc.states=NA)
    phy$node.label <- NA
    tip.states <- NA
  }
  
  # finalize the output
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
  np <- model.set.final$np
  index.matrix <- model.set.final$index.matrix
  rownames(solution) <- colnames(solution) <- StateNames
  rownames(index.matrix) <- colnames(index.matrix) <- StateNames
  # rescale phylogeny
  solution <- solution/H
  solution[solution < lower.bound] <- lower.bound
  
  AIC <- -2*loglik+2*np
  AICc <- -2*loglik+(2*np*(nb.tip/(nb.tip-np-1)))
  
  if (is.character(node.states)) {
    if (node.states == "marginal" || node.states == "scaled"){
      colnames(lik.anc$lik.anc.states) <- StateNames
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
             index.mat=index.matrix,
             data=input.data,
             data.legend = data.legend,
             phy=phy_original,
             states=lik.anc$lik.anc.states,
             tip.states=tip.states,
             states.info = lik.anc$info.anc.states,
             iterations=out$iterations,
             collapse=collapse,
             root.p=root.p,
             pen.type=pen.type,
             lambda=lambda)
  class(obj)<-"corhmm"
  return(obj)
}


### The function used to optimize parameters:

dev.corhmm.dredge <- function(p,phy,liks,Q,rate,root.p,rate.cat,order.test,lewis.asc.bias,pen.type="l1",lambda=1,grad=FALSE){
  p = exp(p)
  cp_root.p <- root.p
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  TIPS <- 1:nb.tip
  comp <- numeric(nb.tip + nb.node)
  #Obtain an object of all the unique ancestors
  anc <- unique(phy$edge[,1])
  k.rates <- dim(Q)[2] / 2
  if (any(is.nan(p)) || any(is.infinite(p))) return(1000000)
  
  Q[] <- c(p, 0)[rate]
  diag(Q) <- -rowSums(Q)
  pen_score <- get_penalty_score(Q, p, pen.type, rate, rate.cat)
  # # if the q matrix has columns not estimated, remove them
  # row2rm <- apply(rate, 1, function(x) all(x == max(rate)))
  # col2rm <- apply(rate, 2, function(x) all(x == max(rate)))
  # Q.root <- Q[!row2rm | !col2rm, !row2rm | !col2rm]
  if(is.character(root.p)){
    if(root.p == "yang"){
      root.test <- Null(Q)
      if(dim(root.test)[2]>1){
        return(1000000)
      }
    }      
  }
  
  if(order.test == TRUE){
    # ensure that the rate classes have mean rates in a consistent order (A > B > C > n)
    StateOrderMat <- matrix(1, (dim(Q)/rate.cat)[1], (dim(Q)/rate.cat)[2])
    RateClassOrderMat <- matrix(0, rate.cat, rate.cat)
    diag(RateClassOrderMat) <- 1:rate.cat
    OrderMat <- RateClassOrderMat %x% StateOrderMat
    Rate01 <- vector("numeric", rate.cat)
    for(i in 1:rate.cat){
      tmp <- Q[OrderMat == i]
      Rate01[i] <- tmp[tmp>=0][1]
    }
    OrderTest <- all.equal(Rate01, sort(Rate01, decreasing = TRUE))
    if(OrderTest != TRUE){
      return(1000000)
    }
  }
  
  for (i in seq(from = 1, length.out = nb.node)) {
    #the ancestral node at row i is called focal
    focal <- anc[i]
    #Get descendant information of focal
    desRows <- which(phy$edge[,1]==focal)
    desNodes <- phy$edge[desRows,2]
    v <- 1
    #Loops through all descendants of focal (how we deal with polytomies):
    for (desIndex in sequence(length(desRows))){
      v <- v*expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks[desNodes[desIndex],]
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
    
    #Sum the likelihoods:
    comp[focal] <- sum(v)
    #Divide each likelihood by the sum to obtain probabilities:
    liks[focal, ] <- v/comp[focal]
  }
  
  #Specifies the root:
  root <- nb.tip + 1L
  #If any of the logs have NAs restart search:
  if (is.na(sum(log(comp[-TIPS])))){return(1000000)}
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
    loglik<- -(sum(log(comp[-TIPS])) + log(sum(flat.root * liks[root,])))
  }
  if(is.character(root.p)){
    # root.p==yang will fix root probabilities based on the inferred rates: q10/(q01+q10)
    if(root.p == "yang"){
      root.p <- Null(Q)
      root.p <- c(root.p/sum(root.p))
      loglik <- -(sum(log(comp[-TIPS])) + log(sum(root.p * liks[root,])))
      if(is.infinite(loglik)){
        return(1000000)
      }
    }else{
      # root.p==maddfitz will fix root probabilities according to FitzJohn et al 2009 Eq. 10:
      root.p = liks[root,] / sum(liks[root,])
      loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,])))))
    }
  }else{
    if(is.numeric(root.p[1])){
      loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,])))))
      if(is.infinite(loglik)){
        return(1000000)
      }
    }
  }
  # root.p!==NULL will fix root probabilities based on user supplied vector:
  if(lewis.asc.bias == TRUE){
    p <- log(p)
    dummy.liks.vec <- getLewisLikelihood(p = p, phy = phy, liks = liks, Q = Q, rate = rate, root.p = cp_root.p, rate.cat = rate.cat)
    loglik <- loglik - log(sum(root.p * (1 - exp(dummy.liks.vec))))
  }
  if(grad){
    epsilon <- 1e-6
    base_log_likelihood <- -(loglik + (pen_score * lambda))
    gradient <- numeric(length(p))
    for (i in seq_along(p)) {
      p_plus <- log(p)
      p_plus[i] <- p_plus[i] + epsilon
      log_likelihood_plus <- -dev.corhmm.dredge(p_plus,phy,liks,Q,rate,root.p,rate.cat,order.test,lewis.asc.bias,pen.type,lambda,grad=FALSE)
      gradient[i] <- (log_likelihood_plus - base_log_likelihood) / epsilon
    }
    cat("\rlnLik:", paste0(base_log_likelihood, "..."))
    # print(-gradient)
    return(list(
      objective = -base_log_likelihood,
      gradient = -gradient))
  }else{
    return(loglik + (pen_score * lambda))
  }
}


### The function used to calculate the penalty:

get_penalty_score <- function(Q, p, pen.type, index.mat, rate.cat){
  if(rate.cat == 1){
    if(pen.type == "l1"){
      pen <- mean(-diag(Q))
    }
    if(pen.type == "l2"){
      pen <- mean(diag(Q)^2)
    }
    if(pen.type == "er"){
      # pen <- sd(-diag(Q))
      if(length(p) < 1){
        pen <- mean(dist(p))
      }else{
        pen <- 0
      }
    }
  }else{
    rate_class_names <- paste0("R", 1:rate.cat)
    pen_by_rc <- numeric(rate.cat)
    diag(Q) <- 0
    for(i in seq_along(rate_class_names)){
      rc_index <- grep(rate_class_names[i], colnames(index.mat))
      if(pen.type == "l1"){
        pen_by_rc[i] <- mean(Q[rc_index, rc_index])
      }
      if(pen.type == "l2"){
        pen_by_rc[i] <- mean(Q[rc_index, rc_index]^2)
      }
      if(pen.type == "er"){
        rates <- Q[rc_index, rc_index]
        rates <- rates[rates>0]
        # pen_by_rc[i] <- sd(rates)
        if(length(rates) < 1){
          pen_by_rc[i] <- mean(dist(rates))
        }else{
          pen_by_rc[i] <- 0
        }
      }
    }
    pen <- sum(pen_by_rc)
  }
  return(pen)
}

### Functions for leave one out cross validation:

get_per_tip_faith_PD <- function(phy, fold_vec){
  fold_hist <- setNames(rep(NA, length(unique(fold_vec))), sort(unique(fold_vec)))
  for(i in 1:length(sort(unique(fold_vec)))){
    focal_fold <- sort(unique(fold_vec))[i]
    focal_tips <- names(fold_vec)[fold_vec == focal_fold]
    focal_phy <- keep.tip(phy, focal_tips)
    fold_hist[i] <- sum(focal_phy$edge.length)/length(focal_tips)
  }
  return(fold_hist)
}

# other weight strategies to test
# split_data_k_other <- function(phy, k=5){
#   weights <- rowSums(vcv.phylo(phy))
#   fold_vec <- setNames(rep(NA, length(phy$tip.label)), phy$tip.label)
#   for(i in 1:(length(fold_vec)-1)){
#     curr_k <- i %% k
#     focal <- sample(names(weights), 1, replace = FALSE, prob = weights)
#     fold_vec[focal] <- curr_k
#     phy <- drop.tip(phy, focal)
#     weights <- rowSums(vcv.phylo(phy))
#   }
#   focal <- sample(names(weights), 1, replace = FALSE, prob = weights)
#   curr_k <- length(weights) %% k
#   fold_vec[focal] <- curr_k
#   return(fold_vec)
# }
# 
# split_k_fold_even <- function(phy, k=5){
#   weights <- setNames(rep(1, length(phy$tip.label)), phy$tip.label)
#   fold_vec <- setNames(rep(NA, length(phy$tip.label)), phy$tip.label)
#   for(i in 1:(length(fold_vec)-1)){
#     curr_k <- i %% k
#     focal <- sample(names(weights), 1, replace = FALSE, prob = weights)
#     fold_vec[focal] <- curr_k
#     weights <- weights[!names(weights) %in% focal]
#   }
#   focal <- sample(names(weights), 1, replace = FALSE, prob = weights)
#   curr_k <- length(weights) %% k
#   fold_vec[focal] <- curr_k
#   return(fold_vec)
# }

get_weights <- function(phy){
  # based on rholf 2001
  C <- vcv.phylo(phy)
  C_inv <- solve(C)
  I <- matrix(0, dim(C)[1], dim(C)[2])
  diag(I) <- 1
  weights <- (t(I) %*% C_inv) %*% matrix(1, dim(C)[1], 1) 
  weights <- c(weights/sum(weights))
  weights <- setNames(weights, colnames(C))
  return(weights)
}

split_data_k_folds <- function(phy, k=5){
  # weights <- rowSums(vcv.phylo(phy))
  weights <- get_weights(phy)
  fold_vec <- setNames(rep(NA, length(phy$tip.label)), phy$tip.label)
  for(i in 1:(length(fold_vec)-1)){
    curr_k <- i %% k
    focal <- sample(names(weights), 1, replace = FALSE, prob = weights)
    fold_vec[focal] <- curr_k
    weights <- weights[!names(weights) %in% focal]
  }
  focal <- sample(names(weights), 1, replace = FALSE, prob = weights)
  curr_k <- length(weights) %% k
  fold_vec[focal] <- curr_k
  return(fold_vec)
}

# Function to perform k-fold cross-validation
kFoldCrossValidation <- function(corhmm_obj, k, lambdas=NULL, return_model=TRUE, save_model_dir=NULL, model_name=NULL) {
  scores <- numeric(k)  # Create an empty vector to store scores for each fold
  folds <- split_data_k_folds(corhmm_obj$phy, k)
  ip <- MatrixToPars(corhmm_obj)
  if(is.null(lambdas)){
    if(is.null(corhmm_obj$lambda)){
      lambdas <- 0
    }else{
      lambdas <- corhmm_obj$lambda
    }
  }
  total_model_list <- vector("list", length(lambdas))
  count <- 1
  for(lambda in lambdas){
    model_list <- vector("list", k)
    cat("Evaluating lambda =", lambda, "\n")  # Print the score for the fold
    for (i in 0:(k-1)) {
      # get the original data
      fold_data <- corhmm_obj$data
      # Split data into training and testing sets for the current fold
      fold_data[folds == i, 2:ncol(fold_data)] <- "?"
      # Train the model on training data
      model <- corHMMDredgeBase(phy = corhmm_obj$phy,
                                         data = fold_data,
                                         rate.cat = corhmm_obj$rate.cat, 
                                         pen.type = corhmm_obj$pen.type, 
                                         lambda = lambda, 
                                         rate.mat = corhmm_obj$index.mat, 
                                         node.states = "marginal", 
                                         fixed.nodes=FALSE, 
                                         root.p=corhmm_obj$root.p, 
                                         ip=ip, 
                                         nstarts=0, 
                                         n.cores=1, 
                                         get.tip.states = TRUE, 
                                         lewis.asc.bias = FALSE, 
                                         collapse = corhmm_obj$collapse, 
                                         lower.bound = 1e-10, 
                                         upper.bound = 100, 
                                         opts=NULL, 
                                         p=NULL,
                                         grad=corhmm_obj$grad)
      
      # Evaluate the model on testing data
      score <- evaluateModel(model, corhmm_obj)
      scores[i+1] <- score  # Store the score for this fold
      if(return_model){
        model_list[[i+1]] <- model # Store the model for this fold
      }else{
        model_list[[i+1]] <- NULL
      }
      if(!is.null(save_model_dir)){
        if(is.null(model_name)){
          model_name <- "corhmm.obj"
        }
        saveRDS(model, file = paste0(save_model_dir, "/", model_name, "_lambda", lambda, "_fold", folds, ".RDS"))
      }
      cat("Fold", i, "Score:", score, "\n")  # Print the score for the fold
    }
    averageScore <- mean(scores)  # Calculate the average score across all folds
    cat("Average Cross-Validation Score:", averageScore, "\n")
    total_model_list[[count]] <- list(models = model_list, scores = scores, averageScore = averageScore)
    count <- count + 1
  }
  names(total_model_list) <- lambdas
  class(total_model_list) <- c("corhmm.kfold")
  return(total_model_list)
}

# Function to evaluate the model on testing data
evaluateModel <- function(model, corhmm_obj){
  tip_liks <- get_tip_liks(corhmm_obj)
  scores <- numeric(dim(tip_liks)[1])
  for(i in 1:nrow(tip_liks)){
    scores[i] <- js_divergence(tip_liks[i,], model$tip.states[i,])
  }
  score <- mean(scores)
  return(score)
}

# Calculate Total Variation Distance
total_variation_distance <- function(P, Q) {
  sum(abs(P - Q)) / 2
}

# Calculate KL Divergence
kl_divergence <- function(P, Q) {
  if (any(P > 0 & Q == 0)) {
    return(Inf)  # To handle the case where Q(i) = 0 and P(i) > 0
  }
  return(sum(P * log(P / Q), na.rm = TRUE))
}

# Calculate Jensen-Shannon Divergence
# This is a symmetric and smoothed version of KL divergence. It is defined as the average of the KL divergences between each distribution and the average of both distributions
js_divergence <- function(P, Q) {
  M <- (P + Q) / 2
  return((kl_divergence(P, M) + kl_divergence(Q, M)) / 2)
}

get_tip_liks <- function(corhmm_obj, return_original=TRUE){
  model.set.final <- get_MSF_from_corhm_obj(corhmm_obj)
  pars <- MatrixToPars(corhmm_obj)
  phy <- corhmm_obj$phy
  phy$node.label <- NULL
  ntips <- Ntip(phy)
  nnodes <- Nnode(phy)
  liks <- model.set.final$liks
  rownames(liks) <- c(phy$tip.label, (ntips+1):(ntips+nnodes))
  tip_liks <- liks[1:ntips,]
  if(return_original){
    return(tip_liks)
  }
  for(i in 1:ntips){
    liks_row <- c()
    liks_copy <- liks
    for(j in 1:ncol(liks)){
      liks_copy[i,] <- 0
      liks_copy[i, j] <- 1
      lik_tmp <- dev.corhmm(log(pars),
                            phy=phy,
                            liks=liks_copy,
                            Q=model.set.final$Q,
                            rate=model.set.final$rate,
                            root.p=corhmm_obj$root.p, 
                            rate.cat = corhmm_obj$rate.cat, 
                            order.test = FALSE, 
                            lewis.asc.bias = FALSE)
      liks_row[j] <- -lik_tmp
    }
    best_probs <- max(liks_row)
    liks_tip_rescaled <- liks_row - best_probs
    tip_liks[i,] <- exp(liks_tip_rescaled) / sum(exp(liks_tip_rescaled))
  }
  return(tip_liks)
}

get_MSF_from_corhm_obj <- function(corhmm_obj){
  model.set.final <- rate.cat.set.corHMM.JDB(
    phy=corhmm_obj$phy,data=corhmm_obj$data,
    rate.cat=corhmm_obj$rate.cat,
    ntraits=dim(corhmm_obj$index.mat)[1]/corhmm_obj$rate.cat,
    model=NULL,
    rate.mat=corhmm_obj$index.mat, 
    collapse=corhmm_obj$collapse)
  return(model.set.final)
}

print.corhmm.kfold <- function(x,...){
  score_table <- do.call(rbind, lapply(x, "[[", "scores"))
  colnames(score_table) <- paste0("Fold:", 0:(dim(score_table)[2]-1))
  rownames(score_table) <- paste0("Lambda:", rownames(score_table))
  score_table <- t(score_table)
  avg_scores <- colMeans(score_table)
  cat("\nScores per fold:\n")
  print(score_table)
  cat("\n")
  cat("Average Scores:\n")
  print(avg_scores)
  cat("\n")
  invisible(x)
}

