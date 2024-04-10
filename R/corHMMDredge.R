######################################################################################################################################
######################################################################################################################################
### corHMM -- Generalized hidden Markov Models
######################################################################################################################################
######################################################################################################################################

corHMMDredge <- function(phy, data, max.rate.cat, pen_type = "l1", lambda = 1, node.states = "marginal", fixed.nodes=FALSE, root.p="maddfitz", drop.par = FALSE, drop.tol = 1e-9, ip=NULL, nstarts=0, n.cores=1, get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-10, upper.bound = 100, opts=NULL, p=NULL){
  
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
    phy$edge.length <- phy$edge.length + sqrt(.Machine$double.eps) # changed to add 1e-5 based on suggestion from Hedvig Skirgård (github issue #27)
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
  order.test <- TRUE
  
  obj <- NULL
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  rate.cat <- max.rate.cat
  root.p <- root.p
  nstarts <- nstarts
  ip <- ip
  rate.mat <- NULL
  model = "ARD"
  
    model.set.final <- rate.cat.set.corHMM.JDB(phy=phy,data=input.data,rate.cat=rate.cat,ntraits=nObs,model=model,rate.mat=rate.mat, collapse=collapse)
    # adjusting the matrix for corhmm dredge which allows for independent rate classes
    model.set.final$index.matrix[!is.na(model.set.final$index.matrix)] <- 1:length(model.set.final$index.matrix[!is.na(model.set.final$index.matrix)])
    model.set.final$rate <- model.set.final$index.matrix
    model.set.final$rate[is.na(model.set.final$rate)] <- max(model.set.final$rate, na.rm = TRUE) + 1
    model.set.final$np <- max(model.set.final$index.matrix, na.rm = TRUE)
    
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
  cat("State distribution in data:\n")
  cat("States:",StateNames,"\n",sep="\t")
  cat("Counts:",print_counts,"\n",sep="\t")
  
  lower = rep(lb, model.set.final$np)
  upper = rep(ub, model.set.final$np)
  
  if(is.null(opts)){
    opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
  }
  if(!is.null(p)){
    cat("Calculating likelihood from a set of fixed parameters", "\n")
    out<-NULL
    est.pars<-log(p)
    out$objective<-dev.corhmm.dredge(est.pars,phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p, rate.cat = rate.cat, order.test = order.test, lewis.asc.bias = lewis.asc.bias, pen_type = pen_type, lambda = lambda)
    loglik <- -out$objective
    est.pars <- exp(est.pars)
  }else{
    if(is.null(ip)){
      #If a user-specified starting value(s) is not supplied this begins loop through a set of randomly chosen starting values:
      #Sets parameter settings for random restarts by taking the parsimony score and dividing
      #by the total length of the tree
      cat("Beginning thorough optimization search -- performing", nstarts, "random restarts", "\n")
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
        out = nloptr(x0=log(starts), eval_f=dev.corhmm.dredge, lb=lower, ub=upper, opts=opts, phy=phy, liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p, rate.cat = rate.cat, order.test = order.test, lewis.asc.bias = lewis.asc.bias, pen_type = pen_type, lambda = lambda)
        tmp[,1] = out$objective
        tmp[,2:(model.set.final$np+1)] = out$solution
        tmp
      }
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
      cat("Beginning subplex optimization routine -- Starting value(s):", ip, "\n")
      ip=ip
      out = nloptr(x0=rep(log(ip), length.out = model.set.final$np), eval_f=dev.corhmm, lb=lower, ub=upper, opts=opts, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p, rate.cat = rate.cat, order.test = order.test, lewis.asc.bias = lewis.asc.bias, pen_type = pen_type, lambda = lambda)
      loglik <- -out$objective
      est.pars <- exp(out$solution)
    }
  }
  
  #Starts the ancestral state reconstructions:
  if(node.states != "none") {
    cat("Finished. Inferring ancestral states using", node.states, "reconstruction.","\n")
  }
  TIPS <- 1:nb.tip
  if (node.states == "marginal" || node.states == "scaled"){
    lik.anc <- ancRECON(phy, input.data, est.pars, rate.cat, rate.mat=rate.mat, method=node.states, ntraits=NULL, root.p=root.p, model = model, get.tip.states = get.tip.states, collapse = collapse)
    pr<-apply(lik.anc$lik.anc.states,1,which.max)
    phy$node.label <- pr
    tip.states <- lik.anc$lik.tip.states
    row.names(tip.states) <- phy$tip.label
  }
  if (node.states == "joint"){
    lik.anc <- ancRECON(phy, input.data, est.pars, rate.cat, rate.mat=rate.mat, method=node.states, ntraits=NULL, root.p=root.p, model = model, get.tip.states = get.tip.states, collapse = collapse)
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
  # return new reduced matrix
  if(drop.par){
    solution[solution < drop.tol] <- NA
    np <- length(solution[!is.na(solution)])
    index.matrix <- solution
    index.matrix[!is.na(index.matrix)] <- 1:np
  }else{
    np <- model.set.final$np
    index.matrix <- model.set.final$index.matrix
  }
  
  rownames(solution) <- colnames(solution) <- StateNames
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
             phy=phy,
             states=lik.anc$lik.anc.states,
             tip.states=tip.states,
             states.info = lik.anc$info.anc.states,
             iterations=out$iterations,
             collapse=collapse,
             root.p=root.p,
             pen_type=pen_type,
             lambda=lambda)
  class(obj)<-"corhmm"
  return(obj)
}


######################################################################################################################################
######################################################################################################################################
### The function used to optimize parameters:
######################################################################################################################################
######################################################################################################################################

dev.corhmm.dredge <- function(p,phy,liks,Q,rate,root.p,rate.cat,order.test,lewis.asc.bias,pen_type="logl1",lambda=1) {
  p = exp(p)
  pen_score <- get_penalty_score(p, pen_type)
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
  return(loglik + (pen_score * lambda))
}


######################################################################################################################################
######################################################################################################################################
### The function used to calculate the penalty:
######################################################################################################################################
######################################################################################################################################

get_penalty_score <- function(p, pen_type){
  if(pen_type == "l1"){
    pen <- sum(p)
  }
  if(pen_type == "l2"){
    pen <- sum(p^2)
  }
  return(pen)
}
