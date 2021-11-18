

######################################################################################################################################
######################################################################################################################################
### corHMM -- Generalized hidden Markov Models
######################################################################################################################################
######################################################################################################################################

corHMM <- function(phy, data, rate.cat, rate.mat=NULL, model = "ARD", node.states = "marginal", fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=1, get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9, upper.bound = 100, opts=NULL){
    
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
    # nObs <- length(CorData$ObservedTraits)
    if(length(grep("&", CorData$corData[,2])) > 0){
      non_and_chars <- as.numeric(CorData$corData[,2][-grep("&", CorData$corData[,2])])
      and_chars <- as.numeric(unlist(strsplit(CorData$corData[,2][grep("&", CorData$corData[,2])], "&")))
      nObs <- max(c(non_and_chars, and_chars))
    }else{
      nObs <- max(as.numeric(CorData$corData[,2]))
    }
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
    
    if(any(phy$edge.length<=1e-5)){
      warning("Branch lengths of 0 detected. Adding 1e-5 to these branches.", immediate. = TRUE)
      phy$edge.length[phy$edge.length<=1e-5] <- 1e-5
    }
    #Creates the data structure and orders the rows to match the tree.
    data.sort <- data.frame(data[,2], data[,2],row.names=data[,1])
    data.sort <- data.sort[phy$tip.label,]
    
    counts <- table(data.sort[,1])
    levels <- levels(as.factor(data.sort[,1]))
    cols <- as.factor(data.sort[,1])
    cat("State distribution in data:\n")
    cat("States:",levels,"\n",sep="\t")
    cat("Counts:",counts,"\n",sep="\t")
    
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
    rate.cat <- rate.cat
    root.p <- root.p
    nstarts <- nstarts
    ip <- ip
    
    model.set.final <- rate.cat.set.corHMM.JDB(phy=phy,data=input.data,rate.cat=rate.cat,ntraits=nObs,model=model,rate.mat=rate.mat, collapse=collapse)
    phy <- reorder(phy, "pruningwise")
    
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
    
    if(is.null(opts)){
      opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
    }
    if(!is.null(p)){
        cat("Calculating likelihood from a set of fixed parameters", "\n")
        out<-NULL
        est.pars<-log(p)
        out$objective<-dev.corhmm(est.pars,phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p, rate.cat = rate.cat, order.test = order.test, lewis.asc.bias = lewis.asc.bias)
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
                out = nloptr(x0=log(starts), eval_f=dev.corhmm, lb=lower, ub=upper, opts=opts, phy=phy, liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p, rate.cat = rate.cat, order.test = order.test, lewis.asc.bias = lewis.asc.bias)
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
            out = nloptr(x0=rep(log(ip), length.out = model.set.final$np), eval_f=dev.corhmm, lb=lower, ub=upper, opts=opts, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p, rate.cat = rate.cat, order.test = order.test, lewis.asc.bias = lewis.asc.bias)
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
    StateNames <- paste("(", rep(1:(dim(model.set.final$index.matrix)[1]/rate.cat), rate.cat), ",", rep(paste("R", 1:rate.cat, sep = ""), each = nObs), ")", sep = "")
    rownames(solution) <- colnames(solution) <- StateNames
    AIC <- -2*loglik+2*model.set.final$np
    AICc <- -2*loglik+(2*model.set.final$np*(nb.tip/(nb.tip-model.set.final$np-1)))
    
    if (is.character(node.states)) {
        if (node.states == "marginal" || node.states == "scaled"){
            colnames(lik.anc$lik.anc.states) <- StateNames
        }
    }
    
    obj = list(loglik = loglik,
    AIC = AIC,
    AICc = AICc,
    rate.cat=rate.cat,
    solution=solution,
    index.mat=model.set.final$index.matrix,
    data=input.data,
    data.legend = data.legend,
    phy=phy,
    states=lik.anc$lik.anc.states,
    tip.states=tip.states,
    states.info = lik.anc$info.anc.states,
    iterations=out$iterations,
    root.p=root.p)
    class(obj)<-"corhmm"
    return(obj)
}


######################################################################################################################################
######################################################################################################################################
### The function used to optimize parameters:
######################################################################################################################################
######################################################################################################################################

dev.corhmm <- function(p,phy,liks,Q,rate,root.p,rate.cat,order.test,lewis.asc.bias) {
  
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
  }
  # root.p!==NULL will fix root probabilities based on user supplied vector:
  if(is.numeric(root.p[1])){
      loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,])))))
      if(is.infinite(loglik)){
          return(1000000)
      }
  }
  if(lewis.asc.bias == TRUE){
    p <- log(p)
    dummy.liks.vec <- getLewisLikelihood(p = p, phy = phy, liks = liks, Q = Q, rate = rate, root.p = cp_root.p, rate.cat = rate.cat)
    loglik <- loglik - log(sum(root.p * (1 - exp(dummy.liks.vec))))
  }
  return(loglik)
}




######################################################################################################################################
######################################################################################################################################
### The various utility functions used
######################################################################################################################################
######################################################################################################################################

getLewisLikelihood <- function(p, phy, liks, Q, rate, root.p, rate.cat){
  nTips <- length(phy$tip.label)
  nNodes <- length(unique(phy$edge[,1]))
  nStates <- dim(Q)[1]/rate.cat
  states_structure <- seq(from = 1, by = nStates, length.out = rate.cat)
  dummy.liks.vec <- vector("numeric", nStates)
  for(state_i in 1:nStates){
    lik_structure_i <- states_structure + (state_i - 1) 
    liks[] <- 0
    liks[1:nTips, lik_structure_i] <- 1
    dummy.liks.vec[state_i] <- -dev.corhmm(p = p, phy = phy, liks = liks, Q = Q, rate = rate, root.p = root.p, rate.cat = rate.cat, order.test = FALSE, lewis.asc.bias = FALSE)
  }
  return(dummy.liks.vec)
}

# JDB modified functions
rate.cat.set.corHMM.JDB<-function(phy,data,rate.cat, ntraits, model, rate.mat=NULL, collapse=TRUE){
    
    obj <- NULL
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    obj$rate.cat<-rate.cat
    if(is.null(rate.mat)){
      rate <- getStateMat4Dat(data, model)$rate.mat
      if(rate.cat > 1){
        StateMats <- vector("list", rate.cat)
        for(i in 1:rate.cat){
          StateMats[[i]] <- rate
        }
        rate <- getFullMat(StateMats)
      }
    }else{
      rate <- rate.mat
      ntraits <- dim(rate)[1]/rate.cat
    }
    nTraits <- dim(rate)[1]
    rate[rate == 0] <- NA
    index.matrix<-rate
    rate[is.na(rate)]<-max(rate,na.rm=TRUE)+1
    
    CorData <- corProcessData(data, collapse = collapse)
    data <- CorData$corData
    nObs <- length(CorData$ObservedTraits)
    
    matching <- match.tree.data(phy,data)
    data <- matching$data
    
    # this is no longer needed since corProcessData will produce a dataset of a specific type every time
    # x <- as.numeric(data.sort[,1])
    # TIPS <- 1:nb.tip
    # if(min(x) !=0){
    #     x <- x - min(x)
    # }
    # this is being removed temporarily - this handles NA tip states
    # for(i in 1:nb.tip){
    #     if(is.na(x[i])){x[i]=2}
    # }
    
    tmp <- matrix(0, nb.tip + nb.node, ntraits)
    for(i in 1:nb.tip){
        state_index <- as.numeric(unlist(strsplit(as.character(data[i,2]), "&")))
        tmp[i, state_index] <- 1
    }
    liks <- matrix(rep(tmp, rate.cat), nb.tip + nb.node, ntraits*rate.cat)
    
    Q <- matrix(0, dim(rate)[1], dim(rate)[1])
    
    obj$np<-max(rate)-1
    obj$rate<-rate
    obj$index.matrix<-index.matrix
    obj$liks<-liks
    obj$Q<-Q
    
    return(obj)
}

corProcessData <- function(data, collapse=TRUE){
  nCol <- dim(data)[2]
  LevelList <- StateMats <- vector("list", nCol-1)
  # detect the number of states in each column. & is treated as indicating polymorphism. ? is treated as unknown data.
  for(i in 2:nCol){
    data[,i] <- as.character(data[,i])
    data_i <- as.character(data[,i])
    data_i <- data_i[!data_i == "?"]
    States_i <- unique(unlist(strsplit(data_i, "&")))
    StateMats[[i-1]] <- getRateCatMat(length(States_i))
    LevelList[[i-1]] <- sort(States_i)
  }
  # identify the possible trait combinations
  TraitList <- expand.grid(LevelList)
  Traits <- sort(apply(TraitList, 1, function(x) paste(c(x), collapse = "_")))
  # convert each column into a numeric value associated with a member of the trait combinations. ? are associated with all values of that column, & indicates the combination of two or more
  search.strings <- observed.traits_index <- combined.data <- c()
  for(i in 1:dim(data)[1]){
    data_rowi <- data[i,2:nCol]
    # and symbolizes it can be any of the separated states
    search.string_i <- paste("^",paste(sapply(data_rowi, function(x) paste("(", gsub("&", "|", x), ")", sep = "")),collapse = "_"), "$", sep="")
    # ? means it can be any of the states in that character
    search.string_i <- gsub("(?)", ".*", search.string_i, fixed=TRUE)
    # if the data is polymorphic it will now have ands separating the corHMM states
    combined.data[i] <- paste(grep(search.string_i, Traits), collapse="&")
    observed.traits_index <- c(observed.traits_index, grep(search.string_i, Traits))
    search.strings[i] <- search.string_i
  }
  ObservedTraits <- sort(Traits[unique(observed.traits_index)])
  if(collapse){
    corData <- data.frame(sp = data[,1], d = sapply(search.strings, function(x) paste(grep(x, ObservedTraits), collapse="&")))
  }else{
    corData <- data.frame(sp = data[, 1], d = sapply(search.strings, function(x) paste(grep(x, Traits),collapse = "&")))
  }
  return(list(StateMats = StateMats,  PossibleTraits = Traits, ObservedTraits = ObservedTraits, corData = corData))
}

print.corhmm<-function(x,...){
    
    ntips=Ntip(x$phy)
    output<-data.frame(x$loglik,x$AIC,x$AICc,x$rate.cat,ntips, row.names="")
    names(output)<-c("-lnL","AIC","AICc","Rate.cat","ntax")
    cat("\nFit\n")
    print(output)
    cat("\n")
    
    UserStates <- corProcessData(x$data)$ObservedTraits
    names(UserStates) <- sort(unique(as.numeric(x$data.legend[,2])))
    cat("Legend\n")
    print(UserStates)
    cat("\n")
    
    param.est<- x$solution
    cat("Rates\n")
    print(param.est)
    cat("\n")
    
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


######################################################################################################################################
######################################################################################################################################
### ORIGINAL CODE THAT IS NO LONGER USED
######################################################################################################################################
######################################################################################################################################

#Generalized ace() function that allows analysis to be carried out when there are polytomies:
dev.corhmm.ORIGINAL <- function(p,phy,liks,Q,rate,root.p,rate.cat,order.test) {
    p = exp(p)
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    TIPS <- 1:nb.tip
    comp <- numeric(nb.tip + nb.node)
    phy <- reorder(phy, "pruningwise")
    #Obtain an object of all the unique ancestors
    anc <- unique(phy$edge[,1])
    k.rates <- dim(Q)[2] / 2
    if (any(is.nan(p)) || any(is.infinite(p))) return(1000000)
    
    Q[] <- c(p, 0)[rate]
    diag(Q) <- -rowSums(Q)
    for (i  in seq(from = 1, length.out = nb.node)) {
        #the ancestral node at row i is called focal
        focal <- anc[i]
        #Get descendant information of focal
        desRows<-which(phy$edge[,1]==focal)
        desNodes<-phy$edge[desRows,2]
        v <- 1
        #Loops through all descendants of focal (how we deal with polytomies):
        for (desIndex in sequence(length(desRows))){
            v<-v*expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks[desNodes[desIndex],]
        }
        #Sum the likelihoods:
        comp[focal] <- sum(v)
        #Divide each likelihood by the sum to obtain probabilities:
        liks[focal, ] <- v/comp[focal]
    }
    #Temporary solution for ensuring an ordered Q with respect to the rate classes. If a simpler model is called this feature is automatically turned off:
    par.order<-NA
    if(k.rates == 2){
        try(par.order <- p[3] > p[8])
        if(!is.na(par.order)){
            if(par.order == TRUE){
                return(1000000)
            }
        }
    }
    if(k.rates == 3){
        try(par.order <- p[3] > p[9] | p[9] > p[14])
        if(!is.na(par.order)){
            if(par.order == TRUE){
                return(1000000)
            }
        }
    }
    if(k.rates == 4){
        try(par.order <- p[3] > p[9] | p[9] > p[15] | p[15] > p[20])
        if(!is.na(par.order)){
            if(par.order == TRUE){
                return(1000000)
            }
        }
    }
    if(k.rates == 5){
        try(par.order <- p[3] > p[9] | p[9] > p[15] | p[15] > p[21] | p[21] > p[26])
        if(!is.na(par.order)){
            if(par.order == TRUE){
                return(1000000)
            }
        }
    }
    if(k.rates == 6){
        try(par.order <- p[3] > p[9] | p[9] > p[15] | p[15] > p[21] | p[21] > p[27] | p[27] > p[32])
        if(!is.na(par.order)){
            if(par.order == TRUE){
                return(1000000)
            }
        }
    }
    if(k.rates == 7){
        try(par.order <- p[3] > p[9] | p[9] > p[15] | p[15] > p[21] | p[21] > p[27] | p[27] > p[33] | p[33] > p[38])
        if(!is.na(par.order)){
            if(par.order == TRUE){
                return(1000000)
            }
        }
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
            loglik<- -(sum(log(comp[-TIPS])) + log(sum(flat.root * liks[root,])))
        }
        else{
            if(is.character(root.p)){
                # root.p==yang will fix root probabilities based on the inferred rates: q10/(q01+q10)
                if(root.p == "yang"){
                    root.p <- Null(Q)
                    root.p <- c(root.p/sum(root.p))
                    loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,])))))
                    if(is.infinite(loglik)){
                        return(1000000)
                    }
                }else{
                    # root.p==maddfitz will fix root probabilities according to FitzJohn et al 2009 Eq. 10:
                    root.p = liks[root,] / sum(liks[root,])
                    loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,])))))
                }
            }
            # root.p!==NULL will fix root probabilities based on user supplied vector:
            else{
                loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,])))))
                if(is.infinite(loglik)){
                    return(1000000)
                }
            }
        }
    }
    loglik
}

rate.cat.set.corHMM <- function (phy, data.sort, rate.cat) {
    k = 2
    obj <- NULL
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    obj$rate.cat <- rate.cat
    rate <- rate.mat.maker(hrm = TRUE, rate.cat = rate.cat)
    index.matrix <- rate
    rate[is.na(rate)] <- max(rate, na.rm = TRUE) + 1
    x <- data.sort[, 1]
    TIPS <- 1:nb.tip
    for (i in 1:nb.tip) {
        if (is.na(x[i])) {
            x[i] = 2
        }
    }
    if (rate.cat == 1) {
        liks <- matrix(0, nb.tip + nb.node, k * rate.cat)
        TIPS <- 1:nb.tip
        for (i in 1:nb.tip) {
            if (x[i] == 0) {
                liks[i, 1] = 1
            }
            if (x[i] == 1) {
                liks[i, 2] = 1
            }
            if (x[i] == 2) {
                liks[i, 1:2] = 1
            }
        }
        Q <- matrix(0, k * rate.cat, k * rate.cat)
    }
    if (rate.cat == 2) {
        liks <- matrix(0, nb.tip + nb.node, k * rate.cat)
        for (i in 1:nb.tip) {
            if (x[i] == 0) {
                liks[i, c(1, 3)] = 1
            }
            if (x[i] == 1) {
                liks[i, c(2, 4)] = 1
            }
            if (x[i] == 2) {
                liks[i, 1:4] = 1
            }
        }
        Q <- matrix(0, k * rate.cat, k * rate.cat)
    }
    if (rate.cat == 3) {
        liks <- matrix(0, nb.tip + nb.node, k * rate.cat)
        for (i in 1:nb.tip) {
            if (x[i] == 0) {
                liks[i, c(1, 3, 5)] = 1
            }
            if (x[i] == 1) {
                liks[i, c(2, 4, 6)] = 1
            }
            if (x[i] == 2) {
                liks[i, 1:6] = 1
            }
        }
        Q <- matrix(0, k * rate.cat, k * rate.cat)
    }
    if (rate.cat == 4) {
        liks <- matrix(0, nb.tip + nb.node, k * rate.cat)
        for (i in 1:nb.tip) {
            if (x[i] == 0) {
                liks[i, c(1, 3, 5, 7)] = 1
            }
            if (x[i] == 1) {
                liks[i, c(2, 4, 6, 8)] = 1
            }
            if (x[i] == 2) {
                liks[i, 1:8] = 1
            }
        }
        Q <- matrix(0, k * rate.cat, k * rate.cat)
    }
    if (rate.cat == 5) {
        liks <- matrix(0, nb.tip + nb.node, k * rate.cat)
        for (i in 1:nb.tip) {
            if (x[i] == 0) {
                liks[i, c(1, 3, 5, 7, 9)] = 1
            }
            if (x[i] == 1) {
                liks[i, c(2, 4, 6, 8, 10)] = 1
            }
            if (x[i] == 2) {
                liks[i, 1:10] = 1
            }
        }
        Q <- matrix(0, k * rate.cat, k * rate.cat)
    }
    if (rate.cat == 6) {
        liks <- matrix(0, nb.tip + nb.node, k * rate.cat)
        for (i in 1:nb.tip) {
            if (x[i] == 0) {
                liks[i, c(1, 3, 5, 7, 9, 11)] = 1
            }
            if (x[i] == 1) {
                liks[i, c(2, 4, 6, 8, 10, 12)] = 1
            }
            if (x[i] == 2) {
                liks[i, 1:12] = 1
            }
        }
        Q <- matrix(0, k * rate.cat, k * rate.cat)
    }
    if (rate.cat == 7) {
        liks <- matrix(0, nb.tip + nb.node, k * rate.cat)
        for (i in 1:nb.tip) {
            if (x[i] == 0) {
                liks[i, c(1, 3, 5, 7, 9, 11, 13)] = 1
            }
            if (x[i] == 1) {
                liks[i, c(2, 4, 6, 8, 10, 12, 14)] = 1
            }
            if (x[i] == 2) {
                liks[i, 1:14] = 1
            }
        }
        Q <- matrix(0, k * rate.cat, k * rate.cat)
    }
    obj$np <- max(rate) - 1
    obj$rate <- rate
    obj$index.matrix <- index.matrix
    obj$liks <- liks
    obj$Q <- Q
    obj
}


######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
