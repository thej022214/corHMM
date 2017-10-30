#############################
#EVOLUTION OF DISCRETE TRAITS, ALLOWING POLYMORPHIC AND MISSING STATES
#library(expm)
#library(phangorn)
#############################

#written by Jeremy M. Beaulieu & Jeffrey C. Oliver

rayDISC<-function(phy,data, ntraits=1, charnum=1, rate.mat=NULL, model=c("ER","SYM","ARD"), node.states=c("joint", "marginal", "scaled"), state.recon=c("subsequently"), lewis.asc.bias=FALSE, p=NULL, root.p=NULL, ip=NULL, lb=0, ub=100, verbose=TRUE, diagn=FALSE){

	# Checks to make sure node.states is not NULL.  If it is, just returns a diagnostic message asking for value.
	if(is.null(node.states)){
		obj <- NULL
		obj$loglik <- NULL
		obj$diagnostic <- paste("No model for ancestral states selected.  Please pass one of the following to rayDISC command for parameter \'node.states\': joint, marginal, or scaled.")
		return(obj)
	}
	else { # even if node.states is not NULL, need to make sure its one of the three valid options
		valid.models <- c("joint", "marginal", "scaled")
		if(!any(valid.models == node.states)){
			obj <- NULL
			obj$loglik <- NULL
			obj$diagnostic <- paste("\'",node.states, "\' is not valid for ancestral state reconstruction method.  Please pass one of the following to rayDISC command for parameter \'node.states\': joint, marginal, or scaled.",sep="")
			return(obj)
		}
		if(length(node.states) > 1){ # User did not enter a value, so just pick marginal.
			node.states <- "marginal"
			cat("No model selected for \'node.states\'. Will perform marginal ancestral state estimation.\n")
		}
	}

    if(!state.recon == "subsequently" & node.states == "marginal" | node.states == "scaled"){
        stop("Simultaneous estimation of rates and states using either marginal or scaled probabilities not yet implemented.", call.=FALSE)
    }

    if(!state.recon == "subsequently"){
        if(!is.null(phy$node.label)){
            if(!is.na(phy$node.label[Ntip(phy)+1])){
                #Checking that the root prior is not set twice.
                root.p <- NULL
            }
        }
    }

    if(!state.recon == "estimate" & !state.recon == "given" & !state.recon == "subsequently"){
        stop("Check that you have a supported state.recon analysis. Options are subsequently, estimate, or given.", call.=FALSE)
    }

    if(state.recon == "subsequently"){
        phy$node.label <- NULL
    }else{
        if(state.recon == "given"){
            if(is.null(phy$node.label)){
                stop("You specified you wanted to estimate rates on a given character history, but the tree does not contain node labels.", call.=FALSE)
            }else{
                if(any(is.na(phy$node.label))){
                    cat("Model will assume you want to estimate rates and states, but include state constraints on some but not all nodes.\n")
                }else{
                    cat("Model will assume you want to estimate rates, but include state constraints all nodes.\n")
                }
            }
        }else{
            if(is.null(phy$node.label)){
                cat("Model will assume you want to estimate rates and states simultaneously.\n")
            }else{
                cat("Model will assume you want to estimate rates and states, but include state constraints on some but not all nodes.\n")
                state.recon="given"
            }
        }
    }

    #Ensures that weird root state probabilities that do not sum to 1 are input:
    if(!is.null(root.p)){
        if(!is.character(root.p)){
            root.p <- root.p/sum(root.p)
        }
    }

	#Creates the data structure and orders the rows to match the tree
	phy$edge.length[phy$edge.length==0]=1e-5

	# Checks to make sure phy & data have same taxa.  Fixes conflicts (see match.tree.data function).
	matching <- match.tree.data(phy,data)
	data <- matching$data
	phy <- matching$phy

	# For character of interest, go ahead and convert any "?" to NA
	data[(data[, charnum + 1] == "?"), charnum + 1] <- NA

	# Wont perform reconstructions on invariant characters -- why not? Seems like you should be able to.
	if(nlevels(as.factor(data[,charnum+1])) <= 1){
		obj <- NULL
		obj$loglik <- NULL
		obj$diagnostic <- paste("Character ",charnum," is invariant. Analysis stopped.",sep="")
		return(obj)
	} else {
		# Still need to make sure second level isnt just an ambiguity
		lvls <- as.factor(data[,charnum+1])
		if(nlevels(as.factor(data[,charnum+1])) == 2 && any(lvls %in% c("?", "NA"))){
			obj <- NULL
			obj$loglik <- NULL
			obj$diagnostic <- paste("Character ",charnum," is invariant. Analysis stopped.",sep="")
			return(obj)
		}
	}

	data.sort <- data.frame(data[,charnum+1],data[,charnum+1],row.names=data[,1]) # added character twice, because at least two columns are necessary
	data.sort <- data.sort[phy$tip.label,] # this might have already been done by match.tree.data

	counts <- table(data.sort[,1])
	levels <- levels(as.factor(data.sort[,1]))
	cols <- as.factor(data.sort[,1])
    if(verbose == TRUE){
        cat("State distribution in data:\n")
        cat("States:",levels,"\n",sep="\t")
        cat("Counts:",counts,"\n",sep="\t")
    }

	#Some initial values for use later - will clean up
	k <- 1 # Only one trait allowed
	factored <- factorData(data.sort,charnum=charnum) # just factoring to figure out how many levels (i.e. number of states) in data.
    nl <- ncol(factored)
	state.names <- colnames(factored) # for subsequent reporting
	bound.hit <- FALSE # to keep track of whether min.rate is one of the rate estimates (and thus, potentially a non-optimal rate)
	# Check to make sure values are reasonable (i.e. non-negative)
	if(ub < 0){
		ub <- 100
	}
	if(lb < 0){
		lb <- 0
	}
	if(ub < lb){ # This user really needs help
		ub <- 100
		lb <- 0
	}

	obj <- NULL
	nb.tip<-length(phy$tip.label)
	nb.node <- phy$Nnode

	model=model
	root.p=root.p
	ip=ip
	model.set.final<-rate.cat.set.rayDISC(phy=phy,data=data.sort,model=model,charnum=charnum)
    if(!is.null(rate.mat)){
		rate <- rate.mat
		model.set.final$np <- max(rate, na.rm=TRUE)
		rate[is.na(rate)]=max(rate, na.rm=TRUE)+1
		model.set.final$rate <- rate
		model.set.final$index.matrix <- rate.mat
	}
	lower = rep(lb, model.set.final$np)
	upper = rep(ub, model.set.final$np)

	opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
	if(!is.null(p)){
        if(verbose == TRUE){
            cat("Calculating likelihood from a set of fixed parameters", "\n")
        }
		out<-NULL
		out$solution<-p
        if(state.recon=="subsequently") {
            out$objective <- dev.raydisc(out$solution,phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p, lewis.asc.bias=lewis.asc.bias)
            loglik <- -out$objective
        } else {
            if(lewis.asc.bias == TRUE){
                loglik.num <- ancRECON(phy=phy, data=data, p=p,  hrm=FALSE, rate.cat=NULL, rate.mat=rate.mat, ntraits=ntraits, method=node.states, model=model, charnum=charnum, root.p=root.p, get.likelihood=TRUE)
                phy.dummy <- phy
                data.dummy <- cbind(phy$tip.label, 0)
                phy.dummy$node.label <- rep(1, length(phy.dummy$node.label))
                loglik.dummy <- ancRECON(phy=phy, data=data, p=p, hrm=FALSE, rate.cat=NULL, rate.mat=rate.mat, ntraits=ntraits, method=node.states, model=model, charnum=charnum, root.p=root.p, get.likelihood=TRUE)
                loglik <- (loglik.num  - log(1 - exp(loglik.dummy)))
                loglik <- out$objective
            }else{
                out$objective <- ancRECON(phy=phy, data=data, p=p, hrm=FALSE, rate.cat=NULL, rate.mat=rate.mat, ntraits=ntraits, method=node.states, model=model, charnum=charnum, root.p=root.p, get.likelihood=TRUE)
                loglik <- out$objective
            }
        }
		est.pars<-out$solution
	} else {
		if(is.null(ip)){
            if(verbose==TRUE){
                cat("Initializing...", "\n")
            }
			#Sets parameter settings for random restarts by taking the parsimony score and dividing
			#by the total length of the tree
			model.set.init<-rate.cat.set.rayDISC(phy=phy,data=data.sort,model="ER",charnum=charnum)
			opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)

			taxa.missing.data.drop <- which(is.na(data.sort[,1]))
			if(length(taxa.missing.data.drop) != 0){
			  tip.labs <- names(taxa.missing.data.drop)
			  dat <- as.matrix(data.sort)
			  dat.red <- dat[-taxa.missing.data.drop,]
			  phy.red <- drop.tip(phy, taxa.missing.data.drop)
			  dat.red <- phyDat(dat.red,type="USER", levels=c("0","1"))
			  par.score <- parsimony(phy.red, dat.red, method="fitch")/2
			}else{
			  dat <- as.matrix(data.sort)
			  dat <- phyDat(dat,type="USER", levels=c("0","1"))
			  par.score <- parsimony(phy, dat, method="fitch")/2
			}
			
			tl <- sum(phy$edge.length)
			mean.change = par.score/tl
			if(mean.change==0){
				ip=0.01 + lb
			}else{
				ip<-rexp(1, 1/mean.change)
			}
			if(ip < lb || ip > ub){ # initial parameter value is outside bounds
				ip <- lb
			}
			lower.init = rep(lb, model.set.init$np)
			upper.init = rep(ub, model.set.init$np)
			init = nloptr(x0=rep(ip, length.out = model.set.init$np), eval_f=dev.raydisc, lb=lower.init, ub=upper.init, opts=opts, phy=phy,liks=model.set.init$liks,Q=model.set.init$Q,rate=model.set.init$rate,root.p=root.p, lewis.asc.bias=lewis.asc.bias)
            if(verbose == TRUE){
                cat("Finished. Beginning thorough search...", "\n")
            }
			lower = rep(lb, model.set.final$np)
			upper = rep(ub, model.set.final$np)
            if(state.recon=="subsequently") {
                out <- nloptr(x0=rep(init$solution, length.out = model.set.final$np), eval_f=dev.raydisc, lb=lower, ub=upper, opts=opts, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, lewis.asc.bias=lewis.asc.bias)
            } else {
                out <- nloptr(x0=rep(init$solution, length.out = model.set.final$np), eval_f=dev.raydisc.rates.and.states, lb=lower, ub=upper, opts=opts, phy=phy, data=data, hrm=FALSE, rate.cat=NULL, rate.mat=rate.mat, ntraits=ntraits, method=node.states, model=model, charnum=charnum, root.p=root.p, lewis.asc.bias=lewis.asc.bias, get.likelihood=TRUE)
            }
			loglik <- -out$objective
			est.pars<-out$solution
		}
		#If a user-specified starting value(s) is supplied:
		else{
            if(verbose == TRUE){
                cat("Beginning subplex optimization routine -- Starting value(s):", ip, "\n")
            }
			opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
            if(state.recon=="subsequently") {
                out <- nloptr(x0=rep(init$solution, length.out = model.set.final$np), eval_f=dev.raydisc, lb=lower, ub=upper, opts=opts, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p,lewis.asc.bias=lewis.asc.bias)
            }else{
                out <- nloptr(x0=rep(init$solution, length.out = model.set.final$np), eval_f=dev.raydisc.rates.and.states, lb=lower, ub=upper, opts=opts, phy=phy, data=data, hrm=FALSE, rate.cat=NULL, rate.mat=rate.mat, ntraits=ntraits, method=node.states, model=model, charnum=charnum, root.p=root.p, lewis.asc.bias=lewis.asc.bias, get.likelihood=TRUE)
            }
			loglik <- -out$objective
			est.pars <- out$solution
		}
	}
	#Starts the summarization process:
    if(verbose==TRUE){
        cat("Finished. Inferring ancestral states using", node.states, "reconstruction.","\n")
    }
	TIPS <- 1:nb.tip
	if(node.states == "marginal" || node.states == "scaled"){
        lik.anc <- ancRECON(phy, data, est.pars, hrm=FALSE, rate.cat=NULL, rate.mat=rate.mat, ntraits=ntraits, method=node.states, model=model, charnum=charnum, root.p=root.p)
		pr <- apply(lik.anc$lik.anc.states,1,which.max)
		phy$node.label <- pr
		tip.states <- lik.anc$lik.tip.states
		#row.names(tip.states) <- phy$tip.label
	}
    if(!state.recon == "given"){
        if(node.states == "joint"){
            lik.anc <- ancRECON(phy, data, est.pars, hrm=FALSE, rate.cat=NULL, rate.mat=rate.mat, ntraits=ntraits, method=node.states, model=model, charnum=charnum, root.p=root.p)
            phy$node.label <- lik.anc$lik.anc.states
            tip.states <- lik.anc$lik.tip.states
        }
    }else{
        if(any(is.na(phy$node.label))){
            lik.anc <- ancRECON(phy, data, est.pars, hrm=FALSE, rate.cat=NULL, rate.mat=rate.mat, ntraits=ntraits, method=node.states, model=model, charnum=charnum, root.p=root.p)
            phy$node.label <- lik.anc$lik.anc.states
            tip.states <- lik.anc$lik.tip.states
        }else{
            lik.anc <- NULL
            lik.anc$lik.anc.states <- phy$node.label
            lik.anc$lik.tip.states <- data.sort[,1]
            tip.states <- lik.anc$lik.tip.states
        }
    }
    
    if(diagn==TRUE){
        if(verbose == TRUE){
            cat("Finished. Performing diagnostic tests.", "\n")
        }
        #Approximates the Hessian using the numDeriv function
		h <- hessian(func=dev.raydisc, x=est.pars, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p, lewis.asc.bias=lewis.asc.bias)
		solution <- matrix(est.pars[model.set.final$index.matrix], dim(model.set.final$index.matrix))
		solution.se <- matrix(sqrt(diag(pseudoinverse(h)))[model.set.final$index.matrix], dim(model.set.final$index.matrix))
		hess.eig <- eigen(h,symmetric=TRUE)
		eigval<-signif(hess.eig$values,2)
		eigvect<-round(hess.eig$vectors, 2)
	}
	else{
		solution <- matrix(est.pars[model.set.final$index.matrix], dim(model.set.final$index.matrix))
		solution.se <- matrix(0,dim(solution)[1],dim(solution)[1])
		eigval<-NULL
		eigvect<-NULL
	}

	if((any(solution == lb,na.rm = TRUE) || any(solution == ub,na.rm = TRUE)) && (lb != 0 || ub != 100)){
		bound.hit <- TRUE
	}
	rownames(solution) <- rownames(solution.se) <- state.names
	colnames(solution) <- colnames(solution.se) <- state.names
	if(is.character(node.states)){
		if (node.states == "marginal" || node.states == "scaled"){
            colnames(lik.anc$lik.anc.states) <- state.names
		}
	}
	obj = list(loglik = loglik, AIC = -2*loglik+2*model.set.final$np,AICc = -2*loglik+(2*model.set.final$np*(nb.tip/(nb.tip-model.set.final$np-1))),ntraits=1, solution=solution, solution.se=solution.se, index.mat=model.set.final$index.matrix, lewis.asc.bias=lewis.asc.bias, opts=opts, data=data, phy=phy, states=lik.anc$lik.anc.states, tip.states=tip.states, iterations=out$iterations, eigval=eigval, eigvect=eigvect,bound.hit=bound.hit)
	if(!is.null(matching$message.data)){ # Some taxa were included in data matrix but not not used because they were not in the tree
		obj$message.data <- matching$message.data
		obj$data <- matching$data # Data used for analyses were different than submitted data; return this matrix
	}
	if(!is.null(matching$message.tree)){ # Some taxa were included in tree, but lacked data.  Coded as missing data.
		obj$message.tree <- matching$message.tree
		obj$data <- matching$data # Data used for analyses were different than submitted data; return this matrix
	}
	class(obj)<-"raydisc"
	return(obj)
}


#Print function
print.raydisc<-function(x,...){

	ntips=Ntip(x$phy)
	output<-data.frame(x$loglik,x$AIC,x$AICc,ntips, row.names="")
	names(output)<-c("-lnL","AIC","AICc","ntax")
	cat("\nFit\n")
	print(output)
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
	if(x$bound.hit){
		cat("At least one rate parameter equals the boundary value set by user (lb or ub).  This may be a non-optimal solution.  Try running again or changing boundary values.\n")
	}
	if(!is.null(x$message.data) || !is.null(x$message.tree)){
		cat("\nThere were differences between the tree and matrix; see message.data and/or message.tree attribute of this rayDISC object for details.\n",sep="")
	}
}


dev.raydisc.rates.and.states <- function(p, phy, data, hrm, rate.cat, rate.mat, ntraits, method, model, charnum, root.p, lewis.asc.bias, get.likelihood) {
    if(lewis.asc.bias == TRUE){
        loglik.num <- ancRECON(phy=phy, data=data, p=p, hrm=hrm, rate.cat=rate.cat, rate.mat=rate.mat, ntraits=ntraits, method=method, model=model, charnum=charnum, root.p=root.p, get.likelihood=get.likelihood)
        phy.dummy <- phy
        data.dummy <- cbind(phy$tip.label, 0)
        phy.dummy$node.label <- rep(1, Nnode(phy.dummy))
        loglik.dummy <- ancRECON(phy=phy, data=data, p=p, hrm=hrm, rate.cat=rate.cat, rate.mat=rate.mat, ntraits=ntraits, method=method, model=model, charnum=charnum, root.p=root.p, get.likelihood=get.likelihood)
        loglik <- (loglik.num  - log(1 - exp(loglik.dummy)))
    }else{
        loglik <- ancRECON(phy=phy, data=data, p=p, hrm=hrm, rate.cat=rate.cat, rate.mat=rate.mat, ntraits=ntraits, method=method, model=model, charnum=charnum, root.p=root.p, get.likelihood=get.likelihood)
    }
    return(-loglik)
}


dev.raydisc <- function(p, phy, liks, Q, rate, root.p, lewis.asc.bias){

    nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	TIPS <- 1:nb.tip
	comp <- numeric(nb.tip + nb.node)

	phy <- reorder(phy, "pruningwise")
	#Obtain an object of all the unique ancestors
	anc <- unique(phy$edge[,1])
	#This bit is to allow packages like "selac" the ability to deal with this function directly:
	if(is.null(rate)){
		Q=Q
	}else{
		if (any(is.nan(p)) || any(is.infinite(p))) return(1000000)
		Q[] <- c(p, 0)[rate]
		diag(Q) <- -rowSums(Q)
	}
    
    if(lewis.asc.bias == TRUE) {
        liks.dummy <- liks
        liks.dummy[TIPS,1] = 1
        liks.dummy[TIPS,2:dim(Q)[1]] = 0
        comp.dummy <- comp
        for (i  in seq(from = 1, length.out = nb.node)) {
            #the ancestral node at row i is called focal
            focal <- anc[i]
            #Get descendant information of focal
            desRows<-which(phy$edge[,1]==focal)
            desNodes<-phy$edge[desRows,2]
            v <- 1
            v.dummy <- 1
            for (desIndex in sequence(length(desRows))){
                v <- v*expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks[desNodes[desIndex],]
                v.dummy <- v.dummy * expm(Q * phy$edge.length[desRows[desIndex]]) %*% liks.dummy[desNodes[desIndex],]
            }
            comp[focal] <- sum(v)
            comp.dummy[focal] <- sum(v.dummy)
            liks[focal, ] <- v/comp[focal]
            liks.dummy[focal, ] <- v.dummy/comp.dummy[focal]
        }
    }else{
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
            if(lewis.asc.bias == TRUE){
                loglik.num <- (sum(log(comp[-TIPS])) + log(sum(exp(log(flat.root)+log(liks[root,])))))
                loglik.denom <- (sum(log(comp.dummy[-TIPS])) + log(sum(exp(log(flat.root)+log(liks.dummy[root,])))))
                loglik <- -(loglik.num  - log(1 - exp(loglik.denom)))
            }else{
                loglik<- -(sum(log(comp[-TIPS])) + log(sum(flat.root * liks[root,])))
            }
        }
        else{
            if(is.character(root.p)){
                # root.p==yang will fix root probabilities based on the inferred rates: q10/(q01+q10)
                if(root.p == "yang"){
                    diag(Q) = 0
                    equil.root = colSums(Q) / sum(Q)
                    if(lewis.asc.bias == TRUE){
                        loglik.num <- (sum(log(comp[-TIPS])) + log(sum(exp(log(equil.root)+log(liks[root,])))))
                        loglik.denom <- (sum(log(comp.dummy[-TIPS])) + log(sum(exp(log(equil.root)+log(liks.dummy[root,])))))
                        loglik <- -(loglik.num  - log(1 - exp(loglik.denom)))
                    }else{
                        loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(equil.root)+log(liks[root,])))))
                    }
                    if(is.infinite(loglik)){
                        return(1000000)
                    }
                }else{
                    # root.p==maddfitz will fix root probabilities according to FitzJohn et al 2009 Eq. 10:
                    root.p = liks[root,] / sum(liks[root,])
                    if(lewis.asc.bias == TRUE){
                        loglik.num <- (sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,])))))
                        loglik.denom <- (sum(log(comp.dummy[-TIPS])) + log(sum(exp(log(root.p)+log(liks.dummy[root,])))))
                        loglik <- -(loglik.num  - log(1 - exp(loglik.denom)))
                    }else{
                        loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,])))))
                    }
                }
            }
            # root.p!==NULL will fix root probabilities based on user supplied vector:
            else{
                if(lewis.asc.bias == TRUE){
                    loglik.num <- (sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,])))))
                    loglik.denom <- (sum(log(comp.dummy[-TIPS])) + log(sum(exp(log(root.p)+log(liks.dummy[root,])))))
                    loglik <- -(loglik.num  - log(1 - exp(loglik.denom)))
                }else{
                    loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,])))))
                }
                if(is.infinite(loglik)){
                    return(1000000)
                }
            }
        }
	}
	loglik
}


rate.cat.set.rayDISC<-function(phy,data,model,charnum){
	k <- 1
	factored <- factorData(data, charnum=charnum)
	nl <- ncol(factored)
	obj <- NULL
	nb.tip<-length(phy$tip.label)
	nb.node <- phy$Nnode
	#rate is a matrix of rate categories (not actual rates)
	rate<-rate.mat.maker(hrm=FALSE,ntraits=1,nstates=nl,model=model)
    index.matrix<-rate
	rate[is.na(rate)]<-max(rate,na.rm=T)+1

	stateTable <- NULL # will hold 0s and 1s for likelihoods of each state at tip
	for(column in 1:nl){
		stateTable <- cbind(stateTable,factored[,column])
	}
	colnames(stateTable) <- colnames(factored)

	ancestral <- matrix(0,nb.node,nl) # all likelihoods at ancestral nodes will be 0
	liks <- rbind(stateTable,ancestral) # combine tip likelihoods & ancestral likelihoods
	rownames(liks) <- NULL

	Q <- matrix(0, nl^k, nl^k)

	obj$np<-max(rate)-1
	obj$rate<-rate
	obj$index.matrix<-index.matrix
	obj$liks<-liks
	obj$Q<-Q

	return(obj)
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

##############
#  findAmps  #
##############
# A function to find positions of ampersands for separating different states.
# Will allow character state to be greater than one character long.
findAmps <- function(string, charnum){
  if (is.na(string)) {
    return(NULL)
  }
	if (!is.character(string)) {
	  return(NULL)
	}
	locs <- NULL # Will hold location values
	for(charnum in 1:nchar(as.character(string))){
		if(substr(string,charnum,charnum) == "&"){
			locs <- c(locs,charnum)
		}
	}
	return(locs)
}

##############
# factorData #
##############
# Function to make factored matrix as levels are discovered.
factorData <- function(data,whichchar=1,charnum){
  charcol <- whichchar+1
  factored <- NULL # will become the matrix.  Starts with no data.
  lvls <- NULL
  numrows <- length(data[,charcol])
  missing <- NULL
  
  for(row in 1:numrows){
    currlvl <- NULL
    currdata <- data[row,charcol]
    if (is.na(currdata)) {
      missing <- c(missing, row)
    } else {
      levelstring <- as.character(currdata)
      ampLocs <- findAmps(levelstring, charnum)
      if(length(ampLocs) == 0) { #No ampersands, character is monomorphic
        currlvl <- levelstring
        if(currlvl == "?" || currlvl == "-" || currlvl == "NA"){ # Check for missing data
          missing <- c(missing,row) # add to list of taxa with missing values, will fill in entire row later
        } else { # Not missing data
          if(length(which(lvls == currlvl)) == 0){# encountered a level not seen yet
            if(length(factored) == 0){ # Matrix is empty, need to create it
              factored <- matrix(0,numrows,1)
              colnames(factored) <- currlvl
              rownames(factored) <- rownames(data)
            } else { # matrix already exists, but need to add a column for the new level
              zerocolumn <- rep(0,numrows)
              factored <- cbind(factored, zerocolumn)
              colnames(factored)[length(factored[1,])] <- currlvl
            }
            lvls <- c(lvls,currlvl) # add that level to the list
          } # already found this level in another state.  Set the value to one
          whichlvl <- which(lvls == currlvl) # this index number should correspond to the column number of the state
          factored[row,whichlvl] <- 1
        }
      } else { #At least one ampersand found, polymorphic character
        start <- 1
        numlvls <- length(ampLocs)+1
        for(part in 1:numlvls){
          # Pull out level from levelstring
          if(part <= length(ampLocs)){ # Havent reached the last state
            currlvl <- substr(levelstring,start,(ampLocs[part]-1)) # pull out value between start and the location-1 of the next ampersand
          } else { # Final state in list
            currlvl <- substr(levelstring,start,nchar(levelstring)) # pull out value between start and the last character of the string
          }
          if(currlvl == "?" || currlvl == "-"){ # Missing data, but polymorphic?
            missing <- c(missing,row) # add to list of taxa with missing values, will fill in entire row later
          }
          else { # Not missing data
            if(length(which(lvls == currlvl)) == 0){# encountered a level not seen yet
              if(length(factored) == 0){ # Matrix is empty, need to create it
                factored <- matrix(0,numrows,1)
                colnames(factored) <- currlvl
                rownames(factored) <- rownames(data)
              } else { # matrix already exists, but need to add a column for the new level
                zerocolumn <- rep(0,numrows)
                factored <- cbind(factored, zerocolumn)
                colnames(factored)[length(factored[1,])] <- currlvl
              }
              lvls <- c(lvls,currlvl) # add that level to the list
            } # already found this level in another state.  Set the value to one
            whichlvl <- which(lvls == currlvl) # this index number should correspond to the column number of the state
            factored[row,whichlvl] <- 1
            start <- ampLocs[part] + 1
          }
        }
      }
      
    }
  }
  #Need to deal with any rows with missing data; fill in NA for all columns for that row
  for(missingrows in 1:length(missing)){
    for(column in 1:length(factored[1,])){
      factored[missing[missingrows],column] <- 1 # All states equally likely
    }
  }
  factored <- factored[,order(colnames(factored))]
  return(factored)
}