

######################################################################################################################################
######################################################################################################################################
### Unit Test Code -- Does Brute force marginal which is time consuming, but provides a check to ensure Maddison 1995 is correct
######################################################################################################################################
######################################################################################################################################

##Does brute force reconstruction -- at node i, fix in state i, then fix in state j, summing over all other states at all other nodes
BruteConditionalProbability <- function(phy, data, p=c(0.005,0.005), root.p=NULL, focal.to.fix=NULL, focal.state=NULL, get.conditional=NULL){
    
    rate <- corHMM:::rate.mat.maker(hrm=FALSE, ntraits=1, nstates=2, model="ARD")
    index.matrix <- rate
    rate[is.na(rate)] <- max(rate,na.rm=T) + 1
    p.new <- p
    Q <- matrix(0, 2, 2)
    
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    data.sort <- data.frame(data[,2],data[,2],row.names=data[,1])
    data.sort <- data.sort[phy$tip.label,]
    liks <- matrix(0, nb.tip + nb.node, 2)
    TIPS <- 1:nb.tip
    x <- data.sort[,1]
    for(i in 1:nb.tip){
        if(x[i]==0){liks[i,1]=1}
        if(x[i]==1){liks[i,2]=1}
        if(x[i]==2){liks[i,1:2]=1}
    }
    
    comp <- numeric(nb.tip + nb.node)
    
    #Obtain an object of all the unique ancestors
    anc <- unique(phy$edge[,1])
    #This bit is to allow packages like "selac" the ability to deal with this function directly:
    Q[] <- c(p.new, 0)[rate]
    diag(Q) <- -rowSums(Q)
    
    if(any(focal.to.fix < nb.tip)){
        which.is.tip <- which(focal.to.fix < nb.tip)
        if(focal.state[which.is.tip] != which.max(liks[focal.to.fix[which.is.tip],])){
            return(-10000000)
        }
    }
    
    for (i  in seq(from = 1, length.out = nb.node)) {
        #the ancestral node at row i is called focal
        focal <- anc[i]
        #Get descendant information of focal
        desRows <- which(phy$edge[,1]==focal)
        desNodes <- phy$edge[desRows,2]
        v <- 1
        for (desIndex in sequence(length(desRows))){
            v <- v * expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks[desNodes[desIndex],]
        }
        
        if(any(focal.to.fix == focal)){
            which.focal <- which(focal.to.fix==focal)
            fixer.focal <- numeric(2)
            fixer.focal[focal.state[which.focal]] <- 1
            v <- v * fixer.focal
        }
        
        comp[focal] <- sum(v)
        liks[focal, ] <- v/comp[focal]
    }
    root <- nb.tip + 1L
    if(!is.null(get.conditional)){
        return(liks[get.conditional,])
    }else{
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
            if(is.infinite(loglik)){
                return(-1000000)
            }else{
                return(loglik)
            }
        }else{
            if(is.character(root.p)){
                # root.p==yang will fix root probabilities based on the inferred rates: q10/(q01+q10)
                if(root.p == "yang"){
                    root.p <- Null(Q)
                    root.p <- c(root.p/sum(root.p))
                    loglik <- sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,]))))
                    if(is.infinite(loglik)){
                        return(-1000000)
                    }else{
                        return(loglik)
                    }
                }else{
                    # root.p==maddfitz will fix root probabilities according to FitzJohn et al 2009 Eq. 10:
                    root.p = liks[root,] / sum(liks[root,])
                    loglik <- sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,]))))
                    if(is.infinite(loglik)){
                        return(-1000000)
                    }else{
                        return(loglik)
                    }
                }
            }else{
                loglik <- sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,]))))
                if(is.infinite(loglik)){
                    return(-1000000)
                }else{
                    return(loglik)
                }
            }
        }
    }
}


GetMarginalBrute <- function(phy, data, p, root.p, n.states, node.fixed, state.fixed){
    phy <- reorder(phy, "pruningwise")
    anc <- unique(phy$edge[,1])
    marginals.x <- matrix(0, length(anc), n.states)
    for(node.index in 1:length(anc)){
        res <- c()
        for(state.index in 1:n.states){
            if(!is.null(node.fixed)){
                res <- c(res, BruteConditionalProbability(phy, data, p=p, root.p=root.p, focal.to.fix=c(node.fixed, anc[node.index]), focal.state=c(state.fixed,state.index), get.conditional=NULL))
            }else{
                res <- c(res, BruteConditionalProbability(phy, data, p=p, root.p=root.p, focal.to.fix=anc[node.index], focal.state=state.index, get.conditional=NULL))
            }
        }
        #res <- res + log(GetRootP(root.p=root.p, p=p, phy=phy, data=data))
        marginals.x[anc[node.index]-Ntip(phy),] <- exp(res - max(res))/sum(exp(res - max(res)))
    }
    if(!is.null(node.fixed)){
        fixer <- numeric(n.states)
        fixer[state.fixed] <- 1
        marginals.x[node.fixed - Ntip(phy),] <- fixer
    }
    return(marginals.x)
}


