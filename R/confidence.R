ComputeCI <- function(corhmm.object, desired.delta = 2, n.points=5000, verbose=TRUE,  print_freq=50, ...) {
	best.neglnl <- -corhmm.object$loglik
	index.mat <- corhmm.object$index.mat
	raw.rates <- corhmm.object$solution
	corhmm.object$node.states <- "none" #we don't want to waste time on this

	par <- MatrixToPars(corhmm.object)
	par.best <- par

	## compute_neglnlikelihood <- function(...) {
	## 	return(-compute_lnlikelihood(...))
        ## }
        if (!requireNamespace("dentist")) {
          stop("please install the 'dentist' package via remotes::install_github('bomeara/dentist'")
        }
	dented_results <- dentist::dent_walk(par=par.best, fn=compute_neglnlikelihood, best_neglnL=best.neglnl,  nsteps=n.points, print_freq=print_freq, corhmm.object=corhmm.object)

	#class(dented_results) <- append("corhmm_confidence", class(dented_results))
    return(dented_results)
}

#' @export
print.corhmm_confidence <- function(x, ...) {
	obj_rates <- as.data.frame(x[,-1])
	best_index <- which.max(x$lnL)[1]
	cat("Range of log likelihoods is ", max(x$lnL), " to ", min(x$lnL), " a difference of ", max(x$lnL)- min(x$lnL), "\n", sep="")
	result <- data.frame(min=apply(obj_rates, 2, min), best=simplify2array(unname(obj_rates[best_index,,drop=TRUE])), max=apply(obj_rates, 2, max))
	rownames(result) <- colnames(obj_rates)
	for (i in sequence(ncol(obj_rates))) {
		for (j in sequence(ncol(obj_rates))) {
			if(i!=j) {
				ratios <- obj_rates[,i] / obj_rates[,j]
				result <- rbind(result, c(min(ratios), obj_rates[best_index,i] / obj_rates[best_index,j], max(ratios)))
				rownames(result)[nrow(result)] <- paste0(colnames(obj_rates)[i], " / ", colnames(obj_rates)[j])
			}
		}
	}
	print(result)
}


HessianConfidence <- function(corhmm.object) {
	index.mat <- corhmm.object$index.mat
	raw.rates <- corhmm.object$solution
	corhmm.object$node.states <- "none" #we don't want to waste time on this

	par <- MatrixToPars(corhmm.object)
	hess <- numDeriv::hessian(compute_neglnlikelihood, par, "Richardson", method.args=list(), corhmm.object)
	fisher_info<-solve(-hess)
	prop_sigma<-sqrt(diag(fisher_info))
	return.matrix <- rbind(par-1.96*prop_sigma, par+1.96*prop_sigma)
	return.matrix[which(return.matrix<0)] <- 0
	rownames(return.matrix) <- c("lower", "upper")
	colnames(return.matrix) <- names(par)
	return(return.matrix)
}


# simple function; returns log likelihood
compute_neglnlikelihood <- function(par, corhmm.object) {

  # legacy code. all previous verions of corhmm don't give the collapse as an output, but it's now needed to create the matrices correctly. 
  if(is.null(corhmm.object$collapse)){
    collapse = TRUE # default setting
  }else{
    collapse = corhmm.object$collapse
  }
	corhmm.object$order.test <- FALSE
	corhmm.object$phy$node.label <- NULL
	nObs <- dim(corhmm.object$index.mat)[1]/corhmm.object$rate.cat
	model.set.final <- rate.cat.set.corHMM.JDB(phy = corhmm.object$phy, data = corhmm.object$data, rate.cat = corhmm.object$rate.cat, ntraits = nObs, model = "ARD", collapse=collapse)
	phy <- reorder(corhmm.object$phy, "pruningwise")
	rate.mat <- corhmm.object$index.mat
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
	result <- dev.corhmm(
		p = log(par),
		phy = corhmm.object$phy,
		liks = model.set.final$liks,
		Q = model.set.final$Q,
		rate = model.set.final$rate,
		root.p = corhmm.object$root.p,
		rate.cat = corhmm.object$rate.cat,
		order.test = corhmm.object$order.test,
		lewis.asc.bias = ifelse(any(grepl("lewis.asc.bias", names(corhmm.object))), corhmm.object$lewis.asc.bias, FALSE),
	  set.fog = FALSE, 
	  fog.vec = model.set.final$fog.vec
	)

	#return(dev.corhmm(log(par), corhmm_object$phy, liks=42, Q=42, rate=42, root.p=corhmm.object$root.p, rate.cat=42, order.test=42, lewis.asc.bias=ifelse(any(grepl("lewis.asc.bias", names(corhmm.object))), corhmm.object$lewis.asc.bias, FALSE)))
	return(result)
}


GenerateValues <- function(par, lower, upper, max.tries=100, expand.prob=0, examined.max, examined.min) {
    if(is.null(lower)) {
        lower <- 0.1*par
    }
    if(is.null(upper)) {
        upper <- 10*par
    }
    pass=FALSE
    tries=0
    while(!pass && tries<=max.tries) {
        tries <- tries+1
        pass=TRUE
        new.vals <- rep(NA, length(par))
        for(i in sequence(length(par))) {
            examined.max[i]<-max(0.001, examined.max[i])
            new.vals.bounds <- sort(c(max(lower[i], 0.9*examined.min[i]), min(upper[i], 1.1*examined.max[i])), decreasing=FALSE)
            new.vals[i]<-stats::runif(1, min=ifelse(is.finite(new.vals.bounds[1]),new.vals.bounds[1], 0.000001) , max=ifelse(is.finite(new.vals.bounds[2]), new.vals.bounds[2], 10000))

            if(new.vals[i]<lower[i]) {
                pass=FALSE
            }
            if(new.vals[i]>upper[i]) {
                pass=FALSE
            }
        }
    }
    if(tries>max.tries) {
        return(NA)
    }
    names(new.vals) <- names(par)
    return(new.vals)
}


GetGridPoints <- function(lower, upper, n.points) {
	n.grains <- ceiling(n.points^(1/length(lower)))
	vector.list <- list()
	for (i in seq_along(lower)) {
		vector.list[[i]] <- seq(from=lower[i], to=upper[i], length.out=n.grains)
	}
	parameter.matrix <- expand.grid(vector.list)
	return(parameter.matrix)
}


GetLHSPoints <- function(lower, upper, n.points) {
	raw.points <- randomLHS(n=n.points, k=length(lower))
	parameter.matrix <- NA*raw.points
	for (i in seq_along(lower)) {
		parameter.matrix[,i] <- qunif(raw.points[,i], min=lower[i], max=upper[i])
	}
	colnames(parameter.matrix) <- names(lower)
	return(parameter.matrix)
}

# A which for multidimensional arrays.
# Mark van der Loo 16.09.2011
#
# A Array of booleans
# returns a sum(A) x length(dim(A)) array of multi-indices where A == TRUE
#
multi.which <- function(A){
    if ( is.vector(A) ) return(which(A))
    d <- dim(A)
    T <- which(A) - 1
    nd <- length(d)
    t( sapply(T, function(t){
        I <- integer(nd)
        I[1] <- t %% d[1]
        sapply(2:nd, function(j){
            I[j] <<- (t %/% prod(d[1:(j-1)])) %% d[j]
        })
        I
    }) + 1 )
}

ParsToMatrix <- function(pars, corhmm.object) {
	return_mat <- matrix(pars[corhmm.object$index.mat], dim(corhmm.object$index.mat))
	rownames(return_mat) <- rownames(corhmm.object$solution)
	colnames(return_mat) <- colnames(corhmm.object$solution)
	return(return_mat)
}

MatrixToPars <- function(corhmm.object) {
	index.mat <- corhmm.object$index.mat
	raw.rates <- corhmm.object$solution

	par <- rep(NA,max(index.mat, na.rm=TRUE))
	for (i in seq_along(par)) {
		par[i] <- raw.rates[which(index.mat==i)][1]
		relevant_indices <- multi.which(index.mat==i)[1,]
		names(par)[i] <- paste0(rownames(raw.rates)[relevant_indices[1]]," -> ", colnames(raw.rates)[relevant_indices[2]])
	}
	return(par)
}


getParamVector <- function(corhmm.object) {
  par <- corhmm.object$solution[!is.na(corhmm.object$solution)]
  position <- 0
  for (col_index in sequence(nrow(corhmm.object$solution))) {
    for(row_index in sequence(nrow(corhmm.object$solution))) {
      if(!is.na(corhmm.object$solution[row_index, col_index])) {
        position <- position+1
        par[position] <- corhmm.object$solution[row_index, col_index]
        names(par)[position] <- paste0(rownames(corhmm.object$solution)[row_index], " to ", colnames(corhmm.object$solution)[col_index])
      }
    }
  }
  return(par)
}
