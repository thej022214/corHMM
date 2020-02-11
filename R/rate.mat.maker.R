#Rate matrix maker and manipulating functions

#written by Jeremy M. Beaulieu and Jeffrey C. Oliver, modified by James D. Boyko
rate.mat.maker<- function (rate.cat, hrm = TRUE, ntraits = NULL, nstates = NULL, 
          model = c("ER", "SYM", "ARD")) 
{
  if (hrm == TRUE) {
    k = 2
    mat1 <- matrix(NA, k * rate.cat, k * rate.cat)
    mat2 <- matrix(NA, k * rate.cat, k * rate.cat)
    vec.tmp1 <- rep(c(0, 1), rate.cat)
    vec.tmp2 <- rep(1:rate.cat, rep(2, rate.cat)) - 1
    for (i in 1:(k * rate.cat)) {
      mat1[i, ] <- abs(vec.tmp1 - vec.tmp1[i])
      mat2[i, ] <- abs(vec.tmp2 - vec.tmp2[i])
    }
    matFINAL <- mat1 + mat2
    rate.mat.index <- matrix(NA, k * rate.cat, k * rate.cat)
    np <- k + (rate.cat - 1) * 6
    index <- matFINAL == 1
    rate.mat.index[index] <- 1:np
    if (rate.cat == 1) {
      rownames(rate.mat.index) <- c("(0)", "(1)")
      colnames(rate.mat.index) <- c("(0)", "(1)")
    }
    if (rate.cat == 2) {
      rownames(rate.mat.index) <- c("(0,R1)", "(1,R1)", 
                                    "(0,R2)", "(1,R2)")
      colnames(rate.mat.index) <- c("(0,R1)", "(1,R1)", 
                                    "(0,R2)", "(1,R2)")
    }
    if (rate.cat == 3) {
      rownames(rate.mat.index) <- c("(0,R1)", "(1,R1)", 
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)")
      colnames(rate.mat.index) <- c("(0,R1)", "(1,R1)", 
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)")
    }
    if (rate.cat == 4) {
      rownames(rate.mat.index) <- c("(0,R1)", "(1,R1)", 
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)", 
                                    "(1,R4)")
      colnames(rate.mat.index) <- c("(0,R1)", "(1,R1)", 
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)", 
                                    "(1,R4)")
    }
    if (rate.cat == 5) {
      rownames(rate.mat.index) <- c("(0,R1)", "(1,R1)", 
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)", 
                                    "(1,R4)", "(0,R5)", "(1,R5)")
      colnames(rate.mat.index) <- c("(0,R1)", "(1,R1)", 
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)", 
                                    "(1,R4)", "(0,R5)", "(1,R5)")
    }
    if (rate.cat == 6) {
      rownames(rate.mat.index) <- c("(0,R1)", "(1,R1)", 
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)", 
                                    "(1,R4)", "(0,R5)", "(1,R5)", "(0,R6)", "(1,R6)")
      colnames(rate.mat.index) <- c("(0,R1)", "(1,R1)", 
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)", 
                                    "(1,R4)", "(0,R5)", "(1,R5)", "(0,R6)", "(1,R6)")
    }
    if (rate.cat == 7) {
      rownames(rate.mat.index) <- c("(0,R1)", "(1,R1)", 
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)", 
                                    "(1,R4)", "(0,R5)", "(1,R5)", "(0,R6)", "(1,R6)", 
                                    "(0,R7)", "(1,R7)")
      colnames(rate.mat.index) <- c("(0,R1)", "(1,R1)", 
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)", 
                                    "(1,R4)", "(0,R5)", "(1,R5)", "(0,R6)", "(1,R6)", 
                                    "(0,R7)", "(1,R7)")
    }
  }
  if (hrm == FALSE) {
    k = ntraits
    nl = 2
    if (ntraits == 1) {
      k <- 1
      nl <- nstates
      if (is.character(model)) {
        rate.mat.index <- matrix(NA, nl, nl)
        tmp2 <- cbind(1:(nl^k), 1:(nl^k))
        index <- matrix(TRUE, nl^k, nl^k)
        diag(index) <- FALSE
        if (model == "ER") {
          np <- 1
          rate.mat.index[index] <- 1:np
        }
        if (model == "SYM") {
          np <- nl * (nl - 1)/2
          sel <- col(rate.mat.index) < row(rate.mat.index)
          rate.mat.index <- t(rate.mat.index)
          rate.mat.index[sel] <- 1:np
          rate.mat.index[upper.tri(rate.mat.index)] = t(rate.mat.index)[upper.tri(rate.mat.index)]
        }
        if (model == "ARD") {
          np <- nl * (nl - 1)
          rate.mat.index[index] <- 1:np
        }
      }
    }
    if (ntraits == 2) {
      mat1 <- matrix(, nl^k, nl^k)
      mat2 <- matrix(, nl^k, nl^k)
      vec.tmp1 <- c(0, 0, 1, 1)
      vec.tmp2 <- c(0, 1, 0, 1)
      for (i in 1:(nl^k)) {
        mat1[i, ] <- abs(vec.tmp1 - vec.tmp1[i])
        mat2[i, ] <- abs(vec.tmp2 - vec.tmp2[i])
      }
      matFINAL <- mat1 + mat2
      if (is.character(model)) {
        rate.mat.index <- matrix(NA, nl^k, nl^k)
        if (model == "ER") {
          np <- 1
          index <- matFINAL == 1
          rate.mat.index[index] <- 1:np
        }
        if (model == "SYM") {
          np <- 4
          index <- matFINAL == 1
          rate.mat.index[index][c(1, 2, 4, 6)] <- rate.mat.index[index][c(3, 
                                                                          5, 7, 8)] <- 1:np
        }
        if (model == "ARD") {
          np <- 8
          index <- matFINAL == 1
          rate.mat.index[index] <- 1:np
        }
      }
    }
    if (ntraits == 3) {
      mat1 <- matrix(, nl^k, nl^k)
      mat2 <- matrix(, nl^k, nl^k)
      mat3 <- matrix(, nl^k, nl^k)
      vec.tmp1 <- c(0, 1, 0, 0, 1, 1, 0, 1)
      vec.tmp2 <- c(0, 0, 1, 0, 1, 0, 1, 1)
      vec.tmp3 <- c(0, 0, 0, 1, 0, 1, 1, 1)
      for (i in 1:(nl^k)) {
        mat1[i, ] <- abs(vec.tmp1 - vec.tmp1[i])
        mat2[i, ] <- abs(vec.tmp2 - vec.tmp2[i])
        mat3[i, ] <- abs(vec.tmp3 - vec.tmp3[i])
      }
      matFINAL <- mat1 + mat2 + mat3
      if (is.character(model)) {
        rate.mat.index <- matrix(NA, nl^k, nl^k)
        if (model == "ER") {
          np <- 1
          index <- matFINAL == 1
          rate.mat.index[index] <- 1:np
        }
        if (model == "SYM") {
          np <- 12
          index <- matFINAL == 1
          rate.mat.index[index][c(1, 2, 3, 5, 6, 8, 9, 
                                  11, 12, 15, 18, 21)] <- rate.mat.index[index][c(4, 
                                                                                  7, 10, 13, 16, 14, 19, 17, 20, 22, 23, 24)] <- 1:np
        }
        if (model == "ARD") {
          np <- 24
          index <- matFINAL == 1
          rate.mat.index[index] <- 1:np
        }
      }
    }
    if (ntraits == 1) {
      rownames(rate.mat.index) <- as.character(1:nl)
      colnames(rate.mat.index) <- as.character(1:nl)
    }
    if (ntraits == 2) {
      rownames(rate.mat.index) <- c("(0,0)", "(0,1)", "(1,0)", 
                                    "(1,1)")
      colnames(rate.mat.index) <- c("(0,0)", "(0,1)", "(1,0)", 
                                    "(1,1)")
    }
    if (ntraits == 3) {
      rownames(rate.mat.index) <- c("(0,0,0)", "(1,0,0)", 
                                    "(0,1,0)", "(0,0,1)", "(1,1,0)", "(1,0,1)", "(0,1,1)", 
                                    "(1,1,1)")
      colnames(rate.mat.index) <- c("(0,0,0)", "(1,0,0)", 
                                    "(0,1,0)", "(0,0,1)", "(1,1,0)", "(1,0,1)", "(0,1,1)", 
                                    "(1,1,1)")
    }
  }
  return(rate.mat.index)
}

rate.par.drop <- function(rate.mat.index=NULL,drop.par=NULL){
	if(is.null(rate.mat.index)){
		stop("Rate matrix needed.  See mat.maker to create one.\n")
	}
	if(is.null(drop.par)){
		cat("No parameters indicated to drop.  Original matrix returned.\n")
		return(rate.mat.index)
	}
	if(max(rate.mat.index,na.rm=TRUE) < max(drop.par,na.rm=TRUE)){
		cat("Some parameters selected for dropping were not in the original matrix.\n")
	}
	drop.par <- unique(drop.par) # in case parameters listed more than once in drop vector
	drop.par <- drop.par[order(drop.par)]
	max <- max(rate.mat.index,na.rm=TRUE)
	for(drop.which in 1:length(drop.par)){
		drop.locs <- which(rate.mat.index == drop.par[drop.which],arr.ind=TRUE)
		rate.mat.index[drop.locs] <- NA
	}
	max <- max - length(drop.par)
	exclude <- which(is.na(rate.mat.index))
	rate.mat.index[-exclude] <- 1:max
	
	return(rate.mat.index)
}

rate.par.eq <- function(rate.mat.index=NULL,eq.par=NULL){
	if(is.null(rate.mat.index)){
		stop("Rate matrix needed.  See mat.maker to create one.\n")
	}
	if(is.null(drop) || length(eq.par) < 2){
		cat("Fewer than two parameters indicated to equalize.  Original matrix returned.\n")
		return(rate.mat.index)
	}
	too.big <- which(eq.par > max(rate.mat.index,na.rm=TRUE))
	if(length(too.big) > 0){
		cat("Some parameters selected for equalizing were not in the original matrix:\n")
		cat("Not in original rate.mat.index:",eq.par[too.big],"\n")
		cat("Original matrix returned.\n")
		return(rate.mat.index)
	}
	eq.par <- unique(eq.par)
	eq.par <- eq.par[order(eq.par)]	
	min <- min(eq.par) # rm.na unnecessary?

	# the decrement index will hold counters to decrement rate index
	dec.index <- matrix(0,length(rate.mat.index[,1]),length(rate.mat.index[1,]))
	for(eq.which in 2:length(eq.par)){
		to.eq <- which(rate.mat.index == eq.par[eq.which],arr.ind=TRUE)
		rate.mat.index[to.eq] <- min
	}
	# the decrement index will hold counters to decrement rate index
	dec.index <- matrix(0,length(rate.mat.index[,1]),length(rate.mat.index[1,]))
	for(eq.which in 2:length(eq.par)){
		to.dec <- which(rate.mat.index > eq.par[eq.which],arr.ind=TRUE) #greater than current decrementer
		dec.index[to.dec] <- dec.index[to.dec] + 1
	}
	rate.mat.index <- rate.mat.index - dec.index
		
	return(rate.mat.index)
}

# JDB modifications - alternative index matrix construction. instructions below.

convert2I <- function(Q){
  I <- matrix(0, dim(Q)[1], dim(Q)[2])
  diag(I) <- 1
  return(I)
}

getStateMat <- function(nState){
  StateMat <- matrix(1:(nState^2),nState,nState)
  diag(StateMat) <- 0
  StateMat[StateMat>0] <- 1:length(StateMat[StateMat>0])
  return(StateMat)
}

dropStateMatPars <- function(StateMat, Pars){
  for(i in Pars){
    StateMat[StateMat==i] <- 0
  }
  pars <- unique(as.vector(StateMat))[-1]
  for(i in 1:length(pars)){
    StateMat[StateMat == pars[i]] <- i
  }
  return(StateMat)
}

keepStateMatPars <- function(StateMat, Pars){
  tmp <- matrix(0, dim(StateMat)[1], dim(StateMat)[2])
  for(i in Pars){
    tmp[StateMat==i] <- 1
  }
  tmp[tmp>0] <- 1:length(tmp[tmp>0])
  return(tmp)
}

equateStateMatPars <- function(StateMat, ParsList){
  if(!class(ParsList) == "list"){
    ParsList <- list(ParsList)
  }
  TestMat <- do.call(rbind, ParsList)
  RowLow <- apply(TestMat, 1, min)
  ParsList <- ParsList[match(sort(RowLow), RowLow)]
  newMat <- StateMat
  pars <- StateMat[StateMat > 0]
  for(i in 1:length(ParsList)){
    newMat[StateMat %in% ParsList[[i]]] <- min(ParsList[[i]])
  }
  pars <- newMat[newMat > 0]
  for(i in 1:length(unique(pars))){
    if(!i %in% pars){
      pars[pars == min(pars[pars > i])] <- i
    }
  }
  newMat[newMat > 0] <- pars
  return(newMat)
}

updateStateMats <- function(StateMats){
  tmp <- unlist(lapply(StateMats, max)) 
  newIndMax <- sapply(1:length(tmp), function(x) sum(tmp[1:x]))
  for(i in 2:length(StateMats)){
    StateMats[[i]][StateMats[[i]] > 0] <- ((newIndMax[i-1]+1):newIndMax[i])[(StateMats[[i]])]
  }
  return(StateMats)
}

getFullMat <- function(StateMats, RateClassMat){
  FullMat <- convert2I(RateClassMat) %x% matrix(0, dim(StateMats[[1]])[1], dim(StateMats[[1]])[2])
  IndMat <- matrix(0, length(StateMats), length(StateMats))
  for(i in 1:length(StateMats)){
    tmpIndMat <- IndMat
    diag(tmpIndMat)[i] <-  1
    FullMat <- FullMat + tmpIndMat %x% StateMats[[i]]
  }
  StartingHiddenInd <- max(FullMat)
  NonDiagonal <- c(which(lower.tri(RateClassMat) & RateClassMat > 0), 
                   which(upper.tri(RateClassMat) & RateClassMat > 0))
  # a matrix for determining where in the full mat we're going
  IndMat <- matrix(0, dim(RateClassMat)[1], dim(RateClassMat)[2])
  # a matrix that contains the diagonals being implemented
  HRMat <- matrix(0, dim(StateMats[[1]])[1], dim(StateMats[[1]])[2])
  for(i in 1:length(NonDiagonal)){
    tmpIndMat <- IndMat
    tmpHRMat <- HRMat
    diag(tmpHRMat) <- StartingHiddenInd + RateClassMat[NonDiagonal[i]]
    tmpIndMat[NonDiagonal[i]] <-  1
    FullMat <- FullMat + tmpIndMat %x% tmpHRMat
  }
  StateNames <- paste("(", rep(1:(dim(StateMats[[1]])[1]), dim(RateClassMat)[1]), ",", rep(paste("R", 1:dim(RateClassMat)[1], sep = ""), each = dim(StateMats[[1]])[1]), ")", sep = "")
  colnames(FullMat) <- rownames(FullMat) <- StateNames
  return(FullMat)
}

rate.mat.maker.JDB <-function(rate.cat, hrm=TRUE, ntraits=2, nstates=NULL, model="ARD"){
  
  if(rate.cat == 1){
    FullMat <- getStateMat(ntraits)
    StateNames <- paste("(", rep(1:ntraits, rate.cat), ",", rep(paste("R", 1:rate.cat, sep = ""), each = ntraits), ")", sep = "")
    rownames(FullMat) <- colnames(FullMat) <- StateNames
    if(model == "ER"){
      FullMat[FullMat > 0] <- 1
    }
    if(model == "SYM"){
      FullMat[upper.tri(FullMat)] <- 1:length(FullMat[upper.tri(FullMat)])
      FullMat <- t(FullMat) 
      FullMat[upper.tri(FullMat)] <- 1:length(FullMat[upper.tri(FullMat)])
    }
    FullMat[FullMat == 0] <- NA
    return(FullMat)
  }
  
  StateMats <- vector("list", rate.cat)
  
  #i should have put the if statements outside... but it's w/e
  for(i in 1:rate.cat){
    if(model == "ARD"){
      StateMats[[i]] <- getStateMat(ntraits)
    }
    
    if(model == "ER"){
      FullMat <- getStateMat(ntraits)
      FullMat[FullMat > 0] <- 1
      StateMats[[i]] <- FullMat
    }
    
    if(model == "SYM"){
      FullMat <- getStateMat(ntraits)
      FullMat[upper.tri(FullMat)] <- 1:length(FullMat[upper.tri(FullMat)])
      FullMat <- t(FullMat) 
      FullMat[upper.tri(FullMat)] <- 1:length(FullMat[upper.tri(FullMat)])
      StateMats[[i]] <- FullMat
    }
  }
  
  RateClassMat <- getStateMat(rate.cat)
  StateMats <- updateStateMats(StateMats) 
  FullMat <- getFullMat(StateMats, RateClassMat)
  
  StateNames <- paste("(", rep(1:ntraits, rate.cat), ",", rep(paste("R", 1:rate.cat, sep = ""), each = ntraits), ")", sep = "")
  
  rownames(FullMat) <- colnames(FullMat) <- StateNames
  FullMat[FullMat == 0] <- NA
  return(FullMat)
}

getRateMat4Dat <- function(data, phy, rate.cat = 1, model = "ARD"){
  
  data.legend <- input.data <- data
  nCol <- dim(data)[2]
  # convert data to numeric
  for(i in 2:nCol){
    data[,i] <- as.factor(data[,i])
  }
  
  # will automatically detect if the input data has multiple columns and convert it to corHMM format.
  if(nCol > 2){
    old.data <- apply(data[,2:nCol], 1, function(x) paste(c(x), collapse = "_"))
    Traits <- unique(old.data)
    nTraits <- length(Traits)
    data <- data.frame(sp = data[,1], d = match(old.data, Traits))
    names(Traits) <- 1:nTraits
  }
  
  rate.mat <- rate.mat.maker.JDB(rate.cat = rate.cat, ntraits = nTraits, model = model)
  rate.mat[is.na(rate.mat)] <- 0
  legend <- gsub("_", " & ", Traits)
  res <- list(legend = legend, rate.mat = rate.mat)
  return(res)
}


# step 1: create the individual processes that you are trying to model (the number of matricies you create is the number of hidden states you're interested in modeling)
# StateMatA <- getStateMat(nState = 3)
# StateMatB <- getStateMat(nState = 3)
# StateMatC <- getStateMat(nState = 3)
# StateMatD <- getStateMat(nState = 3)
# 
# StateMats <- list(StateMatA, StateMatB, StateMatC, StateMatD)
# 
# # step 2: determine how those processes are related to one another
# RateClassMat <- getStateMat(4)
# 
# # this updates all of the matricies to have estimated paramaters
# 
# #StateMats <- updateStateMats(StateMats) 
# 
# # step 3: combine the processes and their relations into a single matrix
# FullMat <- getFullMat(StateMats, RateClassMat)
# 
# require(hisse)
# trans.rate <- TransMatMakerMuHiSSE(hidden.traits=3, include.diagonals = TRUE)
# 
# # hardcoded beacuse lazy 3 state
# hiMat <- rbind(cbind(FullMat, 0, 0, 0, 0),0,0,0,0)
# nObsStates <- dim(StateMats[[1]])[1]
# nHiStates <- length(StateMats)
# nObsParams <- nHiStates*nObsStates
# IndMat <- matrix(1:nObsParams, nObsStates, nHiStates)
# IndVec <- as.vector(rbind(IndMat, (nObsParams+1):(nObsParams+dim(IndMat)[2])))
# 
# hiMat <- hiMat[IndVec,]
# hiMat <- hiMat[,IndVec]
# 
# trans.rate[1:dim(trans.rate)[1], 1:dim(trans.rate)[2]] <- hiMat


## Common function used by multiple methods.

#Note that James Boyko came up with the FindGenerations code, so the insanity of this Rumsfeldian speak is all him....
FindGenerations <- function(phy){
  generation <- list()
  known <- 1:Ntip(phy)
  unknown <- phy$edge[,1]
  needed <- phy$edge[,2]
  root <- min(unknown)
  i <- 1
  repeat{
    knowable <- unknown[needed %in% known]
    knowable <- knowable[duplicated(knowable)]
    generation[[i]] <-  knowable
    
    known <- c(known, knowable)
    needed <- needed[!unknown %in% knowable]
    unknown <- unknown[!unknown %in% knowable]
    i <- i + 1
    if (any(root == knowable)) break
  }
  res <- generation
  return(res)
}
