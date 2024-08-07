
######################################################################################################################################
######################################################################################################################################
### Simple tests of the Marginal reconstruction algorithm
######################################################################################################################################
######################################################################################################################################

test_that("Equal rates, no nodes fixed",{
    skip_on_cran()

    require(expm)
    require(MASS)

    data(primates)
    phy <- primates[[1]]
    phy <- multi2di(phy)
    data <- primates[[2]]
    label.vector <- rep(NA, Ntip(phy) + Nnode(phy))
    homo_gorilla <- getMRCA(phy,tip=c("Homo_sapiens", "Gorilla_gorilla"))
    label.vector[homo_gorilla] <- 2
    phy$node.label <- label.vector[-c(1:Ntip(phy))]

    corHMM.new <- corHMM(phy, data[,c(1,2)], model="ER", rate.cat=1, fixed.nodes=FALSE, p=0.01080903, tip.fog=0)
    corHMM.brute <- corHMM:::GetMarginalBrute(phy=phy, data=data, p=c(0.01080903,0.01080903), root.p="yang", n.states=2, node.fixed=NULL, state.fixed=2)
    comparison <- identical(sum(round(corHMM.brute - corHMM.new$states, 5)), 0)
    expect_true(comparison)
})


test_that("Equal rates, one node fixed",{
    skip_on_cran()
    
    require(expm)
    require(MASS)
    
    data(primates)
    phy <- primates[[1]]
    phy <- multi2di(phy)
    data <- primates[[2]]
    label.vector <- rep(NA, Ntip(phy) + Nnode(phy))
    homo_gorilla <- getMRCA(phy,tip=c("Homo_sapiens", "Gorilla_gorilla"))
    label.vector[homo_gorilla] <- 2
    phy$node.label <- label.vector[-c(1:Ntip(phy))]
    
    corHMM.new <- corHMM(phy, data[,c(1,2)], model="ER", rate.cat=1, fixed.nodes=TRUE, p=0.01080903, tip.fog=0)
    corHMM.brute <- corHMM:::GetMarginalBrute(phy=phy, data=data, p=c(0.01080903,0.01080903), root.p="yang", n.states=2, node.fixed=homo_gorilla, state.fixed=2)
    comparison <- identical(sum(round(corHMM.brute - corHMM.new$states, 5)), 0)
    expect_true(comparison)
})


test_that("All rates different, no nodes fixed",{
    skip_on_cran()
    
    require(expm)
    require(MASS)
    
    data(primates)
    phy <- primates[[1]]
    phy <- multi2di(phy)
    data <- primates[[2]]
    label.vector <- rep(NA, Ntip(phy) + Nnode(phy))
    homo_gorilla <- getMRCA(phy,tip=c("Homo_sapiens", "Gorilla_gorilla"))
    label.vector[homo_gorilla] <- 2
    phy$node.label <- label.vector[-c(1:Ntip(phy))]
    
    corHMM.new <- corHMM(phy, data[,c(1,2)], model="ARD", rate.cat=1, fixed.nodes=FALSE, p=c(0.02165875,0.005681116), tip.fog=0)
    corHMM.brute <- corHMM:::GetMarginalBrute(phy=phy, data=data, p=c(0.02165875,0.005681116), root.p="yang", n.states=2, node.fixed=NULL, state.fixed=2)
    comparison <- identical(sum(round(corHMM.brute - corHMM.new$states, 5)), 0)
    expect_true(comparison)
})


test_that("All rates different, one node fixed",{
    skip_on_cran()
    
    require(expm)
    require(MASS)
    
    data(primates)
    phy <- primates[[1]]
    phy <- multi2di(phy)
    data <- primates[[2]]
    label.vector <- rep(NA, Ntip(phy) + Nnode(phy))
    homo_gorilla <- getMRCA(phy,tip=c("Homo_sapiens", "Gorilla_gorilla"))
    label.vector[homo_gorilla] <- 2
    phy$node.label <- label.vector[-c(1:Ntip(phy))]
    
    corHMM.new <- corHMM(phy, data[,c(1,2)], model="ARD", rate.cat=1, fixed.nodes=TRUE, p=c(0.02165875,0.005681116), tip.fog=0)
    corHMM.brute <- corHMM:::GetMarginalBrute(phy=phy, data=data, p=c(0.02165875,0.005681116), root.p="yang", n.states=2, node.fixed=homo_gorilla, state.fixed=2)
    comparison <- identical(sum(round(corHMM.brute - corHMM.new$states, 5)), 0)
    expect_true(comparison)
})


######################################################################################################################################
######################################################################################################################################
### Simple tests of Old corHMM vs. New corHMM
######################################################################################################################################
######################################################################################################################################


test_that("Simple tests of Old corHMM vs. New corHMM",{
  skip_on_cran()

  load("res_old.corhmm.Rsave")
  phy <- res$old.corhmm$phy
  data <- data.frame(rownames(res$old.corhmm$data), res$old.corhmm$data[,1])
  p <- c(na.omit(as.vector(res$old.corhmm$solution)))
  rate.mat <- res$old.corhmm$index.mat

  corHMM.new <- corHMM(phy, data, model="ARD", rate.cat=2, rate.mat = rate.mat, p = p, root.p = NULL, tip.fog=0)

  comparison <- identical(round(corHMM.new$loglik - res$old.corhmm$loglik, 4), 0)
  expect_true(comparison)
})

######################################################################################################################################
######################################################################################################################################
### Simple tests of corDISC vs. new corHMM
######################################################################################################################################
######################################################################################################################################

test_that("Simple tests of corDISC vs. new corHMM",{
  skip_on_cran()
    
  load("res_old.corhmm.Rsave")
  phy <- res$old.cordisc$phy
  
  data <- data.frame(sp = rownames(res$old.cordisc$data), d1 = res$old.cordisc$data[,1], d2 = res$old.cordisc$data[,2])
  p <- c(na.omit(as.vector(res$old.cordisc$solution)))
  rate.mat <- getStateMat4Dat(data, "ARD")$rate.mat

  corHMM.new <- corHMM(phy, data, model="ARD", rate.cat=1, p = p, root.p = NULL, tip.fog=0)
  
  comparison <- identical(round(corHMM.new$loglik - res$old.cordisc$loglik, 4), 0)
  expect_true(comparison)
})


######################################################################################################################################
######################################################################################################################################
### Simple tests of rayDISC vs. new corHMM
######################################################################################################################################
######################################################################################################################################

test_that("Simple tests of rayDISC vs. new corHMM",{
  skip_on_cran()
  
  load("res_old.corhmm.Rsave")
  phy <- res$old.cordisc$phy

  data <- res$old.raydisc$data
  p <- c(na.omit(as.vector(res$old.raydisc$solution)))
  rate.mat <- getStateMat4Dat(data, "ARD")$rate.mat
  
  corHMM.new <- corHMM(phy, data, model="ARD", rate.cat=1, rate.mat = rate.mat, p = p, root.p = NULL, tip.fog=0)
  
  comparison <- identical(round(corHMM.new$loglik - res$old.raydisc$loglik, 4), 0)
  expect_true(comparison)
})


######################################################################################################################################
######################################################################################################################################
### Simple tests of tip fog
######################################################################################################################################
######################################################################################################################################

test_that("Simple test of no fog vs fog estimates",{
  skip_on_cran()
  
  set.seed(1980)
  phy <- read.tree("balanced.tree")
  markov.rate <- 0.05
  Q <- rbind(c(-markov.rate, markov.rate), c(markov.rate, -markov.rate))
  ages <- 7
  phy$edge.length <- phy$edge.length*ages/max(branching.times(phy))
  traits <- corHMM:::simMarkov(phy, Q=Q, root.freqs=c(1,0))
  phy$node.label <- NULL
  phydat <- data.frame(taxon=phy$tip.label, Reg=traits$TipStates)

  corHMM.nofog <- corHMM(phy, phydat, model="ER", rate.cat=1, tip.fog=0)
  corHMM.fogest <- corHMM(phy, phydat, model="ER", rate.cat=1, tip.fog=c(1,1))

  comparison <- identical(round(corHMM.nofog$loglik - corHMM.fogest$loglik, 4), 0)
  expect_true(comparison)
})


test_that("Simple test of estimated fog vs fixed fog estimate",{
  skip_on_cran()
  
  set.seed(1980)
  phy <- read.tree("balanced.tree")
  markov.rate <- 0.05
  Q <- rbind(c(-markov.rate, markov.rate), c(markov.rate, -markov.rate))
  ages <- 7
  phy$edge.length <- phy$edge.length*ages/max(branching.times(phy))
  traits <- corHMM:::simMarkov(phy, Q=Q, root.freqs=c(1,0))
  phy$node.label <- NULL
  phydat <- data.frame(taxon=phy$tip.label, Reg=traits$TipStates)

  error <- 0.10
  error.absolute <- round(dim(phydat)[1]*error)
  taxa.sample <- sample(1:dim(phydat)[1], error.absolute)
  phydat.new <- phydat
  for(index in 1:length(taxa.sample)){
	  state <- phydat[taxa.sample[index],2]
	  if(state == 1){
		  state <- 2
	  }else{
		  state <- 1
	  }
	  phydat.new[taxa.sample[index],2] <- state
  }

  corHMM.fogest <- corHMM(phy, phydat.new, model="ER", rate.cat=1, tip.fog=c(1,1))
  corHMM.fogfixed <- corHMM(phy, phydat.new, model="ER", rate.cat=1, tip.fog=corHMM.fogest$tip.fog.probs)

  comparison <- identical(round(corHMM.fogest$loglik - corHMM.fogfixed$loglik, 4), 0)
  expect_true(comparison)
})


test_that("Test of ER model with fog vs HRM ER model with fog",{
  skip_on_cran()
  
  set.seed(1980)
  phy <- read.tree("balanced.tree")
  markov.rate <- 0.05
  Q <- rbind(c(-markov.rate, markov.rate), c(markov.rate, -markov.rate))
  ages <- 7
  phy$edge.length <- phy$edge.length*ages/max(branching.times(phy))
  traits <- corHMM:::simMarkov(phy, Q=Q, root.freqs=c(1,0))
  phy$node.label <- NULL
  phydat <- data.frame(taxon=phy$tip.label, Reg=traits$TipStates)

  error <- 0.10
  error.absolute <- round(dim(phydat)[1]*error)
  taxa.sample <- sample(1:dim(phydat)[1], error.absolute)
		  phydat.new <- phydat
		  for(index in 1:length(taxa.sample)){
			  state <- phydat[taxa.sample[index],2]
			  if(state == 1){
				  state <- 2
			  }else{
				  state <- 1
			  }
			  phydat.new[taxa.sample[index],2] <- state
  }

  rate.mat <- corHMM:::rate.mat.maker(rate.cat=2, model="ER")
  rate.mat[which(!is.na(rate.mat))] <- 1
  
  corHMM.fogest <- corHMM(phy, phydat.new, model="ER", rate.cat=1, tip.fog=c(1,1))
  corHMM.fogest.hrm <- corHMM(phy, phydat.new, model="ER", rate.cat=2, rate.mat=rate.mat, tip.fog=c(1,1))

  comparison <- identical(round(corHMM.fogest$loglik - corHMM.fogest.hrm$loglik, 4), 0)
  expect_true(comparison)
})
