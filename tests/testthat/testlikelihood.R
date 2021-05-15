
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

    corHMM.new <- corHMM(phy, data[,c(1,2)], model="ER", rate.cat=1, fixed.nodes=FALSE, p=0.01080903)
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
    
    corHMM.new <- corHMM(phy, data[,c(1,2)], model="ER", rate.cat=1, fixed.nodes=TRUE, p=0.01080903)
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
    
    corHMM.new <- corHMM(phy, data[,c(1,2)], model="ARD", rate.cat=1, fixed.nodes=FALSE, p=c(0.02165875,0.005681116))
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
    
    corHMM.new <- corHMM(phy, data[,c(1,2)], model="ARD", rate.cat=1, fixed.nodes=TRUE, p=c(0.02165875,0.005681116))
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

  corHMM.new <- corHMM(phy, data, model="ARD", rate.cat=2, rate.mat = rate.mat, p = p, root.p = NULL)

  comparison <- identical(corHMM.new$loglik - res$old.corhmm$loglik, 0)
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

  corHMM.new <- corHMM(phy, data, model="ARD", rate.cat=1, p = p, root.p = NULL)
  
  comparison <- identical(corHMM.new$loglik - res$old.cordisc$loglik, 0)
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
  
  corHMM.new <- corHMM(phy, data, model="ARD", rate.cat=1, rate.mat = rate.mat, p = p, root.p = NULL)
  
  comparison <- identical(corHMM.new$loglik - res$old.raydisc$loglik, 0)
  expect_true(comparison)
})
