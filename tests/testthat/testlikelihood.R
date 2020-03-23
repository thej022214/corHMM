
######################################################################################################################################
######################################################################################################################################
### Simple tests of the Marginal reconstruction algorithm
######################################################################################################################################
######################################################################################################################################

test_that("Equal rates, no nodes fixed",{
    skip_on_cran()

    require(expm)
    require(corHMM)
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
    require(corHMM)
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
    require(corHMM)
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
    require(corHMM)
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


#HERE


######################################################################################################################################
######################################################################################################################################
### Simple tests of corDISC vs. new corHMM
######################################################################################################################################
######################################################################################################################################



##HERE


######################################################################################################################################
######################################################################################################################################
### Simple tests of rayDISC vs. new corHMM
######################################################################################################################################
######################################################################################################################################

#<!-- #Section 3: Unit tests of corHMM -->
#<!-- ##3.1: rayDISC-_like_ models in corHMM -->

#<!-- rayDISC is a function that runs Markov models where all states are observed. I have implemented this functionality within corHMM so that users only need to use one function when working with corHMM. Let's simulate a fresh 3-state dataset. -->

#<!-- ```{r} -->
#<!-- phy <- sim.bdtree(b = 1, d = 0.8, n = 100, stop = "taxa", extinct = FALSE) -->
#<!-- phy <- drop.extinct(phy) -->
#<!-- mean.change <- 1/sum(phy$edge.length)*20 -->

#<!-- Q <- matrix(abs(rnorm(9, mean.change, mean.change/2)), 3, 3) -->
#<!-- diag(Q) <- 0 -->
#<!-- diag(Q) <- -rowSums(Q) -->

#<!-- dat <- sim.char(phy, Q, 1, "discrete") -->
#<!-- dat <- dat[,,1] -->
#<!-- data <- data.frame(sp = names(dat), dat = dat) -->
#<!-- ``` -->

#<!-- We can run rayDISC easily.  -->

#<!-- ```{r} -->
#<!-- rayRes <- rayDISC(phy = phy, data = data, model = "ARD", root.p = "flat") -->
#<!-- ``` -->

#<!-- I will take the parameters from rayDISC and use corHMM to evaluate their likelihood. If all goes well the likelihood for these parameters will be the same in rayDISC and corHMM. -->

#<!-- ```{r} -->
#<!-- data.sort <- data.frame(data[,2], data[,2],row.names=data[,1]) -->
#<!-- data.sort <- data.sort[phy$tip.label,] -->
#<!-- model.set.final <- corHMM:::rate.cat.set.corHMM.JDB(phy=phy,data.sort=data.sort,rate.cat=1,ntraits=3, model = "ARD") -->
#<!-- p = c(na.omit(as.vector(rayRes$solution))) -->
#<!-- corRes <- corHMM(phy = phy, data = data, rate.cat = 1, rate.mat = model.set.final$rate, p = p, root.p = "flat") -->
#<!-- c(rayRes$loglik, corRes$loglik) -->
#<!-- rayRes$loglik == corRes$loglik -->
#<!-- ``` -->

#<!-- ##3.2: corHMM works for n States -->

#<!-- To test that corHMM is working for multiple states and hidden rates we need to construct a pseudo-hidden rate model. This pseudo hidden rate model will be evaluated as any rate.cat > 1 model is, but it will have the same likelihood as a model that could be run in rayDISC. The key is distributing those parameters within corHMM such that we don't expect a different likelihood.  -->

#<!-- ```{r} -->
#<!-- data.sort <- data.frame(data[,2], data[,2],row.names=data[,1]) -->
#<!-- data.sort <- data.sort[phy$tip.label,] -->
#<!-- model.set.final<- corHMM:::rate.cat.set.corHMM.JDB(phy=phy,data.sort=data.sort,rate.cat=2,ntraits=3,model ="ARD") -->
#<!-- model.set.final$index.matrix -->
#<!-- ``` -->

#<!-- This is our 2 rate class model. However, we will assign the same parameters to R1 and R2. There should be no difference in likelihood between this model and rayDISC - *we can use the rayDISC result from section 3.1.*. -->

#<!-- ```{r} -->
#<!-- p = c(na.omit(as.vector(rayRes$solution)), na.omit(as.vector(rayRes$solution)), 1,1,0) -->
#<!-- corRes_2rate <- corHMM(phy = phy, data = data, rate.cat = 2, mV = TRUE, rate.mat = model.set.final$rate, p = p, root.p = "flat") -->
#<!-- c(corRes_2rate$loglik, rayRes$loglik) -->
#<!-- ``` -->

#<!-- ##3.3: corHMMv2.0 is the same as previous versions -->

#<!-- For this section I've run some models in the previous version of corHMM. I am going to use those corhmm objects to demonstrate that new corHMM and old corHMM exist in harmony. -->

#<!-- Is old rayDISC the same as new rayDISC? -->

#<!-- ```{r} -->
#<!-- load("~/Desktop/oldCorRes.Rsave") -->
#<!-- p <- na.omit(as.vector(obj$def.rayDISC$solution)) -->
#<!-- newRay <- rayDISC(phy = obj$def.rayDISC$phy, data = obj$def.rayDISC$data, model = "ARD", ntraits = 1, p = p, rate.mat = obj$def.rayDISC$index.mat, root.p = "flat") -->
#<!-- newRay$loglik == obj$def.rayDISC$loglik -->
#<!-- ``` -->

#<!-- Is old corHMM the same as new corHMM? -->

#<!-- ```{r} -->
#<!-- p <- na.omit(as.vector(obj$def.corHMM$solution)) -->
#<!-- data <- data.frame(sp = rownames(obj$def.corHMM$data), dat = obj$def.corHMM$data[,1]) -->
#<!-- newCor <- corHMM(phy = obj$def.corHMM$phy, data = data, rate.cat = 2, rate.mat = obj$def.corHMM$index.mat, p = p) -->
#<!-- obj$def.corHMM$loglik == newCor$loglik -->
#<!-- ``` -->

#<!-- Is old corHMM the same as new corHMM when the model is constrained? -->

#<!-- ```{r} -->
#<!-- p <- na.omit(as.vector(obj$const.corHMM$solution)) -->
#<!-- p <- p[c(-4, -7)] -->
#<!-- data <- data.frame(sp = rownames(obj$const.corHMM$data), dat = obj$const.corHMM$data[,1]) -->
#<!-- newCor <- corHMM(phy = obj$const.corHMM$phy, data = data, rate.cat = 2, rate.mat = obj$const.corHMM$index.mat, p = p) -->
#<!-- obj$const.corHMM$loglik == newCor$loglik -->
#<!-- ``` -->
