# Exact probabilities by enumeration of the internal assignments of a four-tip
# tree. Tip evidence is integrated, not treated as an inferred ancestral state.
fixed_fixture <- function(rc=1L, labels=c(1,NA,NA), prior=NULL) {
  phy <- ape::read.tree(text="((a:.23,b:.91):.31,(c:.52,d:.78):1.06);")
  phy$node.label <- labels
  dat <- data.frame(taxon=letters[1:4], X=c(1L,2L,1L,1L))
  k <- 2L*rc
  idx <- matrix(NA_integer_,k,k)
  idx[row(idx)!=col(idx)] <- seq_len(k*(k-1L))
  rates <- .059 + .041*seq_len(k*(k-1L))
  Q <- matrix(0,k,k)
  Q[!is.na(idx)] <- rates[idx[!is.na(idx)]]
  diag(Q) <- -rowSums(Q)
  if(is.null(prior)) prior <- seq_len(k)/sum(seq_len(k))
  list(phy=phy, data=dat, rc=rc, k=k, idx=idx, rates=rates, Q=Q, prior=prior)
}

fixed_oracle <- function(f, joint.node=integer(), joint.state=integer()) {
  phy <- f$phy
  ntip <- length(phy$tip.label)
  assignments <- as.matrix(expand.grid(rep(list(seq_len(f$k)),phy$Nnode)))
  allowed <- rep(TRUE,nrow(assignments))
  obs <- rep(1:2,f$rc)
  for(i in which(!is.na(phy$node.label)))
    allowed <- allowed & obs[assignments[,i]] == as.integer(phy$node.label[i])
  for(i in which(joint.node>ntip))
    allowed <- allowed & assignments[,joint.node[i]-ntip] == joint.state[i]
  L <- as.numeric(allowed)
  tips <- vector("list",ntip)
  for(e in seq_len(nrow(phy$edge))) {
    parent <- assignments[,phy$edge[e,1]-ntip]
    child <- phy$edge[e,2]
    P <- expm::expm(f$Q * phy$edge.length[e])
    if(child>ntip) {
      L <- L*P[cbind(parent,assignments[,child-ntip])]
    } else {
      tip <- obs == f$data$X[child]
      for(i in which(joint.node==child)) tip <- tip & seq_len(f$k)==joint.state[i]
      weights <- sweep(P[parent,,drop=FALSE],2,tip,"*")
      mass <- rowSums(weights)
      tips[[child]] <- weights/ifelse(mass>0,mass,1)
      L <- L*mass
    }
  }
  prior <- f$prior
  if(identical(prior,"maddfitz")) {
    prior <- vapply(seq_len(f$k),function(s) sum(L[assignments[,1]==s]),0)
    prior <- prior/sum(prior)
  }
  w <- L * prior[assignments[,1]]
  marginal <- vapply(seq_len(f$k),function(s)
    colSums((assignments==s)*w)/sum(w), numeric(phy$Nnode))
  tip.probs <- t(vapply(tips,function(x) colSums(x*w)/sum(w),numeric(f$k)))
  list(loglik=log(sum(w)), marginal=marginal, tips=tip.probs,
       joint=unname(assignments[which.max(w),]))
}

fixed_fit <- function(f,method) {
  corHMM(f$phy,f$data,f$rc,rate.mat=f$idx,p=f$rates,root.p=f$prior,
         node.states=method,fixed.nodes=TRUE,get.tip.states=TRUE,verbose=FALSE)
}

test_that("ancestral constraints condition the whole reconstruction", {
  for(rc in 1:2) for(labels in list(rep(NA,3),c(1,NA,NA),c(NA,2,NA),c(1,2,NA))) {
    f <- fixed_fixture(rc,labels)
    expected <- fixed_oracle(f)
    marginal <- fixed_fit(f,"marginal")
    expect_equal(marginal$loglik,expected$loglik,tolerance=1e-10)
    expect_equal(unname(marginal$states),unname(expected$marginal),tolerance=1e-10)
    expect_equal(unname(marginal$tip.states),unname(expected$tips),tolerance=1e-10)
    joint <- fixed_fit(f,"joint")
    expect_equal(as.integer(joint$states),expected$joint)
    expect_equal(rep(1:2,rc)[joint$tip.states],f$data$X)
    scaled <- fixed_fit(f,"scaled")
    for(i in which(!is.na(labels)))
      expect_equal(sum(scaled$states[i,rep(1:2,rc)!=labels[i]]),0)
  }
  f <- fixed_fixture(2,c(NA,2,NA),"maddfitz")
  expected <- fixed_oracle(f)
  f$phy$node.label <- as.character(f$phy$node.label)
  fitted <- fixed_fit(f,"marginal")
  expect_equal(unname(fitted$states),unname(expected$marginal),tolerance=1e-10)
  expect_equal(unname(fitted$tip.states),unname(expected$tips),tolerance=1e-10)
  expect_equal(as.integer(fixed_fit(f,"joint")$states),expected$joint)
  expect_equal(ancRECON(f$phy,f$data,f$rates,method="marginal",rate.cat=f$rc,
                       rate.mat=f$idx,root.p=f$prior,get.likelihood=TRUE),
               expected$loglik,tolerance=1e-10)
  bad <- fixed_fixture(2,c(3,NA,NA))
  expect_error(fixed_fit(bad,"none"),"observed state")
  f <- fixed_fixture(1,c(1,NA,NA),c(0,1))
  expect_error(fixed_fit(f,"marginal"),"incompatible|zero probability")
  expect_error(fixed_fit(f,"joint"),"incompatible|zero probability")
})

test_that("simmap uses original node observations, not reconstructed node labels", {
  f <- fixed_fixture(2,c(NA,2,NA))
  fit <- fixed_fit(f,"marginal")
  # A fitted observed-state observation and an additional exact tip state.
  expected <- fixed_oracle(f,1L,3L)
  set.seed(143)
  histories <- makeSimmap(model=fit,fix.node=1L,fix.state=3L,nSim=256,method="uniformization")
  endpoints <- function(m, node) {
    if(node<=length(m$tip.label)) tail(names(m$maps[[which(m$edge[,2]==node)]]),1)
    else names(m$maps[[which(m$edge[,1]==node)[1]]])[1]
  }
  ids <- function(node) match(vapply(histories,endpoints,"",node=node),colnames(fit$solution))
  expect_true(all(ids(1)==3L))
  expect_true(all(ids(6) %in% c(2L,4L)))
  expect_equal(sort(unique(ids(6))),c(2L,4L))
  # A downstream observation must change the distribution upstream as well.
  freq <- tabulate(ids(5),nbins=f$k)/length(histories)
  expect_true(all(abs(freq-expected$marginal[1,]) < .11))
  expect_error(makeSimmap(model=fit,fix.node=6L,fix.state=1L),"[Ii]ncompatible|[Cc]onflict")
  expect_error(makeSimmap(model=fit,fix.node=1L,fix.state=2L),"incompatible|conflict")
})
