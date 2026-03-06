
## compare basic pruning example
## what is different about the Q-setting machinery?
hdir <- "~/students/haoyu/pruning/R/" 
source(file.path(hdir, "pruning_funs.R"))
source(file.path(hdir, "Q_template.R"))
source(file.path(hdir, "loglik.R"))

Q2 <- model.set.final$rate
Q2[Q2==5] <- 0

traitList <- multi_to_single(primates$trait[,-1])
traitList[traitList==4] <- 3 ## collapse
Phylodata <- list(Q_template = Q2,
                  tree = phy, trait_values = traitList,
                  traitMatrix = primates$trait)
f0 <- cmb(prune_nll, Phylodata)
p0 <- list(log_trans_rates = starts)
debug(f0)
f0(p0)
ff <- MakeADFun(f0,
                list(log_trans_rates = starts), silent = TRUE)
ff$fn()


