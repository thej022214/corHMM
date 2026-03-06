library(devtools)
library(bbmle)
library(MCMCpack)
library(coda)
library(lattice)
library(numDeriv)
library(GGally)

load_all("~/R/pkgs/corHMM")

data(primates)
phy <- multi2di(primates[[1]])
data <- primates[[2]]
MK_3state <- corHMM(phy = phy, data = data, rate.cat = 1)

## use returned args list
do.call(dev.corhmm,MK_3state$args.list)

## setup for computing likelihood
nllfun <- function(p) {
  a <- MK_3state$args.list
  a$p <- p
  do.call(dev.corhmm,a)
}
p <- MK_3state$args.list$p
names(p) <- parnames(nllfun) <- paste("p",1:4)
nllfun(p)

m1 <- mle2(minuslogl=nllfun, start=p, vecpar=TRUE, lower=log(1e-9), upper=log(100),
           method="L-BFGS-B")
pp <- profile(m1,trace=TRUE)
plot(pp,show.points=TRUE)

H <- hessian(nllfun, p)
sds <- sqrt(diag(solve(H)))
wald.ci <- sweep(qnorm(0.975)*outer(sds,c(-1,1)),1,FUN="+",STATS=p)
prof.ci <- confint(pp)

logpostfun <- function(p,lb=log(1e-9),ub=log(1e2),range=3) {
  prior.mean <- (lb+ub)/2
  prior.sd <- (ub-lb)/(2*range)
  loglik <- -1*nllfun(p)
  log.prior <- sum(dnorm(p,mean=prior.mean,sd=prior.sd,log=TRUE))
  if (is.na(loglik)) browser()
  return(loglik+log.prior)
}

logpostfun(p)

m1 <- MCMCmetrop1R(logpostfun,p,verbose=1000, mcmc=100000)
saveRDS(m1, file="MK_3state_mcmc.rds")
xyplot(m1)
pairs(as.matrix(m1),gap=0,pch=".")
bayes.ci <- t(apply(m1,2,quantile,c(0.025,0.975)))
summary(m1)
raftery.diag(m1)

## should do something that uses multiple chains and finds G-R statistic ...

my_mcmc <- function(data, mapping, ...) {
  ggplot(data = data, mapping=mapping) +
    ## geom_point(..., alpha = 0.2)
    geom_hex() + scale_fill_viridis_c()
}

theme_set(theme_bw())
ggpairs(as.data.frame(m1),
        lower=list(continuous=my_mcmc))

