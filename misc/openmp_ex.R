## let's try OpenMP with an RTMB function and see what happens

set.seed(102)
x <- rnorm(10000)
library(RTMB)

f <- function(pars) {
  lik <- 0
  for (i in 1:10000) {
    lik <- lik + (x[i]-pars$mu)^2
  }
  lik
}

f(list(mu=1))
f(list(mu=0))

ff <- MakeADFun(f, list(mu = 1))
ff$fn()

## ask Kasper: is this meant to be private??
## remotes::install_github("kaskr/rtmbp")
## or
install.packages('RTMBp', repos = c('https://kaskr.r-universe.dev', 'https://cloud.r-project.org'))

system.time(replicate(10000, ff$fn()))
openmp(n = 20, autopar=TRUE)
system.time(replicate(10000, ff$fn()))
