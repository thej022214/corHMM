\name{rate.mat.maker}
\alias{rate.mat.maker}
\alias{rate.par.drop}
\alias{rate.par.eq}
\title{Rate matrix maker}
\description{Generates and manipulates the index of the rate parameters to be optimized}
\usage{
rate.mat.maker(rate.cat, hrm=TRUE, ntraits=NULL, nstates=NULL, 
model=c("ER", "SYM", "ARD"))
rate.par.drop(rate.mat.index=NULL,drop.par=NULL)
rate.par.eq(rate.mat.index=NULL, eq.par=NULL)

}       
\arguments{
\item{rate.cat}{specifies the number of rate categories in the HRM.}
\item{hrm}{a logical indicating whether the underlying model is the hidden rates model (HRM). The default is \code{FALSE}.} 
\item{ntraits}{specifies the number of traits in the data file if the underlying model is not the HRM.}
\item{nstates}{specifies the number of characters in the data file used in rayDISC.}
\item{model}{if the model is not HRM, specifies the underlying model.}
\item{rate.mat.index}{A user-supplied rate matrix index to be manipulated.}
\item{drop.par}{a vector of transitions to be dropped from the model. Use \code{rate.mat.index} to see what correspond to which transition.}
\item{eq.par}{a vector of transitions pairs to be set equal. Use \code{rate.mat.index} to see what correspond to which transition.}
}
\details{
Outputs the full index of the rate parameters that are to be optimized. The intention is that a user might want to see how the matrix is designed prior to an analysis and perhaps drops a few parameters beforehand due to some hypothesis that he or she might have. The resulting matrix can then be plugged directly into corHMM, corDISC, or rayDISC.
}
\value{
Returns a rate matrix index
}
\examples{
#Generate a matrix for two binary traits:
rate.mat<-rate.mat.maker(hrm=FALSE,ntraits=2,model="ARD")
#Drop parameter 8 from the model
rate.mat<-rate.par.drop(rate.mat, drop.par=c(8))
#Set parameters 1 and 2 equal to one another:
rate.mat<-rate.par.eq(rate.mat, eq.par=c(1,2))

#Precursor model. There are many ways to do this, but here is one way 
rate.mat<-rate.mat.maker(hrm=TRUE,rate.cat=2)
rate.mat<-rate.par.drop(rate.mat,c(1,3,4,6,7,8))
rate.mat<-rate.par.eq(rate.mat,c(1,2))
#Now add in a couple more connections:
rate.mat[3,2]<-1
rate.mat[2,3]<-1
#Now just use this matrix when using the corHMM function

#Here is a one way of doing a more complicated precursor:
rate.mat[3,2]<-2
rate.mat[1,3]<-3
rate.mat[2,3]<-4
#Again, just use this matrix when using the corHMM function

#Finally, here is an easier way of doing the precursor:
rate.mat<-rate.mat.maker(hrm=TRUE,rate.cat=2)
rate.mat<-rate.par.drop(rate.mat,c(1,3,4,7))
rate.mat[!is.na(rate.mat)]<-1

#Not run
# pp<-corHMM(primates$tree,primates$trait,rate.cat=2,rate.mat=rate.mat,
# node.states="marginal",diagn=FALSE)

}
\author{Jeremy M. Beaulieu and Jeffrey C. Oliver}
\keyword{models}
