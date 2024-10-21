\name{ancRECON}
\alias{ancRECON}
\title{Ancestral state reconstruction}
\description{Infers ancestral states based on a set of model parameters}
\usage{
ancRECON(phy,data, p, method=c("joint", "marginal", "scaled"),
rate.cat, ntraits=NULL, rate.mat=NULL,
model="ARD", root.p=NULL, get.likelihood=FALSE, get.tip.states = FALSE,
tip.fog=NULL, get.info=FALSE, collapse = TRUE)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
\item{data}{a data matrix containing species information (see Details).}
\item{p}{a vector of transition rates to be used to estimate ancestral states.}
\item{method}{method used to calculate ancestral states at internal nodes. Can be one of: "joint", "marginal", or "scaled" (see Details).}
\item{rate.cat}{specifies the number of rate categories in the HRM.} 
\item{ntraits}{currently, this is automaticall detected and can always be set to NULL.}
\item{rate.mat}{a user-supplied rate matrix index of parameters to be optimized.}
\item{model}{specifies the underlying model if a rate.mat is not provided ("ER", SYM", or "ARD").}
\item{root.p}{a vector used to fix the probabilities at the root, but \dQuote{yang} and \dQuote{maddfitz} can also be supplied to use the method of Yang (2006) and FitzJohn et al (2009) respectively (see details).}
\item{get.likelihood}{a logical indicating whether to obtain the likelihood of the rates and states. The default is \code{FALSE}.}
\item{get.tip.states}{a logical indicating whether just tip reconstructions should be output. The default is \code{FALSE}.}
\item{tip.fog}{provides the probability that an observed state is not actually in the state it is assigned to the reconstruction algorithm. These values are assumed either optimized in \dQuote{corHMM} or supplied by the user.
\code{tip.fog=NULL} when using this as a standalone function.}
\item{get.info}{Whether to return information measures at nodes (Boyko and Beaulieu 2021).}
\item{collapse}{a boolean indicating whether to collapse multiple character combinations into only the observed states. For example, if true a two character dataset contained (0,0), (1,0), and (1,1), this would be collapsed into 1,2,3. However, if set to false it would 1,2,4. In combination with a custom rate matrix this allows for the estimation of transitions between the unobserved character combination. The default is \code{TRUE}}.
}
\details{
This is a stand alone function for computing the marginal, joint, or scaled likelihoods of internal nodes for a given set of transition rates. Like all other functions contained in corHMM, the tree does not have to be bifurcating in order for analyses to be carried out. IMPORTANT: If the corDISC, corHMM, and rayDISC functions are used they automatically provide a tree with the likeliest states as internal node labels. This function is intended for circumstances where the user would like to reconstruct states based on rates estimated elsewhere (e.g. BayesTraits, Mesquite, \code{ape}).

The algorithm based on Pupko et al. (2000, 2002) is used to calculate the \code{joint} estimates of ancestral states. The \code{marginal} method was originally implemented based on a description of an algorithm by Yang (2006). The basic idea is that the tree is rerooted on each internal node, with the marginal likelihood being the probabilities of observing the tips states given that the focal node is the root. However, this takes a ton of time as the number of nodes increase. But, importantly, this does not work easily when the model contains asymmetric rates. Here, we use the same dynamic programming algorithm as Mesquite (Maddison and Maddison, 2011), which is time linear with the number of species and calculates the marginal probability at a node using an additional up and down pass of the tree.  If \code{scaled}, the function uses the same algorithm from ace(). Note that the \code{scaled} method of ace() is simply the conditional likelihoods of observing everything at or above the focal node and these should NOT be used for ancestral state estimation.

The user can fix the root state probabilities by supplying a vector to \code{root.p}. For example, in the two trait case, if the hypothesis is that the root is 00, then the root vector would be \code{root.p=c(1,0,0,0)} for state combinations 00, 01, 10, and 11, respectively. If analyzing a binary or multistate character, the order of root.p is the same order as the traits -- e.g., for states 1, 2, 3, a \code{root.p=c(0,1,0)} would fix the root to be in state 2. If the user supplies the flag \code{root.p}=\dQuote{yang}, then the estimated transition rates are used to set the weights at the root (see pg. 124 Yang 2006), whereas specifying \code{root.p}=\dQuote{maddfitz} employs the same procedure described by Maddison et al. (2007) and FitzJohn et al. (2009). Note that the default \code{root.p=NULL} assumes equal weighting among all possible states. 

Setting get.likelihood=TRUE will provide the user the joint likelihood of the rates and states.
}
\value{
\item{$lik.tip.states}{A matrix of the reconstructed tip values. If the number of rate.cats is greater than 2 then the probability that each observed state is in a particular hidden state is given.}
\item{$lik.anc.states}{For \code{joint}, a vector of likeliest states at internal nodes and tips. For either \code{marginal} or \code{$scaled}, a matrix of the probabilities of each state for each internal node are returned.}
\item{$info.anc.states}{A vector containing the amount of information (in bits) that the tip states and model gives to each node. See Boyko and Beaulieu (2021).}
}
\examples{
\donttest{
data(primates)
phy <- multi2di(primates[[1]])
data <- primates[[2]]
MK_3state <- corHMM(phy = phy, data = data, rate.cat = 1)

# # one way to get the parameters from your corHMM object in the correct order
p <- sapply(1:max(MK_3state$index.mat, na.rm = TRUE), function(x) 
	na.omit(c(MK_3state$solution))[na.omit(c(MK_3state$index.mat) == x)][1])

# using custom params
states_1 <- ancRECON(phy = phy, data = MK_3state$data, p = p, method = "marginal", 
	rate.cat <- MK_3state$rate.cat, ntraits = NULL, rate.mat = MK_3state$index.mat, 
	root.p = MK_3state$root.p)

}
}
\references{
FitzJohn, R.G., W.P. Maddison, and S.P. Otto. 2009. Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Systematic Biology 58:595-611.

Maddison, W.P. and D.R. Maddison. 2011. Mesquite: a modular system for evolutionary analysis. Version 2.75 http://mesquiteproject.org

Pupko, T., I. Pe'er, R. Shamir, and D. Graur. 2000. A fast algorithm for joint reconstruction of ancestral amino-acid sequences. Molecular Biology and Evolution 17:890-896.

Pupko, T., I. Pe'er, D. Graur, M. Hasegawa, and N Friedman N. 2002. A branch-and-bound algorithm for the inference of ancestral amino-acid sequences when the replacement rate varies among sites: application to the evolution of five gene families. Bioinformatics 18:1116-1123.

Yang, Z. 2006. Computational Molecular Evolution. London:Oxford.

Boyko, J. D., and J. M. Beaulieu. 2021. Generalized hidden Markov models for phylogenetic comparative datasets. Methods in Ecology and Evolution 12:468-478.
}

\author{Jeremy M. Beaulieu and Jeffrey C. Oliver}
\keyword{reconstructions}

