\name{corHMM}
\alias{corHMM}
\title{Hidden Rates Model}
\description{Estimates hidden rates underlying the evolution of any number of multi-state characters}
\usage{
corHMM(phy, data, rate.cat, rate.mat=NULL, model = "ARD", 
node.states = "marginal", fixed.nodes=FALSE, p=NULL, root.p="yang", 
tip.fog=NULL, ip=NULL, nstarts=0, n.cores=1, get.tip.states = FALSE, 
lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9, upper.bound = 100, 
opts=NULL)
}
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
\item{data}{a data.frame containing species information. The first column must be species names matching the phylogeny. Additional columns contain discrete character data. All columns are converted to class factor unless already is.factor. Order of states is determined by the levels of column. Unobserved levels may be modeled if collapse=FALSE.}
\item{rate.cat}{specifies the number of hidden rate categories (if 1, there are no hidden states, see Details).}
\item{rate.mat}{a user-supplied matrix containing indexes of parameters to be optimized (see Details).}
\item{model}{One of "ARD", "SYM", or "ER". ARD: all rates differ. SYM: rates between any two states do not differ. ER: all rates are equal.} 
\item{node.states}{method used to calculate ancestral states at internal nodes (see Details).}
\item{fixed.nodes}{specifies that states for nodes in the phylogeny are assumed fixed. These are supplied as node labels in the \dQuote{phylo} object.}
\item{p}{a vector of transition rates. Allows the user to calculate the likelihood given a specified set of parameter values to specified as fixed and calculate the likelihood (see Details).}
\item{root.p}{a vector used to fix the probabilities at the root, but \dQuote{yang} and \dQuote{maddfitz} can also be supplied to use the method of Yang (2006) and FitzJohn et al (2009) respectively (see details).}
\item{tip.fog}{a fixed value or vector of free parameters to estimate the probability that an observed state is not actually in the state it is assigned. By default \code{tip.fog=NULL}. To estimate tip fog probabilities a vector of integers designating the number of free parameters is needed as input. See vignette for details on how to carry this out.}
\item{ip}{initial values used for the likelihood search. Can be a single value or a vector of unique values for each parameter. The default is \code{ip=1}.}
\item{nstarts}{the number of random restarts to be performed. The default is \code{nstarts=0}.}
\item{n.cores}{the number of processor cores to spread out the random restarts.}
\item{get.tip.states}{a boolean indicating whether tip reconstructions should be output. The default is \code{FALSE}.}
\item{lewis.asc.bias}{a boolean indicating whether to correct for observing a dataset that is not univariate. The default is \code{FALSE}}.
\item{collapse}{a boolean indicating whether to collapse multiple character combinations into only the observed states. For example, if true a two character dataset contained (0,0), (1,0), and (1,1), this would be collapsed into 1,2,3. However, if set to false it would 1,2,4. In combination with a custom rate matrix this allows for the estimation of transitions between the unobserved character combination. The default is \code{TRUE}}.
\item{lower.bound}{lower bound for the likelihood search. The default is \code{lower.bound=1e-9}.}
\item{upper.bound}{upper bound for the likelihood search. The default is \code{upper.bound=100}.}
\item{opts}{options to pass to nloptr. default is \code{NULL}.}
}
\details{
This function takes a tree and a trait file and estimates transition rates and ancestral states for any number of discrete characters using a Markov model with or without "hidden" states. Users are advised to read the _Generalized corHMM_ vignette for details on how to make full use of \code{corHMM}'s new functionality. In general, these models describe evolution as discrete transitions between observed states. If \code{ rate.cat} > 1, then the model is a hidden Markov model (HMM; also known as a hidden rates model (HRM)). The HRM is a generalization of the covarion model that allows different rate classes to be treated as "hidden" states. Essentially a hidden Markov model allows for multiple processes to describe the evolution of your observed character. This could be another (hidden) state or a large group of them. Regardless of the reason, an HMM is saying that not all observed characters are expected to act the same way. 

The first column of the input data must be species names (as in the previous version), but there can be any number of data columns. If your dataset does have 2 or more columns of trait information, each column is taken to describe a separate character. The separation of character and state is an important one because corHMM will automatically remove dual transitions from your model. For example, say you had 3 characters each with 2 states (0 or 1), but only three of these combinations were ever observed 0_0_1, 0_1_0, or 1_0_0. With dual transitions disallowed, it is impossible to move between these combinations because it would mean simultaneously losing and gaining a state (0_0_1 -> 0_0_0 -> 0_1_0 in one step.) One way around this is to provide a custom rate matrix to corHMM where transitions are allowed between these states. However, this is also a case where it would seem appropriate to code the data as a single character with 3 states.

Ambiguities (polymorphic taxa or taxa missing data) are assigned likelihoods following Felsenstein (2004, p. 255). Taxa with missing data are coded \dQuote{?} with all states observed at a tip. Polymorphic taxa are coded with states separated by an \dQuote{&}. For example, if a trait has four states and taxonA is observed to be in state 1 and 3, the character would be coded as \dQuote{1&3}. corHMM then uses this information to assign a likelihood of 1.0 to both states. Missing data are treated as ambiguous for all states, thus all states for taxa missing data are assigned a likelihood of 1.0. For example, for a four-state character (i.e. DNA), a taxon missing data will have likelihoods of all four states equal to 1.0 [e.g. L(A)=1.0, L(C)=1.0, L(G)=1.0, L(T)=1.0].

Rate matrices are usually generated via the specified \code{model}, but custom matrices can be used by supplying \code{rate.mat}. Similar to the output of \code{model}, the required format is a matrix of dimensions nTraitStates x nTraitStates, in which each transition is given an index value. Transitions with the same value will be optimised to have the same rate, so e.g. a vector with all 1's would correspond to the ER model. Transitions which are not allowed (e.g. representing two changes at once, like 00 -> 11) or not possible (e.g. the diagonal) can be set as 0 or \code{NA}.

If a user wants to calculate the likelihood for a fixed set of transition rates, without optimising them, a vector \code{p} of those values has to be provided. The positon of those values in the vector should correspond the the index numbers in \code{rate.mat}, if provided, or to those resulting from the specified \code{model}. If for any reason, any index numbers are not present in \code{rate.mat}, those positions in the vector will have to be populated with a placeholder (any number or \code{NA} will do), to prevent rates from shifting when being assigned to the transitions.

The likelihood function is maximized using the bounded subplex optimization routine implemented in the R package \code{nloptr}, which provides a common interface to NLopt, an open-source library for nonlinear optimization. In the former case, however, it is recommended that \code{nstarts} is set to a large value (e.g. 100) to ensure that the maximum likelihood solution is found. Users can set \code{n.cores} to parse the random restarts onto multiple processors.

The user can fix the root state probabilities by supplying a vector to \code{root.p}. For example, if the hypothesis is that the root is 0_S in a model with two hidden rates (their categories being denoted "S" and "F" respectively), then the root vector would be \code{root.p=c(1,0,0,0)} for state combinations 0_S, 1_S, 0_F, and 1_F, respectively. If the user supplies the flag \code{root.p}=\dQuote{NULL}, then there is equal weighting among all possible states in the  model. If the user supplies the flag \code{root.p}=\dQuote{yang}, then the estimated transition rates are used to set the weights at the root (see pg. 124 Yang 2006), whereas specifying \code{root.p}=\dQuote{maddfitz} employs the same procedure described by Maddison et al. (2007) and FitzJohn et al. (2009). Note that the default \code{root.p="yang"}.

Ancestral states can be estimated using marginal, joint, scaled, or none approaches. Marginal gives the likelihood of state at each node, integrating over the states at other nodes. Joint gives the optimal state at each node for the entire tree at once (it can only return the most likely state, i.e. it is not a probability like the marginal reconstruction). Scaled is included for compatibility with ape's ace() function. None suppresses calculation of ancestral states, which can dramatically speed up calculations if you're comparing models but make plotting difficult.
}
\value{
\code{corHMM} returns an object of class \code{corHMM}. This is a list with elements:
\item{$loglik}{the maximum negative log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{Akaike information criterion corrected for sample size.}
\item{$rate.cat}{The number of rate categories specified.}
\item{$solution}{a matrix containing the maximum likelihood estimates of the transition rates. Note that the rate classes are ordered from slowest (R1) to fastest (Rn) with respect to state 0.}
\item{$index.mat}{The indices of the parameters being estimated are returned. This also is a way to allow the estimation of transition rates for parameters not oberved in the dataset. Say you have 2 traits X and Y, where the combinations 00, 01, and 11 are observed (10 is not). A 4 by 4 index matrix could be used to force 10 into the model.}
\item{$data}{User-supplied dataset.}
\item{$data.legend}{User-supplied dataset with an extra column of trait values corresponding to how corHMM calls the user data.}
\item{$phy}{User-supplied tree.}
\item{$states}{The likeliest states at each internal node. The state and rates reconstructed at internal nodes are in the order of the column headings of the rates matrix.}
\item{$tip.states}{The likeliest state at each tip. The state and rates reconstructed at the tips are in the order of the column headings of the rates matrix.}
\item{$states.info}{a vector containing the amount of information (in bits) that the tip states and model gives to each node.}
\item{$iterations}{The number of iterations used by the optimization routine.}
\item{$root.p}{The root prior used in model estimation.}
}
\examples{
\donttest{
data(primates)
phy <- multi2di(primates[[1]])
data <- primates[[2]]
MK_3state <- corHMM(phy = phy, data = data, rate.cat = 1)
MK_3state


### custom rate.mat and user supplied transition rates not to be optimised

# this is the rate matrix used in the above example, which represents the ARD 
#model with four different rates for all four possible transitions
MK_3state$index.mat

# now we specify our own one where the transitions in and out of state 3 are 
#forced to be the same:
rate.mat <- matrix(c(NA, 1, NA, 2, NA, 3, NA, 3, NA), 3, 3)

# we supply the fixed transition rates for which we want the likelihood 
#calculated. Since we optimise 3 states in the above matrix, we have to 
#give a vector with three rates, which will be applied to the transition 
#with the corresponding rate above (e.g. the transition 2 -> 1 happens at 
#rate 0.05)
p <- c(0.05, 0.02, 0.03)

# run
MK_3state.fixed <- corHMM(phy = phy, data = data, rate.cat = 1, rate.mat = rate.mat, p = p)
MK_3state.fixed
}
}
\references{
Beaulieu J.M., B.C. O'Meara, and M.J. Donoghue. 2013. Identifying hidden rate changes in the evolution of a binary morphological character: the evolution of plant habit in campanulid angiosperms. Systematic Biology 62:725-737.

Boyko, J. D., and J. M. Beaulieu. 2021. Generalized hidden Markov models for phylogenetic comparative datasets. Methods in Ecology and Evolution 12:468-478.

Felsenstein, J. 1981. A likelihood approach to character weighting and what it tells us about parsimony and compatibility. Biological Journal of the Linnean Society 16: 183-196.

Felsenstein J. 2004. Inferring phylogenies. Sunderland MA: Sinauer Associates.

FitzJohn, R.G., W.P. Maddison, and S.P. Otto. 2009. Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Systematic Biology 58:595-611.

Maddison, W.P., P.E. Midford, and S.P. Otto. 2007. Estimating a binary characters effect on speciation and extinction. Systematic Biology 56:701-710.

Pagel, M. 1994. Detecting correlated evolution on phylogenies: a gneeral method for the comparative analysis of discrete characters. Proc. R. Soc. Lond. B 255:37-45.

Yang, Z. 2006. Computational Molecular Evolution. Oxford Press:London.
}
\author{Jeremy M. Beaulieu and James D. Boyko}
\keyword{models}
