\name{plot_transition_summary}
\alias{plot_transition_summary}

\title{
Plot Summary Statistics of Stochastic Character Map Transitions
}

\description{
Creates a two-panel plot to visualize the results from a list of stochastic character mapping simulations. The function generates a consolidated boxplot for transition counts and a staggered density plot (ridgeline plot) for transition ages.
}

\usage{
plot_transition_summary(simmap_summaries, cols = NULL)
}

\arguments{
  \item{simmap_summaries}{A list of data frames, where each data frame is the result of applying \code{summarize_single_simmap} to a "simmap" object.}
  \item{cols}{An optional vector of character strings specifying colors for the plots. If \code{NULL}, a default color palette is used. Colors are recycled if there are more transition types than colors provided.}
}

\details{
This function produces a visualization of simulation results from stochastic character mapping. It is designed to be used after processing a list of "simmap" objects with \code{lapply(simmaps, summarize_single_simmap)}.

The graphical output consists of two main panels:
\itemize{
  \item \strong{Left Panel:} A single plot containing vertical boxplots for the count of each transition type, summarized across all simulations. Semi-transparent points are overlaid with a small amount of jitter to show the distribution of individual data points.
  \item \strong{Right Panel:} A ridgeline plot showing the distribution of transition ages.
}

}
\examples{
\dontrun{
data(primates)
phy <- primates[[1]]
phy <- multi2di(phy)
data <- primates[[2]]

##run corhmm
MK <- corHMM(phy, data, 1)

##get simmap from corhmm solution
model <- MK$solution
simmaps <- makeSimmap(tree=phy, data=data, model=model, rate.cat=1, nSim=100, nCores=1)

# get the summary
simmap_summaries <- lapply(simmaps, summarize_single_simmap)
summary_df <- summarize_transition_stats(simmap_summaries)
print(summary_df)

plot_transition_summary(simmap_summaries)
}
}
