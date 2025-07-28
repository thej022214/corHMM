\name{plotDredgeTrace}
\alias{plotDredgeTrace}
\title{Plot the Trace of a corHMMDredge Simulated Annealing Search}
\description{
This function visualizes the search history of a \code{\link{corHMMDredge}} run. It creates a trace plot showing the information criterion score at each iteration, distinguishes between accepted and rejected models, and displays the type of move proposed at each step.
}
\usage{
plotDredgeTrace(dredge_fits, break_size = 5, 
                palette = c("drop" = "#A23B72", "merge" = "#2E86AB", 
                            "free" = "#F18F01", "restart" = "#7209B7", 
                            "none" = "grey60"), 
                legend = TRUE, legend.pos = "topright", ...)
}
\arguments{
  \item{dredge_fits}{An object of class \code{corhmm.dredge} returned by \code{\link{corHMMDredge}} with the argument \code{return.all=TRUE}.}
  \item{break_size}{An integer specifying the size of the empty gap in the plot between the traces of different rate categories.}
  \item{palette}{A named vector of colors for the different move types ("drop", "merge", "free", "restart", "none"). The function will add default colors for any unpecified move types found in the data.}
  \item{legend}{Logical. If TRUE, a legend is added to the plot.}
  \item{legend.pos}{A character string specifying the position of the legend (e.g., "topright", "bottomleft").}
  \item{...}{Additional arguments to be passed to the main \code{plot} function, such as \code{main}, \code{xlab}, \code{ylab}, etc.}
}
\details{
The plot provides a comprehensive overview of the simulated annealing search process. The main plot area shows the information criterion score (e.g., AIC) on the y-axis against the iteration number on the x-axis. Each point represents a model that was evaluated.

\itemize{
  \item \strong{Green circles} denote proposed models that were accepted.
  \item \strong{Red 'x's} denote proposed models that were rejected.
  \item A \strong{grey line} connects the sequence of accepted models to show the path of the search.
  \item A \strong{colored bar} at the bottom of the plot indicates the type of move proposed at each iteration ("drop", "merge", "free", or "restart").
  \item If the dredge search was performed over multiple rate categories (i.e., \code{max.rate.cat > 1}), the traces for each rate category run are separated by a \strong{dashed vertical line}.
}

This visualization is useful for diagnosing the performance of the search, such as assessing whether the temperature and cooling schedule were appropriate, or if the search became trapped in a local optimum.
}
\author{
James D. Boyko
}
\seealso{
\code{\link{corHMMDredge}}
}