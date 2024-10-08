\name{plot_batch_profile_lik}
\alias{plot_batch_profile_lik}
\title{Plot Batch Profile Likelihoods}
\description{
Plots profile likelihoods for each parameter in a batch from the output of a profile corHMM object, displaying the log-likelihoods across parameter values and indicating the maximum likelihood estimate (MLE) and the 95\% confidence interval.
}
\usage{
plot_batch_profile_lik( corhmm_profile, n_cols = NULL, n_rows = NULL, mar = c(5, 4, 4, 1) + 0.1, ci_level = 1.96, polygon_col = "lightgrey", line_col = "black",  line_type = "l", mle_col = "blue", ci_line_col = "black", ci_line_type = "dashed", axis_tick_length = -0.2, label_cex = 0.7, ylim=NULL)
}
\arguments{
\item{corhmm_profile}{a list containing profile likelihood tables for each parameter and the corHMM object with MLE and loglik attributes.}

\item{n_cols}{optional; number of columns for the plotting layout. If NULL, automatically calculated based on the number of parameters.}

\item{n_rows}{optional; number of rows for the plotting layout. If NULL, automatically calculated based on the number of parameters.}

\item{mar}{margins around the plot. Defaults to c(5, 4, 4, 1) + 0.1.}

\item{ci_level}{z-value for the 95\% confidence interval. Defaults to 1.96.}

\item{polygon_col}{color of the polygon under the curve. Defaults to "lightgrey".}

\item{line_col}{color of the profile likelihood curve. Defaults to "black".}

\item{line_type}{type of the profile likelihood curve. Defaults to "l".}

\item{mle_col}{color of the MLE point. Defaults to "blue".}

\item{ci_line_col}{color of the 95\% CI line. Defaults to "black".}

\item{ci_line_type}{line type for the 95\% CI line. Defaults to "dashed".}

\item{axis_tick_length}{length of the axis ticks, with negative values indicating ticks pointing inwards. Defaults to -0.2.}

\item{label_cex}{character expansion size for labels, affecting the size of text labels. Defaults to 0.7.}

\item{ylim}{a user-specified upper and lower limit to the y-axis.}
}
\details{
This function generates a series of plots for the profile likelihood of each parameter in the input corHMM model object. It visualizes the log-likelihood across the range of parameter values, highlights the maximum likelihood estimate (MLE), and denotes the 95\% confidence interval with a horizontal dashed line. The function is designed to accommodate the batch analysis of multiple parameters, organizing the plots in a specified layout and allowing for extensive customization of plot aesthetics.
}
\value{
Generates and displays a grid of profile likelihood plots, one for each parameter in the input model, with visual cues for MLE and confidence intervals.
}
\author{James D. Boyko}
\keyword{plot}