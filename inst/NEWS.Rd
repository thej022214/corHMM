\name{NEWS}
\title{corHMM News}
\encoding{UTF-8}

\section{CHANGES IN VERSION 2.10.5}{
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{corHMMDredge} no longer proposes structurally invalid
  models during simulated annealing. The "drop" move could remove the last
  remaining free rate, producing an all-NA index matrix; this gave a
  \code{no non-missing arguments to max} warning and an
  \code{invalid 'times' argument} error, which was caught internally and so
  did not affect the returned models. Proposals are now checked at every
  rate category: a model must retain at least one free parameter and no
  state may be isolated (NA in both its row and its column). Irreversible
  models, in which a state has rates in only one direction, remain in the
  search space. \code{corHMMDredgeBase} now also rejects a degenerate
  \code{rate.mat} with an informative error.
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Because invalid proposals are no longer generated, the random
  number stream during a dredge has changed. Runs using \code{seed} will not
  reproduce the trajectories of earlier versions. Dredges with hidden rate
  categories may also skip degenerate models they previously fit and
  accepted.
    }
  }
}

\section{CHANGES IN VERSION 2.7.1.4}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item \code{Equate/DropStateMatPars} will take empty lists (and
  return the original matrix)
    }
  }
}
