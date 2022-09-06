#' Sample points from along a ridge
#' This "dents" the likelihood surface by reflecting points better than a threshold back across the threshold (think of taking a hollow plastic model of a mountain and punching the top so it's a volcano). It then uses essentially a Metropolis-Hastings walk to wander around the new rim. It adjusts the proposal width so that it samples points around the desired likelihood. 
#' This is better than using the curvature at the maximum likelihood estimate since it can actually sample points in case the assumptions of the curvature method do not hold. It is better than varying one parameter at a time while holding others constant because that could miss ridges: if I am fitting 5=x+y, and get a point estimate of (3,2), the reality is that there are an infinite range of values of x and y that will sum to 5, but if I hold x constant it looks like y is estimated very precisely. Of course, one could just fully embrace the Metropolis-Hastings lifestyle and use a full Bayesian approach.
#' 
#' While running, it will display the current range of likelihoods in the desired range (by default, the best negative log likelihood + 2 negative log likelihood units) and the parameter values falling in that range. If things are working well, the range of values will stabilize during a search.
#' 
#' The algorithm tunes: if it is moving too far away from the desired likelihoods, it will decrease the proposal width; if it staying in areas better than the desired likelihood, it will increase the proposal width. It will also expand the proposal width for parameters where the extreme values still appear good enough to try to find out the full range for these values. 
#' 
#' In general, the idea of this is not to give you a pleasingly narrow range of possible values -- it is to try to find the actual uncertainty, including finding any ridges that would not be seen in univariate space.
#' @param par Starting parameter vector, generally at the optimum. If named, the vector names are used to label output parameters.
#' @param fn The likelihood function, assumed to return negative log likelihoods
#' @param best_neglnL The negative log likelihood at the optimum; other values will be greater than this.
#' @param delta How far from the optimal negative log likelihood to focus samples
#' @param nsteps How many steps to take in the analysis
#' @param print_freq Output progress every print_freq steps.
#' @param lower_bound Minimum parameter values to try. One for all or a vector of the length of par.
#' @param upper_bound Maximum parameter values to try. One for all or a vector of the length of par.
#' @param adjust_width_interval When to try automatically adjusting proposal widths
#' @param badval Bad negative log likelihood to return if a non-finite likelihood is returned
#' @param sd_vector Vector of the standard deviations to use for proposals. Generated automatically if NULL
#' @param restart_after Sometimes the search can get stuck outside the good region but still accept moves. After this many steps without being inside the good region, restart from one of the past good points
#' @param debug If TRUE, prints out much more information during a run
#' @param ... Other arguments to fn. 
#' @return A dentist object containing results, the data.frame of negative log likelihoods and the parameters associated with them; acceptances, the vector of whether a proposed move was accepted each step; best_neglnL, the best value passed into the analysis; delta, the desired offset; all_ranges, a summary of the results.
#' @export
#' @examples
#' # Univariate case
#' sims <- stats::rnorm(100, mean=17)
#' possible_means <- seq(from=16, to=18, length.out=100) # for optimize
#' 
#' # Make sure we have a function that takes in a parameters vector, other arguments if needed,
#' # and returns the negative log likelihood
#' dnorm_to_run <- function(par, sims) {
#'   return(-sum(stats::dnorm(x=sims, mean=par, log=TRUE)))
#' }
#' 
#' optimized_results <- stats::optimize(dnorm_to_run,interval=range(possible_means), 
#'   sims=sims, maximum=FALSE)
#' best_par <- optimized_results$minimum
#' names(best_par) <- "mean"
#' best_neglnL <- optimized_results$objective
#' 
#' dented_results <- dent_walk(par=best_par, fn=dnorm_to_run, best_neglnL=best_neglnL, sims=sims)
#' plot(dented_results)
#' 
#' # Multivariate case
#' sims <- stats::rlnorm(100, meanlog=1, sdlog=3)
#' 
#' dlnorm_to_run <- function(par, sims) {
#'   return(-sum(stats::dlnorm(sims, meanlog=par[1], sdlog=par[2], log=TRUE)))
#' }
#' 
#' optimized_results <- stats::optim(c(meanlog=.9, sdlog=2.9), dlnorm_to_run, sims=sims)
#' best_par <- optimized_results$par
#' best_neglnL <- optimized_results$value
#' 
#' dented_results <- dent_walk(par=best_par, fn=dlnorm_to_run, best_neglnL=best_neglnL, sims=sims)
#' plot(dented_results)
dent_walk <- function(par, fn, best_neglnL, delta=2, nsteps=1000, print_freq=50, lower_bound=0, upper_bound=Inf, adjust_width_interval=100, badval=1e9, sd_vector=NULL, debug=FALSE, restart_after=50, ...) {
  results <- data.frame(matrix(NA, nrow=nsteps+1, ncol=length(par)+1))
  results[1,] <- c(best_neglnL, par)
  if(is.null(sd_vector[1])) {
    sd_vector <- 0.1*abs(par)
  }
  sd_vector_positive <- sd_vector[which(sd_vector>0)]
  sd_vector[sd_vector==0] <- 0.5*min(sd_vector_positive)
  acceptances <- rep(NA, nsteps)
  old_params <- par
  old_dented_neglnL <- dent_likelihood(best_neglnL, best_neglnL, delta)
  rep_index <- 0
  nsteps_original <- nsteps
  steps_since_in_region <- 0
  while (rep_index <= nsteps) {
    rep_index <- rep_index+1
    if(steps_since_in_region>restart_after) {
      if(debug) {
        print("Have not been inside target region, so starting again by sampling a random point within there")
      }
      steps_since_in_region <- 0
      good_enough <- which(results[1:rep_index,1]<=best_neglnL+delta)
      chosen_good <- sample(good_enough, 1)
      param_names <- names(old_params)
      old_params <- as.numeric(results[chosen_good, -1])
      names(old_params) <- param_names
      old_dented_neglnL <- dent_likelihood(results[chosen_good,1], best_neglnL, delta)
    }
    
    new_params <- dent_propose(old_params, lower_bound=lower_bound, upper_bound=upper_bound, sd=sd_vector) 
    
    new_neglnL <- fn(par=new_params, ...)
    
    if(!is.finite(new_neglnL)) {
      new_neglnL <- max(badval, 10+10*abs(best_neglnL))
    }
    
    if((new_neglnL-best_neglnL)>delta) {
      steps_since_in_region <- steps_since_in_region+1	
    } else {
      steps_since_in_region <- 0	
    }
    
    if((new_neglnL-best_neglnL)<(-0.01)) {
      warning(paste0("Found an undented likelihood that was ", round(best_neglnL-new_neglnL,2), " log likelihood units BETTER than the best likelihood. Now using this as the new best, but suggests your initial search was flawed. Doing more steps starting from here."), immediate.=TRUE)
      best_neglnL <- new_neglnL
      old_dented_neglnL <- dent_likelihood(best_neglnL, best_neglnL, delta)
      nsteps <- nsteps_original+nsteps
    }
    new_dented_neglnL <- dent_likelihood(new_neglnL, best_neglnL, delta)
    
    
    results[rep_index+1,] <- c(new_neglnL, new_params)
    if(debug) {
      print(results[rep_index+1,] )
      print(paste0("old_dented_neglnL: ", old_dented_neglnL, " new_dented_neglnL: ", new_dented_neglnL))
    }
    if(new_dented_neglnL<=old_dented_neglnL) {
      old_params <- new_params
      acceptances[rep_index]<-TRUE
      old_dented_neglnL <- new_dented_neglnL
    } else {
      if(new_dented_neglnL-old_dented_neglnL < stats::runif(1)) {
        old_params <- new_params
        acceptances[rep_index]<-TRUE
        old_dented_neglnL <- new_dented_neglnL
      } else {
        acceptances[rep_index]<-FALSE
      }
    }
    
    if(rep_index%%adjust_width_interval==0) { # adaptively change proposal width for all params at once
      sd_original <- sd_vector
      acceptances_run <- utils::tail(acceptances[!is.na(acceptances)],adjust_width_interval)
      if(sum(acceptances_run)/length(acceptances_run) > 0.3) {
        sd_vector <- sd_vector * 1.5
        print("increasing proposal width for all parameters")
        if(debug) {
          print("changed proposals")
          print(data.frame(old=sd_original, new=sd_vector))	
        }
      }
      if(sum(acceptances_run)/length(acceptances_run) < 0.1) {
        
        sd_vector <- sd_vector * 0.8
        print("decreasing proposal width for all parameters")
        if(debug) {
          print("changed proposals")
          print(data.frame(old=sd_original, new=sd_vector))	
        }
      }     
    }
    if(rep_index%%(2*adjust_width_interval)==0) { # if we haven't found values of some parameters that are outside the CI, widen the search for just those
      sd_original <- sd_vector
      good_enough_results_range <- apply(results[which(results[1:(rep_index+1),1]-min(results[1:(rep_index+1),1])<=delta),], 2, range)[,-1]
      all_results_range <- apply(results[1:(rep_index+1),], 2, range)[,-1]
      not_past_bounds <- rep(FALSE, ncol(results)-1)
      if(!is.null(dim(all_results_range))) {
        not_past_bounds <- apply(all_results_range==good_enough_results_range, 2, any)
      } else {
        not_past_bounds <- any(all_results_range==good_enough_results_range) # single parameter
      }
      if(any(not_past_bounds)) {
        print("increasing proposal width for some parameters")
        sd_vector[not_past_bounds] <- sd_vector[not_past_bounds]*1.5
        if(debug) {
          print("changed proposals")
          print(data.frame(old=sd_original[not_past_bounds], new=sd_vector[not_past_bounds]))	
        }
      }
    }
    if(rep_index%%print_freq==0) {
      print(paste("Done replicate",rep_index))
      intermediate_results <- apply(results[which(results[1:(rep_index+1),1]-min(results[1:(rep_index+1),1])<=delta),], 2, range)
      print(paste0("CI of values (the ", length(which(results[1:(rep_index+1),1]-min(results[1:(rep_index+1),1])<=delta)), " replicates within ", delta, " neglnL of the optimum)"))
      try(colnames(intermediate_results) <- c("neglnL", names(par)))
      print(intermediate_results)
    }
  }
  colnames(results) <- c("neglnL", names(par))
  CI_range <- apply(results[which(results[1:(rep_index+1),1]-min(results[1:(rep_index+1),1])<=delta),], 2, range)[,-1]
  
  total_range <- apply(results, 2, range)[,-1]
  best_vals <- results[which.min(results[,1]),][,-1]
  if(is.null(dim(CI_range))) { #single parameter, so rearrange
    CI_range <- t(t(CI_range))
    total_range <- t(t(total_range))
  }
  all_ranges <- rbind(best_vals, CI_range, total_range)
  rownames(all_ranges) <- c("best", "lower.CI", "upper.CI", "lowest.examined", "highest.examined")
  if(min(results[,1])<best_neglnL) {
    warning(paste0("The best negative log likelihood found during sampling was ",  best_neglnL-min(results[,1]), " negative log likelihood units BETTER than the starting best value. This can indicate your original optimization search failed to find the global optimum, especially if this number is large (>0.1)"))
  }
  
  final_results <- list(results=results, acceptances=acceptances, best_neglnL=best_neglnL, delta=delta, all_ranges=all_ranges)
  class(final_results) <- c("dentist", "list")
  return(final_results)
}

#' Propose new values
#' This proposes new values using a normal distribution centered on the original parameter values, with desired standard deviation. If any proposed values are outside the bounds, it will propose again.
#' @param old_params The original parameter values
#' @param lower_bound Minimum parameter values to try. One for all or a vector of the length of par.
#' @param upper_bound Maximum parameter values to try. One for all or a vector of the length of par.
#' @param sd Standard deviation to use for the proposals. One for all or a vector of the length of par.
#' @return A vector of the new parameter values
dent_propose <- function(old_params, lower_bound=-Inf, upper_bound=Inf, sd=1) {
  sd <- abs(sd)
  if(runif(1)<0.1) { #try changing all
    new_params <- stats::rnorm(length(old_params), old_params, sd)
  } else { #try sampling some but not all. Usually just one.
    new_params <- old_params
    focal <- sample.int(length(old_params),min(length(old_params), ceiling(stats::rexp(1, 1/2))))
    new_params[focal] <- stats::rnorm(1, old_params[focal], ifelse(length(sd)==1,sd, sd[focal]))  
  }
  while(any(new_params<lower_bound) | any(new_params>upper_bound)) {
    sd <- sd*0.1
    new_params <- dent_propose(old_params, lower_bound=lower_bound, upper_bound=upper_bound, sd=sd)
  }
  return(new_params)
}

#' Dents the likelihood surface
#' This takes any values that are better (lower) than the desired negative log likelihood and reflects them across the best_neglnL + delta line, "denting" the likelihood surface.
#' @param neglnL The original negative log likelihood
#' @param best_neglnL The negative log likelihood at the optimum; other values will be greater than this.
#' @param delta How far from the optimal negative log likelihood to focus samples
#' @return The transformed negative log likelihood
#' @export
#' @examples
#' sims <- stats::rnorm(100, mean=17)
#' possible_means <- seq(from=16, to=18, length.out=100)
#' results_normal <- rep(NA, length(possible_means))
#' results_dented <- results_normal
#' dnorm_to_run <- function(par, sims) {
#'   return(-sum(dnorm(x=sims, mean=par, log=TRUE)))
#' }
#' best_neglnL <- optimize(dnorm_to_run,interval=range(possible_means), 
#'   sims=sims, maximum=FALSE)$objective
#' for (i in seq_along(possible_means)) {
#'   results_normal[i] <- dnorm_to_run(possible_means[i], sims=sims)
#'   results_dented[i] <- dent_likelihood(results_normal[i], best_neglnL=best_neglnL, delta=2)
#' }
#' 
#' plot(possible_means, results_normal)
#' lines(possible_means, results_dented)
dent_likelihood <- function(neglnL, best_neglnL, delta=2) {
  difference <- neglnL - (best_neglnL+delta)
  if(difference<0) {
    neglnL <- (best_neglnL+delta) - difference
  }
  return(neglnL)
}

#' Plot the dented samples
#' This will show the univariate plots of the parameter values versus the likelihood as well as bivariate plots of pairs of parameters to look for ridges.
#' @param x An object of class dentist
#' @param ... Other arguments to pass to plot
#' @export
plot.dentist <- function(x, ...) {
  nparams <- ncol(x$results)-1
  nplots <- nparams + (nparams^2 - nparams)/2
  results <- x$results
  threshold <- x$best_neglnL + x$delta
  results$color <- ifelse(results[,1]<=threshold, "black", "gray")
  results_outside <- subset(results, results$color=="gray")
  results_inside <- subset(results, results$color=="black")
  graphics::par(mfrow=c(ceiling(nplots/nparams), nparams))
  for (i in sequence(nparams)) {
    plot(results[,i+1], results[,1], pch=20, col=results$color, main=colnames(results)[i+1], xlab=colnames(results)[i+1], ylab="Negative Log Likelihood", bty="n", ...)
    graphics::abline(h=threshold, col="blue")
    graphics::points(results[which.min(results[,1]), i+1], results[which.min(results[,1]),1], pch=21, col="red")
  }
  for (i in sequence(nparams)) {
    for (j in sequence(nparams)) {
      if(j>i) {
        plot(results_outside[,i+1], results_outside[,j+1], pch=20, col=results_outside$color, xlab=colnames(results)[i+1], ylab=colnames(results)[j+1], bty="n", main=paste0(colnames(results)[j+1], " vs. ", colnames(results)[i+1]), ...)
        graphics::points(results_inside[,i+1], results_inside[,j+1], pch=20, col=results_inside$color)
        graphics::points(results[which.min(results[,1]), i+1], results[which.min(results[,1]),j+1], pch=21, col="red")
      }
    }	
  }
}

#' Summarize dentist
#' Display summary of output from dent_walk
#' @param object An object of class dentist
#' @param ... Other arguments (not used)
#' @export
summary.dentist <- function(object, ...) {
  cat(paste0("This ran ", nrow(object$results)-1, " steps looking for all points within ", object$delta, " negative log likelihood units of the best parameter values.\n"))
  cat("\nParameters: \n")
  print(object$all_ranges)
}

#' Print dentist
#' print summary of output from dent_walk
#' @param x An object of class dentist
#' @param ... Other arguments (not used)
#' @export
print.dentist <- function(x, ...) {
  summary.dentist(x,...)	
}
