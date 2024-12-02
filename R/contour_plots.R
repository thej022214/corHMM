# Function to generate logarithmically spaced points
generate_log_points <- function(mle, range_factor, n_points = 10){
  # mle: Maximum Likelihood Estimate of the parameter
  # range_factor: Multiplicative factor to define the range around the MLE
  # n_points: Number of points to generate
  min_val <- sapply(mle, function(x) x / range_factor)
  max_val <- sapply(mle, function(x) x * range_factor)
  log_space <- mapply(function(x,y) 
    exp(seq(log(x), log(y), length.out = n_points)), x = min_val, y = max_val)
  return(log_space)
}

fixed_corhmm <- function(par_free_0, par_fixed, par_fixed_index, 
                         dredge, pen_type, lambda, corhmm_obj){
  pars <- rep(NA, length.out = length(c(par_free_0, par_fixed)))
  pars[par_fixed_index] <- exp(par_fixed)
  pars[-par_fixed_index] <- exp(par_free_0)
  neglnLik <- compute_neglnlikelihood(pars,  corhmm_obj)
  if(dredge){
    # the dredge algorithm rescales trees to H=1, so the pars are actually H times faster
    H <- max(node.depth.edgelength(corhmm_obj$phy))
    Q = corhmm_obj$solution
    Q[] <- pars[corhmm_obj$index.mat]*H
    Q[is.na(Q)] <- 0
    diag(Q) <- -rowSums(Q)
    pen_score <- get_penalty_score(Q, pars*H, corhmm_obj$pen.type, corhmm_obj$index.mat, corhmm_obj$rate.cat)
    neglnLik <- neglnLik + pen_score
  }
  return(neglnLik)
}

optimize_fixed_corhmm <- function(par_free_0, par_fixed, par_fixed_index, 
                                  dredge, pen_type, lambda, corhmm_obj){
  optim_result <- optim(par = log(par_free_0), 
                        fn = fixed_corhmm, 
                        method = "L-BFGS-B", 
                        par_fixed = log(par_fixed), 
                        par_fixed_index=par_fixed_index, 
                        dredge=dredge,
                        pen_type=pen_type,
                        lambda=lambda,
                        corhmm_obj=corhmm_obj, 
                        lower=log(1e-10),
                        upper=log(1e10))
  return(optim_result)
}

get_profile_lik <- function(par_free_0, par_fixed_values, par_fixed_index, ncores=NULL, 
                            dredge=FALSE, pen_type=NULL, lambda=NULL, corhmm_obj=NULL){
  if(is.null(ncores)){
    optim_res_list <- lapply(par_fixed_values, function(x) 
      optimize_fixed_corhmm(par_free_0, x, par_fixed_index, dredge, pen_type, lambda, corhmm_obj))
  }else{
    ncores = min(parallel::detectCores()-1, ncores)
    optim_res_list <- parallel::mclapply(par_fixed_values, function(x) 
      optimize_fixed_corhmm(par_free_0, x, par_fixed_index, dredge, pen_type, lambda, corhmm_obj), 
      mc.cores = ncores)
  }
  profile_table <- data.frame(par_value = par_fixed_values, 
                              lnLik = -unlist(lapply(optim_res_list, "[[", "value")))
  return(list(profile_table=profile_table,
              optim_res=optim_res_list))
}


# Assuming generate_log_points is defined elsewhere and works as intended

# Wrapper function to perform profile likelihood analysis for multiple parameters
get_batch_profile_lik <- function(corhmm_obj, range_factor, n_points, verbose=FALSE, ncores=NULL, dredge=FALSE) {
  # Generate logarithmically spaced points for all parameters
  # mle_pars is expected to be a named list or vector of MLEs for each parameter
  mle_pars <- MatrixToPars(corhmm_obj)
  log_points_list <- generate_log_points(mle_pars, range_factor, n_points)
  profile_lik_results <- list()
  if(dredge){
    pen_type <- corhmm_obj$pen.type
    lambda <- corhmm_obj$lambda
  }else{
    pen_type <- NULL
    lambda <- NULL
  }
  for(i in seq_along(mle_pars)){
    if(verbose){
      cat("\n", i, "of", length(mle_pars), "...")
    }
    param_name <- names(mle_pars)[i]
    par_fixed_values <- log_points_list[, i]
    result <- get_profile_lik(mle_pars[-i], par_fixed_values, i, ncores, 
                              dredge, pen_type, lambda, corhmm_obj)
    profile_lik_results[[param_name]] <- result
  }
  profile_lik_results$corhmm_obj = corhmm_obj
  return(profile_lik_results)
}

plot_batch_profile_lik <- function(corhmm_profile, n_cols = NULL, n_rows = NULL,
                                   mar = c(5, 4, 4, 1) + 0.1, ci_level = 1.96, 
                                   polygon_col = "lightgrey", line_col = "black", line_type = "l", 
                                   mle_col = "blue", ci_line_col = "black", ci_line_type = "dashed", 
                                   axis_tick_length = -0.2, label_cex = 0.7, ylim=NULL, xlab="Parameter Value", ...){
  # Calculate the number of parameters and adjust layout if not manually specified
  n_params <- length(corhmm_profile) - 1
  if(is.null(n_cols) || is.null(n_rows)) {
    n_cols <- ceiling(sqrt(n_params))
    n_rows <- ceiling(n_params / n_cols)
  }
  layout(matrix(1:(n_rows*n_cols), nrow = n_rows, byrow = TRUE))
  par(mar = mar) # Use the user-defined margins
  
  all_lliks <- unlist(lapply(corhmm_profile[1:n_params], function(x) x$profile_table$lnLik))
  if(is.null(ylim)){
    y_min <- min(all_lliks)
    y_max <- max(all_lliks)
  }else{
    y_min <- ylim[1]
    y_max <- ylim[2]
  }
  
  mle_pars <- MatrixToPars(corhmm_profile$corhmm_obj)
  loglik <- corhmm_profile$corhmm_obj$loglik
  ci_limit = loglik - ci_level # Adjusted for user-defined CI level
  param_names = names(mle_pars)
  
  for (i in 1:n_params) {
    plot(corhmm_profile[[i]]$profile_table, type = "n", log = "x", bty = "n", axes = FALSE,
         main = bquote(theta[.(i)]), xlab = "", ylab = "", ylim = c(y_min, y_max), xaxt="n", ...)
    grid()
    axis(side = 2, las = 1, tcl = axis_tick_length, ...) # User-defined tick length
    title(xlab = xlab, ylab = "Log-Likelihood", line = 2.5, ...)
    bottom <- par("usr")[3]
    
    profile_data = corhmm_profile[[i]]$profile_table
    
    polygon(
      x = c(profile_data$par_value, rev(profile_data$par_value)),
      y = c(profile_data$lnLik, rep(bottom, length(profile_data$par_value))),
      col = polygon_col, border = NA
    )
    
    lines(profile_data$par_value, profile_data$lnLik, type = line_type, col = line_col, ...)
    
    points(mle_pars[i], loglik, pch = 19, col = mle_col)
    text(x = mle_pars[i] * 1.1, y = loglik, labels = paste("MLE =", round(mle_pars[i], 3)), pos = 4, cex = label_cex, col = line_col)
    
    abline(h = ci_limit, col = ci_line_col, lty = ci_line_type)
    
    x_range = range(profile_data$par_value)
    log_ticks = exp(seq(log(x_range[1]), log(x_range[2]), length.out = 5))
    axis(side = 1, at = log_ticks, labels = FALSE, tcl = axis_tick_length)
    labels = sapply(log_ticks, function(x) sprintf("%.1e", x))
    text(x = log_ticks, 
         y = par("usr")[3] - (par("usr")[4]-par("usr")[3])*0.05, 
         labels = labels, srt = 45, adj = 1, xpd = TRUE, cex = label_cex)
  }
  par(mar = c(5, 4, 4, 2) + 0.1) # Reset default margins
}

  # library(corHMM)
  # data(primates)
  # phy <- multi2di(primates[[1]])
  # data <- primates[[2]]
  # corhmm_fit <- corHMM(phy = phy, data = data, rate.cat = 1, root.p = "yang")
  # corhmm_fit_l1 <- corHMM:::corHMMDredge(phy = phy, data = data, 1,
  #                                        pen_type = "l1", lambda = 1, root.p = "yang")
  # 
  # corhmm_fit_l2 <- corHMM:::corHMMDredge(phy = phy, data = data, 1,
  #                                        pen_type = "l2", lambda = 1, root.p = "yang")
  # 
  # corhmm_fit_l3 <- corHMM:::corHMMDredge(phy = phy, data = data, 1,
  #                                        pen_type = "l1", lambda = 1, root.p = "maddfitz")
  # 
  # par(mfrow=c(1,2), mar = c(.1,.1,.1,.1))
  # plotRECON(corhmm_fit$phy, corhmm_fit$states, pie.cex = 1, show.tip.label = FALSE)
  # tiplabels(pie = corhmm_fit$tip.states, piecol = c("white", "black", "red"))
  # plotRECON(corhmm_fit_l3$phy, corhmm_fit_l3$states, pie.cex = 1, show.tip.label = FALSE, piecolors = c("white", "black", "red"))
  # tiplabels(pie = corhmm_fit_l3$tip.states, piecol = c("white", "black", "red"))
  # 
  # corhmm_profile <- corHMM:::get_batch_profile_lik(corhmm_obj = corhmm_fit,
  #                                                range_factor = 10000,
  #                                                n_points = 20,
  #                                                ncores = 10,
  #                                                dredge = FALSE)
  # 
  # 
  # corHMM:::plot_batch_profile_lik(corhmm_profile, ylim = c(-55, -40))
  # 
  # dredge_profile <- corHMM:::get_batch_profile_lik(corhmm_obj = corhmm_fit_dredge,
  #                                                range_factor = 10000,
  #                                                n_points = 20,
  #                                                ncores = 10,
  #                                                dredge = TRUE)
  # 
  # corHMM:::plot_batch_profile_lik(dredge_profile, ylim = c(-55, -40))

# a test profile likelihood algorithm.
# data <- rnorm(100, mean = 50, sd = 10) # 100 random normal data points
# 
# # Generate a grid of parameter values
# mu_vals <- seq(40, 60, by = 0.5)
# sigma_vals <- seq(5, 25, by = 0.5)
# 
# # Initialize a matrix to store likelihood values
# lnLikelihood_matrix <- matrix(nrow = length(mu_vals), ncol = length(sigma_vals))
# 
# # Calculate likelihood for each combination of mu and sigma
# for (i in 1:length(mu_vals)) {
#   for (j in 1:length(sigma_vals)) {
#     lnLikelihood_matrix[i, j] <- sum(dnorm(data, mean = mu_vals[i], sd = sigma_vals[j], log = TRUE))
#   }
# }
# 
# logLikelihood <- function(mu, sigma, data) {
#   sum(dnorm(data, mean = mu, sd = sigma, log = TRUE))
# }
# 
# # Profile likelihood for mu
# profile_lik_mu <- sapply(mu_vals, function(mu) {
#   optimize(f = function(sigma) -logLikelihood(mu, sigma, data), 
#            interval = c(1, 20))$objective
# })
# 
# # Profile likelihood for sigma
# profile_lik_sigma <- sapply(sigma_vals, function(sigma) {
#   optimize(f = function(mu) -logLikelihood(mu, sigma, data), 
#            interval = c(30, 70))$objective
# })
# 
# par(mfrow=c(3,1))
# # Plotting the profile likelihood for mu
# plot(mu_vals, -profile_lik_mu, type = 'l', xlab = "Mu", ylab = "-Log Likelihood")
# 
# # Plotting the profile likelihood for sigma
# plot(sigma_vals, -profile_lik_sigma, type = 'l', xlab = "Sigma", ylab = "-Log Likelihood")
# 
# image(mu_vals, sigma_vals, lnLikelihood_matrix)
# points(x = 50, y = 10, pch=21, bg="white", cex = 2)