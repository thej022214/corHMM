### corHMM -- Generalized hidden Markov Models
# automatic fitting 
corHMMDredge <- function(phy, data, max.rate.cat=1, init.rate.cat=1, 
  root.p="maddfitz", pen.type = "l1", lambda = 0, drop.threshold = 1e-7, 
  criterion="AIC", merge.threshold=0, index_mat=NULL, node.states = "marginal", 
  fixed.nodes=FALSE, ip=NULL, nstarts=0, n.cores=1, get.tip.states = FALSE, 
  lewis.asc.bias = FALSE, collapse = FALSE, lower.bound = 1e-10, 
  upper.bound = 100, opts=NULL, verbose=TRUE, p=NULL, rate.cat=NULL, grad=FALSE, 
  max.iterations = 200, initial.temp = 2, cooling.rate = 0.95, 
  temp.schedule = "exponential", seed = NULL, return.all=FALSE) {
  
  if((is.null(p) & !is.null(rate.cat))){
    print("A rate category was given without specifying a parameter vector (p)")
    return(NULL)
  }
  if((!is.null(p) & is.null(rate.cat))){
    print("A parameter vector (p) was given without specifying a rate category")
    return(NULL)
  }
  
  # Set seed for reproducibility
  if(!is.null(seed)) set.seed(seed)
  
  # FIXED FIT
  init.root.p <- root.p
  if(!is.null(p) & !is.null(rate.cat)){
    if(verbose){
      cat("Evaluating fixed parameters p =", p, "\n")
    }
    fixd_fit <- corHMMDredgeBase(phy=phy, data=data, rate.cat=rate.cat, 
      root.p=root.p, pen.type = pen.type, lambda = lambda, rate.mat=rate.mat, 
      node.states = node.states, fixed.nodes=fixed.nodes, ip=ip, 
      nstarts=nstarts, n.cores=n.cores, get.tip.states = get.tip.states, 
      lewis.asc.bias = lewis.asc.bias, collapse = collapse, 
      lower.bound = lower.bound, upper.bound = upper.bound, 
      opts=opts, p=p, grad=FALSE)
    return(fixd_fit)
  }
  
  # Initialize
  if(is.null(index_mat)){
    curr_index_mat <- getStateMat4Dat(data, collapse = FALSE, indep = FALSE)$rate.mat
    max_index_mat <- curr_index_mat
    max_index_mat[max_index_mat == 0] <- NA
    max_index_mat_cp <- max_index_mat
  } else {
    curr_index_mat <- getStateMat4Dat(data, collapse = FALSE, indep = FALSE)$rate.mat
    max_index_mat <- curr_index_mat
    max_index_mat[max_index_mat == 0] <- NA
    max_index_mat_cp <- max_index_mat
    curr_index_mat <- index_mat
  }
  
  fit_set <- list()
  model_improved <- TRUE
  hmm_valid <- TRUE
  count <- 0
  current_rate_category <- init.rate.cat
  
  if(init.rate.cat > 1){
    max_index_mat <- getFullMat(replicate(current_rate_category, max_index_mat_cp, simplify = FALSE), 
      getStateMat(current_rate_category))
    max_index_mat[max_index_mat > 0] <- 1:sum(max_index_mat > 0, na.rm = TRUE)
    max_index_mat[max_index_mat==0] <- NA
  }
  
  if(verbose){
    cat("Beginning SA dredge...\n")
    cat("SA parameters: max_iter =", max.iterations, ", init_temp =", initial.temp, 
      ", cooling =", cooling.rate, "\n")
  }
  
  # Main loop over rate categories
  while(model_improved){
    count <- count + 1
    if(!inherits(root.p, "character")){
      root.p <- rep(init.root.p, current_rate_category)
    }
    
    if(verbose){
      cat("\n=== RATE CATEGORY", current_rate_category, "===\n")
    }
    
    # Initial fit for this rate category
    curr_fit <- try(corHMMDredgeBase(phy=phy, data=data, 
      rate.cat=current_rate_category, root.p=root.p, pen.type = pen.type, 
      lambda = lambda, rate.mat=curr_index_mat, node.states = node.states, 
      fixed.nodes=fixed.nodes, ip=ip, nstarts=nstarts, 
      n.cores=n.cores, get.tip.states = get.tip.states, 
      lewis.asc.bias = lewis.asc.bias, collapse = collapse, 
      lower.bound = lower.bound, upper.bound = upper.bound, 
      opts=opts, p=NULL, grad=grad))
    
    if(inherits(curr_fit, "try-error")){
      warning("Model fitting failed. Stopping dredge.")
      model_improved <- FALSE
      fit_set[[count]] <- curr_fit
      next
    }
    
    curr_info_criterion <- curr_fit[[criterion]]
    if(verbose){
      cat("Initial", criterion, ":", curr_info_criterion, "\n")
      print(curr_fit$index.mat)
      cat("\n")
    }
    
    # SIMULATED ANNEALING WITHIN THIS RATE CATEGORY
    sa_result <- sa_within_rate_category(phy, data, curr_fit, 
      curr_index_mat, max_index_mat,
      current_rate_category, root.p, 
      pen.type, lambda, node.states, fixed.nodes, 
      ip, nstarts, n.cores, get.tip.states, 
      lewis.asc.bias, collapse, lower.bound, 
      upper.bound, opts, grad, criterion,
      drop.threshold, merge.threshold,
      max.iterations, initial.temp, cooling.rate, 
      temp.schedule, verbose)
    fit_set[[count]] <- sa_result
    
    if(verbose){
      cat("Best", criterion, "for Rate Class", current_rate_category, "after SA:", round(sa_result$model_summary[1,2], 3), "\n")
      cat("SA iterations:", sa_result$iterations, "\n")
      cat("SA acceptance rate:", round(sa_result$acceptance_rate, 3), "\n")
    }
    
    if(current_rate_category < max.rate.cat){
      current_rate_category <- current_rate_category + 1
      curr_index_mat <- getFullMat(replicate(current_rate_category, max_index_mat_cp, simplify = FALSE), 
        getStateMat(current_rate_category))
      curr_index_mat[curr_index_mat > 0] <- 1:sum(curr_index_mat > 0, na.rm = TRUE)
      curr_index_mat[curr_index_mat==0] <- NA
      max_index_mat <- curr_index_mat
      model_improved <- TRUE
    }else{
      model_improved <- FALSE
    }
  }
  
  if(verbose){
    cat("\nDone.\n")
  }
  
  all_models <- list()
  for(i in 1:length(fit_set)){
    all_models <- c(all_models, fit_set[[i]]$all_models)
  }
  all_models <- prune_redundant(all_models)
  class(all_models) <- "corhmm.dredge"
  
  if(!return.all){
    return(all_models)
  }else{
    return(list(all_models=all_models, sa_fits = fit_set))
  }
}

# SA within a single rate category
sa_within_rate_category <- function(phy, data, initial_fit, initial_index_mat, 
  max_index_mat, rate_category, root.p, pen.type, lambda, node.states, 
  fixed.nodes, ip, nstarts, n.cores, get.tip.states, lewis.asc.bias, collapse, 
  lower.bound, upper.bound, opts, grad, criterion, drop.threshold, 
  merge.threshold, max.iterations, initial.temp, cooling.rate, temp.schedule, 
  verbose, restart.strategy = "fixed_steps", # "fixed_steps", "energy_threshold", "random", "none"
  restart.interval = 20, restart.threshold = 1.5, restart.probability = 0.05,
  restart.temp.reset = FALSE) {
  
  # Initialize SA state
  current_fit <- initial_fit
  index_id <- paste0(c(current_fit$index.mat), collapse = "_")
  current_index_mat <- initial_index_mat
  current_score <- current_fit[[criterion]]
  
  # Initialize storage for all unique models
  all_models <- list() # all models refers to accpeted models
  every_model <- list() # this is every tested model
  all_index_mats <- list()
  every_score <- all_scores <- numeric()
  model_ids <- character()
  every_move <- character()
  accepted <- numeric()
  restarted <- numeric()
  
  # Store initial model
  initial_id <- paste0(c(initial_index_mat), collapse = "_")
  every_model[[1]] <- all_models[[1]] <- initial_fit
  all_index_mats[[1]] <- initial_index_mat
  every_score[1] <- all_scores[1] <- current_score
  model_ids[1] <- initial_id
  every_move[1] <- "none"
  accepted[1] <- 1
  
  best_fit <- current_fit
  best_index_mat <- current_index_mat
  best_score <- current_score
  
  current_temp <- initial.temp
  accepted_moves <- 0
  total_moves <- 0
  restart_count <- 0
  steps_since_restart <- 0
  steps_since_accept <- 0
  steps_since_best <- 0
  
  if(verbose){
    cat("Starting SA with", criterion, "=", round(current_score, 3), "\n")
    cat("Restart strategy:", restart.strategy, "\n")
  }
  
  for(iteration in 1:max.iterations) {
    steps_since_restart <- steps_since_restart + 1
    steps_since_accept <- steps_since_accept + 1
    steps_since_best <- steps_since_best + 1
    
    # Check restart conditions
    should_restart <- FALSE
    restart_reason <- ""
    
    if(restart.strategy == "fixed_steps" && steps_since_best >= restart.interval 
      && steps_since_restart >= restart.interval
      && steps_since_accept >= restart.interval) {
      should_restart <- TRUE
      restart_reason <- "fixed_steps"
    } else if(restart.strategy == "energy_threshold" && current_score > best_score * restart.threshold) {
      should_restart <- TRUE
      restart_reason <- "energy_threshold"
    } else if(restart.strategy == "random" && runif(1) < restart.probability) {
      should_restart <- TRUE
      restart_reason <- "random"
    }
    
    # Perform restart if needed
    if(should_restart && restart.strategy != "none") {
      every_move[iteration+1] <- "restart"
      available_indices <- which(all_scores != best_score)
      if(length(available_indices) > 0) {
        restart_idx <- sample(available_indices, 1)
        current_fit <- all_models[[restart_idx]]
        current_index_mat <- all_index_mats[[restart_idx]]
        current_score <- all_scores[restart_idx]
        restart_source <- paste0("accepted_model_", restart_idx)
      }else{
        current_fit <- best_fit
        current_index_mat <- best_index_mat
        current_score <- best_score
      }
      if(restart.temp.reset) {
        current_temp <- initial.temp
      }
      restart_count <- restart_count + 1
      steps_since_restart <- 0
      if(verbose) {
        cat("RESTART", restart_count, "at iteration", iteration, 
          "(", restart_reason, ") - Back to a previous model", criterion, ":", round(current_score, 3), "\n")
      }
      next  # Skip to next iteration after restart
    }
    
    # Update temperature
    if(temp.schedule == "exponential") {
      if(restart.temp.reset) {
        # Use steps since restart for temperature calculation
        current_temp <- initial.temp * (cooling.rate ^ steps_since_restart)
      } else {
        # Use total iterations for temperature calculation
        current_temp <- initial.temp * (cooling.rate ^ iteration)
      }
    } else if(temp.schedule == "linear") {
      if(restart.temp.reset) {
        # Linear cooling based on steps since restart
        remaining_steps <- max.iterations - iteration
        restart_window <- min(restart.interval, remaining_steps)
        current_temp <- initial.temp * (1 - steps_since_restart/restart_window)
      } else {
        current_temp <- initial.temp * (1 - iteration/max.iterations)
      }
    }
    
    # Stop if temperature too low
    if(current_temp < 0.001) break
    
    # Propose a move (stochastic drop or merge)
    move_result <- propose_sa_move_within_rate_cat(current_fit, drop.threshold, merge.threshold, max_index_mat)
    
    if(is.null(move_result$new_index_mat)) {
      next  # No valid move available
    }
    
    if(rate_category > 1 && test_hmm(move_result$new_index_mat, rate_category)){
      next # invalid HMM
    }
    
    
    move_id <- paste0(c(move_result$new_index_mat), collapse = "_")
    if(move_id %in% index_id){
      next
    } else {
      index_id <- c(index_id, move_id)
    }
    
    total_moves <- total_moves + 1
    
    # Fit the proposed model
    proposed_fit <- try(corHMMDredgeBase(phy=phy, data=data, rate.cat=rate_category, 
      root.p=root.p, pen.type = pen.type, lambda = lambda, 
      rate.mat=move_result$new_index_mat, node.states = node.states, 
      fixed.nodes=fixed.nodes, ip=ip, nstarts=nstarts, 
      n.cores=n.cores, get.tip.states = get.tip.states, 
      lewis.asc.bias = lewis.asc.bias, collapse = collapse, 
      lower.bound = lower.bound, upper.bound = upper.bound, 
      opts=opts, p=NULL, grad=grad))
    
    if(inherits(proposed_fit, "try-error")) {
      next  # Skip failed fits
    }
    if(proposed_fit$loglik == -1e+06){
      next  # Skip failed fits
    }
    
    proposed_score <- proposed_fit[[criterion]]
    every_model[[iteration+1]] <- proposed_fit
    every_score[iteration+1] <- proposed_score
    every_move[iteration+1] <- move_result$move_type
    
    # Accept/reject decision (lower score is better for AIC/BIC)
    delta <- proposed_score - current_score
    accept_prob <- if(delta <= 0) {
      1.0  # Always accept improvement
    } else {
      exp(-delta / current_temp)  # Metropolis criterion
    }
    
    if(runif(1) < accept_prob) {
      accepted[iteration+1] <- TRUE
      # Accept the move
      current_fit <- proposed_fit
      current_index_mat <- move_result$new_index_mat
      current_score <- proposed_score
      accepted_moves <- accepted_moves + 1
      steps_since_accept <- 0
      
      # Store this accepted model
      model_count <- length(all_models) + 1
      all_models[[model_count]] <- current_fit
      all_index_mats[[model_count]] <- current_index_mat
      all_scores[model_count] <- current_score
      model_ids[model_count] <- move_id
      
      # Update best if this is the best so far
      if(current_score < best_score) {
        best_fit <- current_fit
        best_index_mat <- current_index_mat
        best_score <- current_score
        steps_since_best <- 0
      }
      if(verbose) {
        cat("Iter", iteration, "- New", paste0(criterion, ":"), round(current_score, 3),
          "- Best", paste0(criterion, ":"), round(best_score, 3), "\n",
          "Move:", move_result$move_type, "- Temp:", round(current_temp, 4), 
          "- Steps since restart:", steps_since_restart, "\n",
          "Index Matrix:\n")
        print(move_result$new_index_mat)
        cat("\n")
      }
    }
  }
  
  acceptance_rate <- if(total_moves > 0) accepted_moves / total_moves else 0
  
  # Create summary data frame of all models
  model_summary <- data.frame(
    model_id = model_ids,
    score = all_scores,
    stringsAsFactors = FALSE
  )
  
  # Sort by score (best first)
  model_summary <- model_summary[order(model_summary$score), ]
  
  return(list(
    best_fit = best_fit,
    best_index_mat = best_index_mat,
    best_score = best_score,
    all_models = all_models,
    all_index_mats = all_index_mats,
    all_scores = all_scores,
    model_ids = model_ids,
    model_summary = model_summary,
    iterations = iteration,
    acceptance_rate = acceptance_rate,
    total_unique_models = length(all_models),
    restart_count = restart_count,
    final_steps_since_best = steps_since_best,
    every_model = every_model,
    every_score = every_score,
    every_move = every_move,
    accepted = accepted
  ))
}

# Propose stochastic moves within rate category
propose_sa_move_within_rate_cat <- function(current_fit, drop.threshold, merge.threshold, max_index_mat) {
  
  # Randomly choose between drop and merge
  move_type <- sample(c("drop", "merge", "free"), 1, prob = c(1/3, 1/3, 1/3))
  
  if(move_type == "drop") {
    new_index_mat <- propose_stochastic_drop(current_fit, drop.threshold)
  }
  if(move_type == "merge") {
    new_index_mat <- propose_stochastic_merge(current_fit, merge.threshold)
  } 
  if(move_type == "free") {
    new_index_mat <- propose_stochastic_free(current_fit, max_index_mat)
  } 

  return(list(
    move_type = move_type,
    new_index_mat = new_index_mat
  ))
}

# Stochastic parameter dropping
propose_stochastic_drop <- function(current_fit, drop.threshold) {
  Q <- current_fit$solution
  Q[is.na(Q)] <- Inf
  
  # Find parameters below threshold
  small_pars <- which(Q < drop.threshold)
  
  # If no small parameters, consider dropping from smallest quartile
  if(length(small_pars) == 0) {
    all_pars <- Q[!is.infinite(Q)]
    if(length(all_pars) == 0) return(NULL)
    
    threshold_25 <- quantile(all_pars, 0.25)
    candidate_pars <- which(Q <= threshold_25)
    if(length(candidate_pars) > 0) {
      small_pars <- candidate_pars
    }
  }
  
  if(length(small_pars) == 0) return(NULL)
  
  # Stochastically select parameters to drop
  # Probability inversely related to parameter size
  par_values <- Q[small_pars]
  drop_probs <- 1 / (par_values + 1e-10)
  drop_probs <- drop_probs / sum(drop_probs)
  
  # Select parameter(s) to drop
  n_to_drop <- sample(1:min(3, length(small_pars)), 1)
  if(length(small_pars) == 1){
    to_drop <- small_pars
  }else{
    to_drop <- sample(small_pars, n_to_drop, prob = drop_probs)
  }
  new_index_mat <- current_fit$index.mat
  new_index_mat[to_drop] <- NA
  pars <- sort(unique(na.omit(as.vector(new_index_mat))))
  for(i in 1:length(pars)){
    new_index_mat[new_index_mat == pars[i]] <- i
  }
  return(new_index_mat)
}

# Stochastic parameter merging
propose_stochastic_merge <- function(current_fit, merge.threshold) {
  if(current_fit$rate.cat > 1) {
    # Multi-rate category case
    current_pars <- MatrixToPars(current_fit)
    rate_classes <- paste("R", 1:current_fit$rate.cat, sep = "")
    par_list <- vector("list", current_fit$rate.cat+1)
    index_list <- vector("list", current_fit$rate.cat+1)
    
    for(i in seq(current_fit$rate.cat)){
      search_string <- paste0(rate_classes[i], " .* -> ", rate_classes[i])
      index_list[[i]] <- grep(search_string, names(current_pars))
      par_list[[i]] <- current_pars[index_list[[i]]]
    }
    index_list[[current_fit$rate.cat+1]] <- (1:length(current_pars))[-unlist(index_list)]
    par_list[[current_fit$rate.cat+1]] <- current_pars[index_list[[current_fit$rate.cat+1]]]
    
    # Check if any rate class has enough parameters
    valid_classes <- which(sapply(par_list, length) > 1)
    if(length(valid_classes) == 0) return(NULL)
    
    # Stochastically select which class to merge within
    selected_class <- sample(valid_classes, 1)
    selected_pars <- par_list[[selected_class]]
    selected_indices <- index_list[[selected_class]]
    
    # Stochastic merge within selected class
    merger_indices <- stochastic_merge_pars(selected_pars, merge.threshold)
    if(is.null(merger_indices)) return(NULL)
    
    focal_merger <- selected_indices[merger_indices]
    new_index_mat <- equateStateMatPars(current_fit$index.mat, focal_merger)
    
  } else {
    # Single rate category case
    current_pars <- MatrixToPars(current_fit)
    focal_merger <- stochastic_merge_pars(current_pars, merge.threshold)
    if(is.null(focal_merger)) return(NULL)
    
    new_index_mat <- equateStateMatPars(current_fit$index.mat, focal_merger)
  }
  pars <- sort(unique(na.omit(as.vector(new_index_mat))))
  for(i in 1:length(pars)){
    new_index_mat[new_index_mat == pars[i]] <- i
  }
  return(new_index_mat)
}

propose_stochastic_free <- function(current_fit, max_index_mat) {
  duplicates <- !is.na(current_fit$index.mat) & duplicated(current_fit$index.mat, MARGIN = 0)
  dropped <- is.na(current_fit$index.mat) & !is.na(max_index_mat)
  if(sum(dropped | duplicates) == 0) return(NULL)
  n_free <- sample(1:min(3, sum(dropped | duplicates)), 1)
  focal_free <- sample(which(dropped | duplicates), n_free)
  new_index_mat <- current_fit$index.mat
  new_index_mat[focal_free] <- max(current_fit$index.mat, na.rm = TRUE)+1:n_free
  pars <- sort(unique(na.omit(as.vector(new_index_mat))))
  for(i in 1:length(pars)){
    new_index_mat[new_index_mat == pars[i]] <- i
  }
  return(new_index_mat)
}


# Stochastic version of merge_current_pars
stochastic_merge_pars <- function(current_pars, merge.threshold) {
  if(length(current_pars) <= 1) return(NULL)
  
  # Compute distance matrix
  dist_mat <- as.matrix(dist(current_pars))
  dist_mat[upper.tri(dist_mat, diag = TRUE)] <- Inf
  
  # Find pairs within merge threshold, or use closest pairs
  valid_pairs <- which(dist_mat <= merge.threshold, arr.ind = TRUE)
  if(nrow(valid_pairs) == 0) {
    # No pairs within threshold, consider closest pairs with some randomness
    min_dist <- min(dist_mat[dist_mat != Inf])
    # Allow up to 20% larger distance to add stochasticity
    tolerance <- min_dist * (1 + runif(1) * 0.2)
    valid_pairs <- which(dist_mat <= tolerance, arr.ind = TRUE)
  }
  
  if(nrow(valid_pairs) == 0) return(NULL)
  
  # Stochastically select pair (closer pairs more likely)
  pair_distances <- sapply(1:nrow(valid_pairs), function(i) {
    dist_mat[valid_pairs[i,1], valid_pairs[i,2]]
  })
  
  # Inverse probability (closer pairs more likely)
  merge_probs <- 1 / (pair_distances + 1e-10)
  merge_probs <- merge_probs / sum(merge_probs)
  
  selected_pair_idx <- sample(nrow(valid_pairs), 1, prob = merge_probs)
  focal_merger <- valid_pairs[selected_pair_idx, ]
  
  # Expand cluster as in original
  avg_par <- mean(current_pars[focal_merger])
  additional_mergers <- which(abs(current_pars - avg_par) < merge.threshold)
  
  if(length(additional_mergers) > length(focal_merger)){
    focal_merger <- additional_mergers
  }
  
  return(focal_merger)
}

prune_redundant <- function(model_list){
  model_table <- getModelTable(model_list)
  model_table$rounded_lnLik <- round(model_table$lnLik, 6)
  model_table$n_rates <- unlist(lapply(model_list, function(x) sum(!is.na(x$solution))))
  keep_indices <- integer(0)
  for (ll in unique(model_table$rounded_lnLik)) {
    group_idx <- which(model_table$rounded_lnLik == ll)
    if (length(group_idx) == 1) {
      keep_indices <- c(keep_indices, group_idx)
    } else {
      group <- model_table[group_idx, ]
      min_pars <- min(group$np)
      best_by_pars <- group[group$np == min_pars, ]
      if (nrow(best_by_pars) == 1) {
        keep_idx <- rownames(best_by_pars)
      } else {
        min_rates <- min(best_by_pars$n_rates)
        best_by_both <- best_by_pars[best_by_pars$n_rates == min_rates, ]
        keep_idx <- rownames(best_by_both)[1]
      }
      keep_indices <- c(keep_indices, as.integer(keep_idx))
    }
  }
  model_table_unique <- model_table[keep_indices, ]
  pruned_model_list <- model_list[keep_indices]
  return(pruned_model_list)
}

get_best_info_criterion <- function(corhmm.obj.list, criterion, rate.cat){
  if(length(corhmm.obj.list) == 1){
    return(Inf)
  }
  info_criteria <- unlist(lapply(corhmm.obj.list, "[[", criterion))
  rate_cats <- unlist(lapply(corhmm.obj.list, "[[", "rate.cat"))
  if(any(rate_cats == rate.cat)){
    best_info_criterion <- min(info_criteria[rate_cats == rate.cat])
  }else{
    return(Inf)
  }
  return(best_info_criterion)
}

drop_pars <- function(corhmm.obj, drop.threshold){
  Q <- corhmm.obj$solution
  Q[is.na(Q)] <- Inf
  to_drop <- unique(corhmm.obj$index.mat[which(Q < drop.threshold)])
  if(length(to_drop) == 0){
    return(NULL)
  }
  index_mat <- dropStateMatPars(corhmm.obj$index.mat, to_drop)
  return(index_mat)
}

merge_pars <- function(corhmm.obj, merge.threshold){
  index_mat <- corhmm.obj$index.mat
  if(corhmm.obj$rate.cat > 1){
    current_pars <- MatrixToPars(corhmm.obj)
    rate_classes <- paste("R", 1:corhmm.obj$rate.cat, sep = "")
    par_list <- vector("list", corhmm.obj$rate.cat+1)
    index_list <- vector("list", corhmm.obj$rate.cat+1)
    for(i in seq(corhmm.obj$rate.cat)){
      search_string <- paste0(rate_classes[i], " .* -> ", rate_classes[i])
      index_list[[i]] <- grep(search_string, names(current_pars))
      par_list[[i]] <- current_pars[index_list[[i]]]
    }
    index_list[[corhmm.obj$rate.cat+1]] <- (1:length(current_pars))[-unlist(index_list)]
    par_list[[corhmm.obj$rate.cat+1]] <- current_pars[index_list[[corhmm.obj$rate.cat+1]]]
    er_test <- all(unlist(lapply(index_list, function(x) length(x) <= 1)))
    if(er_test){
      return(NULL)
    }
    potential_mergers <- lapply(par_list, 
      function(x) merge_current_pars(x, merge.threshold))
    diffs <- mapply(function(x,y){abs(diff(y[x]))}, x=potential_mergers, y=par_list)
    focal_merger <- index_list[[which.min(diffs)]][potential_mergers[[which.min(diffs)]]]
    index_mat_merged <- equateStateMatPars(index_mat, focal_merger)
  }else{
    current_pars <- MatrixToPars(corhmm.obj)
    focal_merger <- merge_current_pars(current_pars, merge.threshold)
    if(is.null(focal_merger)){
      return(NULL)
    }
    index_mat_merged <- equateStateMatPars(index_mat, focal_merger)
  }
  return(index_mat_merged)
}

test_hmm <- function(index_mat, rate_cat){
  rate_cat_tests <- vector(length = rate_cat)
  rate_class_labels <- paste0("R", 1:rate_cat)
  for(i in 1:length(rate_cat_tests)){
    to_rc <- index_mat[,grep(rate_class_labels[i], colnames(index_mat))]
    from_rc <- index_mat[grep(rate_class_labels[i], rownames(index_mat)),]
    rate_cat_tests[i]<- !all(is.na(to_rc)) | !all(is.na(from_rc))
  }
  return(any(!rate_cat_tests))
}

test_validity_hmm <- function(corhmm_obj){
  index_mat <- corhmm_obj$index.mat
  rate_cat <- corhmm_obj$rate.cat
  rate_cat_tests <- vector(length = rate_cat)
  rate_class_labels <- paste0("R", 1:rate_cat)
  for(i in 1:length(rate_cat_tests)){
    to_rc <- index_mat[,grep(rate_class_labels[i], colnames(index_mat))]
    from_rc <- index_mat[grep(rate_class_labels[i], rownames(index_mat)),]
    rate_cat_tests[i]<- !all(is.na(to_rc)) | !all(is.na(from_rc))
  }
  return(rate_cat_tests)
}

merge_current_pars <- function(current_pars, merge.threshold){
  if(length(current_pars) <= 1){
    return(NULL)
  }
  dist_mat <- as.matrix(dist(current_pars))
  dist_mat[upper.tri(dist_mat)] <- 0
  focal_merger <- which(dist_mat == min(dist(current_pars)), arr.ind = TRUE)
  avg_par <- mean(current_pars[focal_merger])
  additional_mergers <- which(abs(current_pars - avg_par) < merge.threshold)
  if(length(additional_mergers) > length(focal_merger)){
    focal_merger <- additional_mergers
  }
  return(focal_merger)
}

# this is the function that does most of the heavy lifting
corHMMDredgeBase <- function(phy, data, rate.cat, root.p="maddfitz", pen.type = "l1", lambda = 1, rate.mat=NULL, node.states = "marginal", fixed.nodes=FALSE, ip=NULL, nstarts=0, n.cores=1, get.tip.states = FALSE,lewis.asc.bias = FALSE, collapse = FALSE, lower.bound = 1e-10, upper.bound = 100, opts=NULL, p=NULL, grad=FALSE){
  
  # Checks to make sure node.states is not NULL.  If it is, just returns a diagnostic message asking for value.
  if(is.null(node.states)){
    obj <- NULL
    obj$loglik <- NULL
    obj$diagnostic <- paste("No model for ancestral states selected.  Please pass one of the following to corHMM command for parameter \'node.states\': joint, marginal, scaled, or none.")
    return(obj)
  } else { # even if node.states is not NULL, need to make sure its one of the three valid options
    valid.models <- c("joint", "marginal", "scaled", "none")
    if(!any(valid.models == node.states)){
      obj <- NULL
      obj$loglik <- NULL
      obj$diagnostic <- paste("\'",node.states, "\' is not valid for ancestral state reconstruction method.  Please pass one of the following to corHMM command for parameter \'node.states\': joint, marginal, scaled, or none.",sep="")
      return(obj)
    }
    if(length(node.states) > 1){ # User did not enter a value, so just pick marginal.
      node.states <- "marginal"
      cat("No model selected for \'node.states\'. Will perform marginal ancestral state estimation.\n")
    }
  }
  
  if(fixed.nodes == FALSE){
    if(!is.null(phy$node.label)){
      phy$node.label <- NULL
      cat("You specified \'fixed.nodes=FALSE\' but included a phy object with node labels. These node labels have been removed.\n")
    }
  }
  
  #Ensures that weird root state probabilities that do not sum to 1 are input:
  if(!is.null(root.p)){
    if(!is.character(root.p)){
      root.p <- root.p/sum(root.p)
    }
  }
  
  if(pen.type == "unreg"){
    lambda <- 0
  }
  
  # rescale phy to height of one
  phy_original <- phy
  H <- max(node.depth.edgelength(phy))
  phy$edge.length <- phy$edge.length/H
  upper.bound <- upper.bound * H
  lower.bound <- lower.bound / H
  
  input.data <- data
  
  nCol <- dim(data)[2]
  
  CorData <- corProcessData(data, collapse = collapse)
  data.legend <- data <- CorData$corData
  nObs <- length(CorData$ObservedTraits)

    # Checks to make sure phy & data have same taxa. Fixes conflicts (see match.tree.data function).
  matching <- match.tree.data(phy,data)
  data <- matching$data
  phy <- matching$phy
  
  # Will not perform reconstructions on invariant characters (unless rate params have been given!)
  if(nlevels(as.factor(data[,1])) <= 1 & !is.null(p)){
    obj <- NULL
    obj$loglik <- NULL
    obj$diagnostic <- paste("Character is invariant. Analysis stopped.",sep="")
    return(obj)
  } else {
    # Still need to make sure second level isnt just an ambiguity
    lvls <- as.factor(data[,1])
    if(nlevels(as.factor(data[,1])) == 2 && length(which(lvls == "?"))){
      obj <- NULL
      obj$loglik <- NULL
      obj$diagnostic <- paste("Character is invariant. Analysis stopped.",sep="")
      return(obj)
    }
  }
  
  if(any(phy$edge.length<=.Machine$double.eps)){
    warning(paste0("Branch lengths of 0 detected. Adding ", sqrt(.Machine$double.eps)), immediate. = TRUE)
    #   phy$edge.length[phy$edge.length<=1e-5] <- 1e-5
    phy$edge.length <- phy$edge.length + sqrt(.Machine$double.eps) 
  }
  #Creates the data structure and orders the rows to match the tree.
  data.sort <- data.frame(data[,2], data[,2],row.names=data[,1])
  data.sort <- data.sort[phy$tip.label,]
  
  counts <- table(data.sort[,1])
  levels <- levels(as.factor(data.sort[,1]))
  cols <- as.factor(data.sort[,1])

  #Some initial values for use later
  k=2
  if(upper.bound < lower.bound){
    cat("Your upper bound is smaller than your lower bound.\n")
  }
  lb <- log(lower.bound)
  ub <- log(upper.bound)
  order.test <- FALSE
  
  obj <- NULL
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  root.p <- root.p
  nstarts <- nstarts
  ip <- ip
  model = "ARD"
  
  model.set.final <- rate.cat.set.corHMM.JDB(phy=phy,data=input.data,rate.cat=rate.cat,ntraits=nObs,model=model,rate.mat=rate.mat, collapse=collapse)
  # adjusting the matrix for corhmm dredge which allows for independent rate classes
  if(is.null(rate.mat)){
    model.set.final$index.matrix[!is.na(model.set.final$index.matrix)] <- 1:length(model.set.final$index.matrix[!is.na(model.set.final$index.matrix)])
    model.set.final$rate <- model.set.final$index.matrix
    model.set.final$rate[is.na(model.set.final$rate)] <- max(model.set.final$rate, na.rm = TRUE) + 1
    model.set.final$np <- max(model.set.final$index.matrix, na.rm = TRUE)
  }
  phy <- reorder(phy, "pruningwise")
  
  if(collapse){
    StateNames <- gsub("_", "|", CorData$ObservedTraits)
  }else{
    StateNames <- gsub("_", "|", CorData$PossibleTraits)
  }  
  print_counts <- rep(0, length(StateNames))
  if(length(grep("&", names(counts))) > 0){
    counts <- counts[-grep("&", names(counts))]
  }
  print_counts[as.numeric(names(counts))] <- counts
  # cat("State distribution in data:\n")
  # cat("States:",StateNames,"\n",sep="\t")
  # cat("Counts:",print_counts,"\n",sep="\t")
  
  lower = rep(lb, model.set.final$np)
  upper = rep(ub, model.set.final$np)
  
  if(is.null(opts)){
    if(grad){
      opts <- list("algorithm"="NLOPT_LD_MMA", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
    }else{
      opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
    }
  }
  if(!is.null(p)){
    cat("Calculating likelihood from a set of fixed parameters", "\n")
    out<-NULL
    p <- p*H
    est.pars<-log(p)
    out$objective<-dev.corhmm.dredge(est.pars,phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p, rate.cat = rate.cat, order.test = order.test, lewis.asc.bias = lewis.asc.bias, pen.type = pen.type, lambda = lambda)
    loglik <- -out$objective
    est.pars <- exp(est.pars)
  }else{
    if(is.null(ip)){
      #If a user-specified starting value(s) is not supplied this begins loop through a set of randomly chosen starting values:
      #Sets parameter settings for random restarts by taking the parsimony score and dividing
      #by the total length of the tree
      # cat("Beginning thorough optimization search -- performing", nstarts, "random restarts", "\n")
      taxa.missing.data.drop <- which(is.na(data.sort[,1]))
      if(length(taxa.missing.data.drop) != 0){
        tip.labs <- names(taxa.missing.data.drop)
        dat <- as.matrix(data.sort)
        dat.red <- dat[-taxa.missing.data.drop,]
        phy.red <- drop.tip(phy, taxa.missing.data.drop)
        dat.red <- phyDat(dat.red,type="USER", levels=levels)
        phy.tmp <- multi2di(phy.red)
        par.score <- parsimony(phy.tmp, dat.red, method="fitch")/2
      }else{
        dat <- as.matrix(data.sort)
        dat <- phyDat(dat,type="USER", levels=levels)
        phy.tmp <- multi2di(phy)
        par.score <- parsimony(phy.tmp, dat, method="fitch")/2
      }
      tl <- sum(phy$edge.length)
      mean.change = par.score/tl
      
      random.restart<-function(nstarts){
        tmp = matrix(,1,ncol=(1+model.set.final$np))
        if(mean.change==0){
          starts=rep(0.01+exp(lb), model.set.final$np)
        }else{
          starts<-sort(rexp(model.set.final$np, 1/mean.change), decreasing = TRUE)
        }
        starts[starts < exp(lb)] = exp(lb)
        starts[starts > exp(ub)] = exp(lb)
        out = nloptr(x0=log(starts), eval_f=dev.corhmm.dredge, lb=lower, ub=upper, opts=opts, phy=phy, liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p, rate.cat = rate.cat, order.test = order.test, lewis.asc.bias = lewis.asc.bias, pen.type = pen.type, lambda = lambda, grad=grad)
        tmp[,1] = out$objective
        tmp[,2:(model.set.final$np+1)] = out$solution
        tmp
      }
      # this is the first model pass
      if(n.cores > 1){
        restart.set<-mclapply(1:nstarts, random.restart, mc.cores=n.cores)
      }else{
        restart.set<-lapply(1:nstarts, random.restart)
      }
      #Finds the best fit within the restart.set list
      best.fit<-which.min(unlist(lapply(restart.set, function(x) x[1])))
      #Generates an object to store results from restart algorithm:
      out<-NULL
      out$objective=unlist(restart.set[[best.fit]][,1])
      out$solution=unlist(restart.set[[best.fit]][,2:(model.set.final$np+1)])
      loglik <- -out$objective
      est.pars <- exp(out$solution)
    }else{
      # the user has specified initial params
      # cat("Beginning subplex optimization routine -- Starting value(s):", ip, "\n")
      ip=ip
      out = nloptr(x0=rep(log(ip), length.out = model.set.final$np), eval_f=dev.corhmm.dredge, lb=lower, ub=upper, opts=opts, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p, rate.cat = rate.cat, order.test = order.test, lewis.asc.bias = lewis.asc.bias, pen.type = pen.type, lambda = lambda, grad=grad)
      loglik <- -out$objective
      est.pars <- exp(out$solution)
    }
  }
  
  #Starts the ancestral state reconstructions:
  if(node.states != "none") {
    # cat("Finished. Inferring ancestral states using", node.states, "reconstruction.","\n")
  }
  TIPS <- 1:nb.tip
  if (node.states == "marginal" || node.states == "scaled"){
    lik.anc <- ancRECON(phy, input.data, est.pars, rate.cat, rate.mat=model.set.final$rate, method=node.states, ntraits=NULL, root.p=root.p, model = model, get.tip.states = get.tip.states, collapse = collapse)
    pr<-apply(lik.anc$lik.anc.states,1,which.max)
    phy$node.label <- pr
    tip.states <- lik.anc$lik.tip.states
    row.names(tip.states) <- phy$tip.label
  }
  if (node.states == "joint"){
    lik.anc <- ancRECON(phy, input.data, est.pars, rate.cat, rate.mat=model.set.final$rate, method=node.states, ntraits=NULL, root.p=root.p, model = model, get.tip.states = get.tip.states, collapse = collapse)
    phy$node.label <- lik.anc$lik.anc.states
    tip.states <- lik.anc$lik.tip.states
  }
  if (node.states == "none") {
    lik.anc <- list(lik.tip.states=NA, lik.anc.states=NA, info.anc.states=NA)
    phy$node.label <- NA
    tip.states <- NA
  }
  
  # finalize the output
  solution <- matrix(est.pars[model.set.final$index.matrix], dim(model.set.final$index.matrix))
  if(collapse){
    StateNames <- rep(gsub("_", "|", CorData$ObservedTraits), rate.cat)
    RCNames <- rep(paste("R", 1:rate.cat, sep = ""), each = length(CorData$ObservedTraits))
  }else{
    StateNames <- rep(gsub("_", "|", CorData$PossibleTraits), rate.cat)
    RCNames <- rep(paste("R", 1:rate.cat, sep = ""), each = length(CorData$PossibleTraits))
  }
  if(rate.cat > 1){
    StateNames <- paste(RCNames, StateNames)
  }
  np <- model.set.final$np
  index.matrix <- model.set.final$index.matrix
  rownames(solution) <- colnames(solution) <- StateNames
  rownames(index.matrix) <- colnames(index.matrix) <- StateNames
  # rescale phylogeny
  solution <- solution/H
  solution[solution < lower.bound] <- lower.bound
  
  AIC <- -2*loglik+2*np
  AICc <- -2*loglik+(2*np*(nb.tip/(nb.tip-np-1)))
  
  if (is.character(node.states)) {
    if (node.states == "marginal" || node.states == "scaled"){
      colnames(lik.anc$lik.anc.states) <- StateNames
    }
  }
  
  if(loglik == -1e+06){
    warning("corHMM may have failed to optimize correctly, consider checking inputs and running again.", immediate. = TRUE)
  }
  
  obj = list(loglik = loglik,
             AIC = AIC,
             AICc = AICc,
             rate.cat=rate.cat,
             solution=solution,
             index.mat=index.matrix,
             data=input.data,
             data.legend = data.legend,
             phy=phy_original,
             states=lik.anc$lik.anc.states,
             tip.states=tip.states,
             states.info = lik.anc$info.anc.states,
             iterations=out$iterations,
             collapse=collapse,
             root.p=root.p,
             pen.type=pen.type,
             lambda=lambda)
  class(obj)<-"corhmm"
  return(obj)
}


### The function used to optimize parameters:

dev.corhmm.dredge <- function(p,phy,liks,Q,rate,root.p,rate.cat,order.test,lewis.asc.bias,pen.type="l1",lambda=1,grad=FALSE){
  p = exp(p)
  cp_root.p <- root.p
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  TIPS <- 1:nb.tip
  comp <- numeric(nb.tip + nb.node)
  #Obtain an object of all the unique ancestors
  anc <- unique(phy$edge[,1])
  k.rates <- dim(Q)[2] / 2
  if (any(is.nan(p)) || any(is.infinite(p))) return(1000000)
  
  Q[] <- c(p, 0)[rate]
  diag(Q) <- -rowSums(Q)
  pen_score <- get_penalty_score(Q, p, pen.type, rate, rate.cat)
  # # if the q matrix has columns not estimated, remove them
  # row2rm <- apply(rate, 1, function(x) all(x == max(rate)))
  # col2rm <- apply(rate, 2, function(x) all(x == max(rate)))
  # Q.root <- Q[!row2rm | !col2rm, !row2rm | !col2rm]
  if(is.character(root.p)){
    if(root.p == "yang"){
      root.test <- Null(Q)
      if(dim(root.test)[2]>1){
        return(1000000)
      }
    }      
  }
  
  if(order.test == TRUE){
    # ensure that the rate classes have mean rates in a consistent order (A > B > C > n)
    StateOrderMat <- matrix(1, (dim(Q)/rate.cat)[1], (dim(Q)/rate.cat)[2])
    RateClassOrderMat <- matrix(0, rate.cat, rate.cat)
    diag(RateClassOrderMat) <- 1:rate.cat
    OrderMat <- RateClassOrderMat %x% StateOrderMat
    Rate01 <- vector("numeric", rate.cat)
    for(i in 1:rate.cat){
      tmp <- Q[OrderMat == i]
      Rate01[i] <- tmp[tmp>=0][1]
    }
    OrderTest <- all.equal(Rate01, sort(Rate01, decreasing = TRUE))
    if(OrderTest != TRUE){
      return(1000000)
    }
  }
  
  for (i in seq(from = 1, length.out = nb.node)) {
    #the ancestral node at row i is called focal
    focal <- anc[i]
    #Get descendant information of focal
    desRows <- which(phy$edge[,1]==focal)
    desNodes <- phy$edge[desRows,2]
    v <- 1
    #Loops through all descendants of focal (how we deal with polytomies):
    for (desIndex in sequence(length(desRows))){
      v <- v*expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks[desNodes[desIndex],]
    }
    
    ##Allows for fixed nodes based on user input tree.
    if(!is.null(phy$node.label)){
      if(!is.na(phy$node.label[focal - nb.tip])){
        fixer.tmp = numeric(dim(Q)[2]/rate.cat)
        fixer.tmp[phy$node.label[focal - nb.tip]] = 1
        fixer = rep(fixer.tmp, rate.cat)
        v <- v * fixer
      }
    }
    
    #Sum the likelihoods:
    comp[focal] <- sum(v)
    #Divide each likelihood by the sum to obtain probabilities:
    liks[focal, ] <- v/comp[focal]
  }
  
  #Specifies the root:
  root <- nb.tip + 1L
  #If any of the logs have NAs restart search:
  if (is.na(sum(log(comp[-TIPS])))){return(1000000)}
  equil.root <- NULL
  
  for(i in 1:ncol(Q)){
    posrows <- which(Q[,i] >= 0)
    rowsum <- sum(Q[posrows,i])
    poscols <- which(Q[i,] >= 0)
    colsum <- sum(Q[i,poscols])
    equil.root <- c(equil.root,rowsum/(rowsum+colsum))
  }
  if (is.null(root.p)){
    flat.root = equil.root
    k.rates <- 1/length(which(!is.na(equil.root)))
    flat.root[!is.na(flat.root)] = k.rates
    flat.root[is.na(flat.root)] = 0
    loglik<- -(sum(log(comp[-TIPS])) + log(sum(flat.root * liks[root,])))
  }
  if(is.character(root.p)){
    # root.p==yang will fix root probabilities based on the inferred rates: q10/(q01+q10)
    if(root.p == "yang"){
      root.p <- Null(Q)
      root.p <- c(root.p/sum(root.p))
      loglik <- -(sum(log(comp[-TIPS])) + log(sum(root.p * liks[root,])))
      if(is.infinite(loglik)){
        return(1000000)
      }
    }else{
      # root.p==maddfitz will fix root probabilities according to FitzJohn et al 2009 Eq. 10:
      root.p = liks[root,] / sum(liks[root,])
      loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,])))))
    }
  }else{
    if(is.numeric(root.p[1])){
      loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,])))))
      if(is.infinite(loglik)){
        return(1000000)
      }
    }
  }
  # root.p!==NULL will fix root probabilities based on user supplied vector:
  if(lewis.asc.bias == TRUE){
    p <- log(p)
    dummy.liks.vec <- getLewisLikelihood(p = p, phy = phy, liks = liks, Q = Q, rate = rate, root.p = cp_root.p, rate.cat = rate.cat)
    loglik <- loglik - log(sum(root.p * (1 - exp(dummy.liks.vec))))
  }
  if(grad){
    epsilon <- 1e-6
    base_log_likelihood <- -(loglik + (pen_score * lambda))
    gradient <- numeric(length(p))
    for (i in seq_along(p)) {
      p_plus <- log(p)
      p_plus[i] <- p_plus[i] + epsilon
      log_likelihood_plus <- -dev.corhmm.dredge(p_plus,phy,liks,Q,rate,root.p,rate.cat,order.test,lewis.asc.bias,pen.type,lambda,grad=FALSE)
      gradient[i] <- (log_likelihood_plus - base_log_likelihood) / epsilon
    }
    cat("\rlnLik:", paste0(base_log_likelihood, "..."))
    # print(-gradient)
    return(list(
      objective = -base_log_likelihood,
      gradient = -gradient))
  }else{
    return(loglik + (pen_score * lambda))
  }
}


### The function used to calculate the penalty:

get_penalty_score <- function(Q, p, pen.type, index.mat, rate.cat){
  if(rate.cat == 1){
    if(pen.type == "l1"){
      pen <- mean(-diag(Q))
    }
    if(pen.type == "l2"){
      pen <- mean(diag(Q)^2)
    }
    if(pen.type == "er"){
      # pen <- sd(-diag(Q))
      if(length(p) < 1){
        pen <- mean(dist(p))
      }else{
        pen <- 0
      }
    }
  }else{
    rate_class_names <- paste0("R", 1:rate.cat)
    pen_by_rc <- numeric(rate.cat)
    diag(Q) <- 0
    for(i in seq_along(rate_class_names)){
      rc_index <- grep(rate_class_names[i], colnames(index.mat))
      if(pen.type == "l1"){
        pen_by_rc[i] <- mean(Q[rc_index, rc_index])
      }
      if(pen.type == "l2"){
        pen_by_rc[i] <- mean(Q[rc_index, rc_index]^2)
      }
      if(pen.type == "er"){
        rates <- Q[rc_index, rc_index]
        rates <- rates[rates>0]
        # pen_by_rc[i] <- sd(rates)
        if(length(rates) < 1){
          pen_by_rc[i] <- mean(dist(rates))
        }else{
          pen_by_rc[i] <- 0
        }
      }
    }
    pen <- sum(pen_by_rc)
  }
  return(pen)
}

### Functions for leave one out cross validation:

get_per_tip_faith_PD <- function(phy, fold_vec){
  fold_hist <- setNames(rep(NA, length(unique(fold_vec))), sort(unique(fold_vec)))
  for(i in 1:length(sort(unique(fold_vec)))){
    focal_fold <- sort(unique(fold_vec))[i]
    focal_tips <- names(fold_vec)[fold_vec == focal_fold]
    focal_phy <- keep.tip(phy, focal_tips)
    fold_hist[i] <- sum(focal_phy$edge.length)/length(focal_tips)
  }
  return(fold_hist)
}

# other weight strategies to test
# split_data_k_other <- function(phy, k=5){
#   weights <- rowSums(vcv.phylo(phy))
#   fold_vec <- setNames(rep(NA, length(phy$tip.label)), phy$tip.label)
#   for(i in 1:(length(fold_vec)-1)){
#     curr_k <- i %% k
#     focal <- sample(names(weights), 1, replace = FALSE, prob = weights)
#     fold_vec[focal] <- curr_k
#     phy <- drop.tip(phy, focal)
#     weights <- rowSums(vcv.phylo(phy))
#   }
#   focal <- sample(names(weights), 1, replace = FALSE, prob = weights)
#   curr_k <- length(weights) %% k
#   fold_vec[focal] <- curr_k
#   return(fold_vec)
# }
# 
# split_k_fold_even <- function(phy, k=5){
#   weights <- setNames(rep(1, length(phy$tip.label)), phy$tip.label)
#   fold_vec <- setNames(rep(NA, length(phy$tip.label)), phy$tip.label)
#   for(i in 1:(length(fold_vec)-1)){
#     curr_k <- i %% k
#     focal <- sample(names(weights), 1, replace = FALSE, prob = weights)
#     fold_vec[focal] <- curr_k
#     weights <- weights[!names(weights) %in% focal]
#   }
#   focal <- sample(names(weights), 1, replace = FALSE, prob = weights)
#   curr_k <- length(weights) %% k
#   fold_vec[focal] <- curr_k
#   return(fold_vec)
# }

get_weights <- function(phy){
  # based on rholf 2001
  C <- vcv.phylo(phy)
  C_inv <- solve(C)
  I <- matrix(0, dim(C)[1], dim(C)[2])
  diag(I) <- 1
  weights <- (t(I) %*% C_inv) %*% matrix(1, dim(C)[1], 1) 
  weights <- c(weights/sum(weights))
  weights <- setNames(weights, colnames(C))
  return(weights)
}

split_data_k_folds <- function(phy, k=5){
  # weights <- rowSums(vcv.phylo(phy))
  weights <- get_weights(phy)
  fold_vec <- setNames(rep(NA, length(phy$tip.label)), phy$tip.label)
  for(i in 1:(length(fold_vec)-1)){
    curr_k <- i %% k
    focal <- sample(names(weights), 1, replace = FALSE, prob = weights)
    fold_vec[focal] <- curr_k
    weights <- weights[!names(weights) %in% focal]
  }
  focal <- sample(names(weights), 1, replace = FALSE, prob = weights)
  curr_k <- length(weights) %% k
  fold_vec[focal] <- curr_k
  return(fold_vec)
}

# Function to perform k-fold cross-validation
kFoldCrossValidation <- function(corhmm_obj, k, lambdas=NULL, return_model=TRUE, save_model_dir=NULL, model_name=NULL) {
  scores <- numeric(k)  # Create an empty vector to store scores for each fold
  folds <- split_data_k_folds(corhmm_obj$phy, k)
  ip <- MatrixToPars(corhmm_obj)
  if(is.null(lambdas)){
    if(is.null(corhmm_obj$lambda)){
      lambdas <- 0
    }else{
      lambdas <- corhmm_obj$lambda
    }
  }
  total_model_list <- vector("list", length(lambdas))
  count <- 1
  for(lambda in lambdas){
    model_list <- vector("list", k)
    cat("Evaluating lambda =", lambda, "\n")  # Print the score for the fold
    for (i in 0:(k-1)) {
      # get the original data
      fold_data <- corhmm_obj$data
      # Split data into training and testing sets for the current fold
      fold_data[folds == i, 2:ncol(fold_data)] <- "?"
      # Train the model on training data
      model <- corHMMDredgeBase(phy = corhmm_obj$phy,
                                         data = fold_data,
                                         rate.cat = corhmm_obj$rate.cat, 
                                         pen.type = corhmm_obj$pen.type, 
                                         lambda = lambda, 
                                         rate.mat = corhmm_obj$index.mat, 
                                         node.states = "marginal", 
                                         fixed.nodes=FALSE, 
                                         root.p=corhmm_obj$root.p, 
                                         ip=ip, 
                                         nstarts=0, 
                                         n.cores=1, 
                                         get.tip.states = TRUE, 
                                         lewis.asc.bias = FALSE, 
                                         collapse = corhmm_obj$collapse, 
                                         lower.bound = 1e-10, 
                                         upper.bound = 100, 
                                         opts=NULL, 
                                         p=NULL,
                                         grad=FALSE)
      
      # Evaluate the model on testing data
      score <- evaluateModel(model, corhmm_obj)
      scores[i+1] <- score  # Store the score for this fold
      if(return_model){
        model_list[[i+1]] <- model # Store the model for this fold
      }else{
        model_list[[i+1]] <- NULL
      }
      if(!is.null(save_model_dir)){
        if(is.null(model_name)){
          model_name <- "corhmm.obj"
        }
        saveRDS(model, file = paste0(save_model_dir, "/", model_name, "_lambda", lambda, "_fold", folds, ".RDS"))
      }
      cat("Fold", i, "Score:", score, "\n")  # Print the score for the fold
    }
    averageScore <- mean(scores)  # Calculate the average score across all folds
    cat("Average Cross-Validation Score:", averageScore, "\n")
    total_model_list[[count]] <- list(models = model_list, scores = scores, averageScore = averageScore)
    count <- count + 1
  }
  names(total_model_list) <- lambdas
  class(total_model_list) <- c("corhmm.kfold")
  return(total_model_list)
}

# Function to evaluate the model on testing data
evaluateModel <- function(model, corhmm_obj){
  tip_liks <- get_tip_liks(corhmm_obj)
  scores <- numeric(dim(tip_liks)[1])
  for(i in 1:nrow(tip_liks)){
    scores[i] <- js_divergence(tip_liks[i,], model$tip.states[i,])
  }
  score <- mean(scores)
  return(score)
}

# Calculate Total Variation Distance
total_variation_distance <- function(P, Q) {
  sum(abs(P - Q)) / 2
}

# Calculate KL Divergence
kl_divergence <- function(P, Q) {
  if (any(P > 0 & Q == 0)) {
    return(Inf)  # To handle the case where Q(i) = 0 and P(i) > 0
  }
  return(sum(P * log(P / Q), na.rm = TRUE))
}

# Calculate Jensen-Shannon Divergence
# This is a symmetric and smoothed version of KL divergence. It is defined as the average of the KL divergences between each distribution and the average of both distributions
js_divergence <- function(P, Q) {
  M <- (P + Q) / 2
  return((kl_divergence(P, M) + kl_divergence(Q, M)) / 2)
}

get_tip_liks <- function(corhmm_obj, return_original=TRUE){
  model.set.final <- get_MSF_from_corhm_obj(corhmm_obj)
  pars <- MatrixToPars(corhmm_obj)
  phy <- corhmm_obj$phy
  phy$node.label <- NULL
  ntips <- Ntip(phy)
  nnodes <- Nnode(phy)
  liks <- model.set.final$liks
  rownames(liks) <- c(phy$tip.label, (ntips+1):(ntips+nnodes))
  tip_liks <- liks[1:ntips,]
  if(return_original){
    return(tip_liks)
  }
  for(i in 1:ntips){
    liks_row <- c()
    liks_copy <- liks
    for(j in 1:ncol(liks)){
      liks_copy[i,] <- 0
      liks_copy[i, j] <- 1
      lik_tmp <- dev.corhmm(log(pars),
                            phy=phy,
                            liks=liks_copy,
                            Q=model.set.final$Q,
                            rate=model.set.final$rate,
                            root.p=corhmm_obj$root.p, 
                            rate.cat = corhmm_obj$rate.cat, 
                            order.test = FALSE, 
                            lewis.asc.bias = FALSE)
      liks_row[j] <- -lik_tmp
    }
    best_probs <- max(liks_row)
    liks_tip_rescaled <- liks_row - best_probs
    tip_liks[i,] <- exp(liks_tip_rescaled) / sum(exp(liks_tip_rescaled))
  }
  return(tip_liks)
}

get_MSF_from_corhm_obj <- function(corhmm_obj){
  model.set.final <- rate.cat.set.corHMM.JDB(
    phy=corhmm_obj$phy,data=corhmm_obj$data,
    rate.cat=corhmm_obj$rate.cat,
    ntraits=dim(corhmm_obj$index.mat)[1]/corhmm_obj$rate.cat,
    model=NULL,
    rate.mat=corhmm_obj$index.mat, 
    collapse=corhmm_obj$collapse)
  return(model.set.final)
}

getCVTable <- function(x){
  score_table <- do.call(rbind, lapply(x, "[[", "scores"))
  colnames(score_table) <- paste0("Fold:", 0:(dim(score_table)[2]-1))
  rownames(score_table) <- paste0("Lambda:", rownames(score_table))
  score_table <- t(score_table)
  avg_scores <- colMeans(score_table)
  # cat("\nScores per fold:\n")
  # print(score_table)
  # cat("\n")
  # cat("Average Scores:\n")
  # print(avg_scores)
  # cat("\n")
  return(list(score_table=score_table, avg_scores=avg_scores))
}

plotDredgeTrace <- function(dredge_fits,
  break_size = 5,
  palette = c("drop" = "#A23B72",
    "merge" = "#2E86AB",
    "free" = "#F18F01",
    "restart" = "#7209B7",
    "none" = "grey60"),
  legend = TRUE,
  legend.pos = "topright",
  ...) {
  
  # 1. Input Validation
  if (!is.list(dredge_fits) || !("sa_fits" %in% names(dredge_fits))) {
    stop("Input must be a corHMMDredge object created with return.all = TRUE.")
  }
  
  sa_fits <- dredge_fits$sa_fits
  num_fits <- length(sa_fits)
  if (num_fits == 0) {
    message("No simulated annealing fits found in the object to plot.")
    return(invisible(NULL))
  }
  
  # 2. Data Extraction and Combination
  combined_scores <- list()
  combined_moves <- list()
  combined_accepted <- list()
  fit_lengths <- numeric(num_fits)
  rate_cat_labels <- sapply(sa_fits, function(x) x$best_fit$rate.cat)
  
  for (i in 1:num_fits) {
    fit_data <- sa_fits[[i]]
    
    # Filter out NAs which can occur if the loop ends prematurely
    valid_indices <- !is.na(fit_data$every_move)
    
    scores <- fit_data$every_score[valid_indices]
    moves <- fit_data$every_move[valid_indices]
    accepted <- fit_data$accepted[valid_indices]
    
    # Treat NAs in 'accepted' vector as not accepted
    accepted[is.na(accepted)] <- 0
    # The first step is the initial model, so it's "accepted" by definition
    if(length(accepted) > 0) accepted[1] <- 1
    
    combined_scores[[i]] <- scores
    combined_moves[[i]] <- moves
    combined_accepted[[i]] <- accepted
    fit_lengths[i] <- length(scores)
  }
  
  # Interleave the data with NAs for breaks between rate category runs
  final_scores <- unlist(lapply(1:num_fits, function(i) {
    if (i < num_fits) c(combined_scores[[i]], rep(NA, break_size)) else combined_scores[[i]]
  }))
  final_moves <- unlist(lapply(1:num_fits, function(i) {
    if (i < num_fits) c(combined_moves[[i]], rep(NA, break_size)) else combined_moves[[i]]
  }))
  final_accepted <- unlist(lapply(1:num_fits, function(i) {
    if (i < num_fits) c(combined_accepted[[i]], rep(NA, break_size)) else combined_accepted[[i]]
  }))
  
  # Calculate positions for vertical separator lines
  vline_pos <- cumsum(fit_lengths)[-num_fits] + (1:(num_fits-1)) * break_size - break_size/2
  
  # 3. Plot Setup
  all_possible_moves <- unique(unlist(combined_moves, use.names = FALSE), na.rm = TRUE)
  
  # Ensure all moves in the data have a color assigned from the palette
  missing_moves <- setdiff(all_possible_moves, names(palette))
  if (length(missing_moves) > 0) {
    new_colors <- setNames(rep("grey50", length(missing_moves)), missing_moves)
    palette <- c(palette, new_colors)
    warning("Some move types were not in the default palette and have been assigned a default color.")
  }
  move_colors <- palette[all_possible_moves]
  
  if(all(is.na(final_scores))) {
    message("No scores available to plot.")
    return(invisible(NULL))
  }
  
  y_range <- diff(range(final_scores, na.rm = TRUE))
  y_bottom <- min(final_scores, na.rm = TRUE) - y_range * 0.25
  y_top <- max(final_scores, na.rm = TRUE) + y_range * 0.05
  total_len <- length(final_scores)
  
  # Default plot arguments, can be overridden by user with `...`
  plot_args <- list(x = seq_along(final_scores), y = final_scores, type = "n",
    ylim = c(y_bottom, y_top), xlim = c(0.5, total_len + 0.5),
    xlab = "Iteration", ylab = "Information Criterion Score",
    main = "Simulated Annealing Trace",
    las = 1, bty = "l", xaxs = "i")
  
  user_args <- list(...)
  for (arg_name in names(user_args)) {
    plot_args[[arg_name]] <- user_args[[arg_name]]
  }
  
  # Create the empty plot
  do.call(graphics::plot, plot_args)
  
  # 4. Plotting Elements
  # **FIXED LINE HERE**
  graphics::abline(h = pretty(final_scores, n = 8), col = "gray90", lty = 1, lwd = 0.5)
  
  # Add separator lines and labels for different rate category runs
  if (num_fits > 1) {
    graphics::abline(v = vline_pos, col = "gray40", lty = "dashed", lwd = 1.5)
    text_pos <- c(fit_lengths[1]/2, vline_pos + break_size/2 + fit_lengths[-1]/2)
    # Ensure there are labels for all fits
    if (length(rate_cat_labels) == num_fits){
      labels <- paste(rate_cat_labels, "Rate Class(es)")
      graphics::mtext(labels, side=1, line=2.5, at=text_pos, cex=0.9)
    }
  }
  
  graphics::lines(seq_along(final_scores), final_scores, lwd = 1.5, col = "gray50")
  
  accepted_idx <- which(final_accepted == 1)
  rejected_idx <- which(final_accepted == 0)
  
  graphics::points(accepted_idx, final_scores[accepted_idx], 
    col = "darkgreen", pch = 19, cex = 1.0)
  
  graphics::points(rejected_idx, final_scores[rejected_idx], 
    col = "firebrick", pch = 4, cex = 0.8, lwd = 1.5)
  
  # Draw the bar at the bottom indicating the move type at each iteration
  rect_height_fraction <- 0.15
  rect_bottom <- y_bottom + y_range * 0.05
  rect_top <- rect_bottom + y_range * rect_height_fraction
  
  for(i in seq_along(final_scores)) {
    if(!is.na(final_moves[i])) {
      graphics::rect(xleft = i - 0.5, 
        xright = i + 0.5,
        ybottom = rect_bottom,
        ytop = rect_top,
        col = move_colors[final_moves[i]],
        border = 'grey20')
    }
  }
  
  # Add legend
  if (legend) {
    legend_labels <- c("Accepted", "Rejected", names(move_colors))
    legend_colors <- c("darkgreen", "firebrick", move_colors)
    legend_pch <- c(19, 4, rep(15, length(move_colors)))
    legend_lwd <- c(NA, 1.5, rep(NA, length(move_colors)))
    legend_pt.cex <- c(1.0, 0.8, rep(1.5, length(move_colors)))
    
    graphics::legend(legend.pos,
      legend = legend_labels,
      col = legend_colors,
      pch = legend_pch,
      lwd = legend_lwd,
      pt.cex = legend_pt.cex,
      bty = "n",
      cex = 0.8)
  }
  
  return(invisible(NULL))
}

# old version that doesn't use simm ann
# corHMMDredge <- function(phy, data, max.rate.cat, root.p="maddfitz", 
#   pen.type = "l1", lambda = 0, drop.par = TRUE, drop.threshold = 1e-7,  
#   info.threshold=2, criterion="AIC", merge.params=TRUE, merge.threshold=0, 
#   rate.mat=NULL, node.states = "marginal", fixed.nodes=FALSE, ip=NULL, 
#   nstarts=0, n.cores=1, get.tip.states = FALSE, lewis.asc.bias = FALSE, 
#   collapse = FALSE, lower.bound = 1e-10, upper.bound = 100, opts=NULL, 
#   verbose=TRUE, p=NULL, rate.cat=NULL, grad=FALSE){
#   
#   if((is.null(p) & !is.null(rate.cat))){
#     print("A rate category was given without specifying a parameter vector (p)")
#     return(NULL)
#   }
#   if((!is.null(p) & is.null(rate.cat))){
#     print("A parameter vector (p) was given without specifying a rate category")
#     return(NULL)
#   }
#   
#   # FIXED FIT
#   init.root.p <- root.p
#   if(!is.null(p) & !is.null(rate.cat)){
#     if(verbose){
#       cat("Evaluating fixed parameters p =", p, "\n")
#     }
#     fixd_fit <- corHMMDredgeBase(phy=phy, 
#       data=data, 
#       rate.cat=rate.cat, 
#       root.p=root.p,
#       pen.type = pen.type, 
#       lambda = lambda, 
#       rate.mat=rate.mat, 
#       node.states = node.states, 
#       fixed.nodes=fixed.nodes, 
#       ip=ip, 
#       nstarts=nstarts, 
#       n.cores=n.cores, 
#       get.tip.states = get.tip.states, 
#       lewis.asc.bias = lewis.asc.bias, 
#       collapse = collapse, 
#       lower.bound = lower.bound, 
#       upper.bound = upper.bound, 
#       opts=opts, 
#       p=p,
#       grad=FALSE)
#     return(fixd_fit)
#   }
#   
#   # automatic search starting at rate cat 1
#   curr_index_mat <- rate.mat
#   fit_set <- list()
#   model_improved <- TRUE
#   hmm_valid <- TRUE
#   count <- 1
#   current_rate_category <- 1
#   if(verbose){
#     cat("Begining dredge...\n")
#   }
#   while(model_improved){
#     if(!inherits(root.p, "character")){
#       root.p <- rep(init.root.p, current_rate_category)
#     }
#     curr_fit <- try(corHMMDredgeBase(phy=phy, 
#       data=data, 
#       rate.cat=current_rate_category, 
#       root.p=root.p,
#       pen.type = pen.type, 
#       lambda = lambda, 
#       rate.mat=curr_index_mat, 
#       node.states = node.states, 
#       fixed.nodes=fixed.nodes, 
#       ip=ip, 
#       nstarts=nstarts, 
#       n.cores=n.cores, 
#       get.tip.states = get.tip.states, 
#       lewis.asc.bias = lewis.asc.bias, 
#       collapse = collapse, 
#       lower.bound = lower.bound, 
#       upper.bound = upper.bound, 
#       opts=opts, 
#       p=NULL,
#       grad=grad))
#     if(inherits(curr_fit, "try-error")){
#       warning("Model fitting failed. Stopping dredge.")
#       model_improved <- FALSE
#       fit_set[[count]] <- curr_fit
#       next
#     }
#     curr_info_criterion <- curr_fit[[criterion]]
#     fit_set[[count]] <- curr_fit
#     if(verbose){
#       cat("\n")
#       cat("AIC:", curr_info_criterion)
#       cat("\nMapping matrix:\n")
#       print(curr_fit$index.mat)
#     }
#     best_info_criterion <- get_best_info_criterion(fit_set, criterion, current_rate_category)
#     model_improved <- diff(c(best_info_criterion, curr_info_criterion)) < info.threshold
#     count <- count + 1
#     if(model_improved | current_rate_category < max.rate.cat){
#       # try dropping pars
#       curr_index_mat <- drop_pars(curr_fit, drop.threshold)
#       if(is.null(curr_index_mat)){
#         # try merging pars
#         curr_index_mat <- merge_pars(curr_fit, merge.threshold)
#         if(is.null(curr_index_mat)){
#           current_rate_category <- current_rate_category + 1
#           model_improved <- TRUE
#           if(current_rate_category > max.rate.cat){
#             model_improved <- FALSE
#           }else{
#             if(curr_fit$rate.cat > 1){
#               # test the current rate category
#               hmm_valid <- test_validity_hmm(curr_fit)
#               if(all(hmm_valid)){
#                 cat("\n", "No unique parameters in rate class. Halting dredge", "\n")
#                 model_improved <- FALSE
#               }
#             }else{
#               if(verbose){
#                 cat("\n", rep("*", count-1), "Continuing dredge", rep("*", count-1), "\n")
#               }
#             }
#           }
#         }else{
#           if(verbose){
#             cat("\n", rep("*", count-1), "Continuing dredge", rep("*", count-1), "\n")
#           }
#         }
#       }else{
#         if(verbose){
#           cat("\n", rep("*", count-1), "Continuing dredge", rep("*", count-1), "\n")
#         }
#       }
#     }
#   }
#   if(verbose){
#     cat("\nDone.\n")
#   }
#   tmp <- try(prune_redundant(fit_set))
#   if(class(tmp) != "try-error"){
#     fit_set <- tmp
#   }
#   class(fit_set) <- "corhmm.dredge"
#   return(fit_set)
# }
