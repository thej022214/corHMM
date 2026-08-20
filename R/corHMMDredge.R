### corHMM -- Generalized hidden Markov Models
# automatic fitting 
corHMMDredge <- function(phy, data, max.rate.cat=1, init.rate.cat=1, 
  root.p="maddfitz", tip.fog=NULL, fog.ip = 0.01, pen.type = "l1", lambda = 0, 
  drop.threshold = 1e-7, criterion="AIC", merge.threshold=0, index_mat=NULL, 
  node.states = "none", fixed.nodes=FALSE, ip=NULL, nstarts=1, n.cores=1, 
  get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = FALSE, lower.bound = 1e-10, 
  upper.bound = 100, opts=NULL, verbose=TRUE, p=NULL, rate.cat=NULL, use_RTMB=TRUE, 
  max.iterations = 1000, initial.temp = 2, cooling.rate = 0.95, 
  temp.schedule = "exponential", seed = NULL, require.entry = FALSE,
  stall.limit = 500, polish.restarts = 2, polish.max.fits = 5,
  checkpoint.file = NULL, checkpoint.interval = 50){

  if((is.null(p) & !is.null(rate.cat))){
    print("A rate category was given without specifying a parameter vector (p)")
    return(NULL)
  }
  if((!is.null(p) & is.null(rate.cat))){
    print("A parameter vector (p) was given without specifying a rate category")
    return(NULL)
  }
  criterion <- match.arg(criterion, c("AIC", "AICc", "BIC"))
  if(length(polish.max.fits) != 1 || !is.finite(polish.max.fits) ||
    polish.max.fits < 1) stop("polish.max.fits must be a positive number.",
      call. = FALSE)
  polish.max.fits <- as.integer(polish.max.fits)
  
  # Set seed for reproducibility
  if(!is.null(seed)) set.seed(seed)
  
  # FIXED FIT
  init.root.p <- root.p
  if(!is.null(p) & !is.null(rate.cat)){
    if(verbose){
      cat("Evaluating fixed parameters p =", p, "\n")
    }
    fixd_fit <- corHMMDredgeBase(phy=phy, data=data, rate.cat=rate.cat, 
      root.p=root.p, pen.type = pen.type, lambda = lambda, rate.mat=index_mat, 
      node.states = node.states, fixed.nodes=fixed.nodes, ip=ip, 
      nstarts=nstarts, n.cores=n.cores, get.tip.states = get.tip.states, 
      lewis.asc.bias = lewis.asc.bias, collapse = collapse, 
      lower.bound = lower.bound, upper.bound = upper.bound, 
      tip.fog=tip.fog, fog.ip=fog.ip, opts=opts, p=p, use_RTMB=use_RTMB)
    return(fixd_fit)
  }
  
  # Initialize
  default_index_mat <- getStateMat4Dat(data, collapse = collapse,
    indep = FALSE)$rate.mat
  default_index_mat[default_index_mat == 0] <- NA
  if(is.null(index_mat)){
    curr_index_mat <- default_index_mat
  } else {
    curr_index_mat <- index_mat
  }
  curr_index_mat[curr_index_mat == 0] <- NA
  initial_reason <- validate_index_mat(curr_index_mat, 1, require.entry)
  if(!is.null(initial_reason)) {
    stop("index_mat is outside the legal model space: ", initial_reason,
      call. = FALSE)
  }
  max_index_mat_cp <- curr_index_mat
  max_index_mat <- curr_index_mat
  max_index_mat[max_index_mat == 0] <- NA
  
  fit_set <- list()
  model_improved <- TRUE
  hmm_valid <- TRUE
  count <- 0
  # always fit the simplest model
  if(init.rate.cat > 1){
    warning("init.rate.cat ignored: the dredge always starts at rate.cat = 1")
  }
  current_rate_category <- 1
  
  
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
    
    curr_prep <- corHMMDredgePrep(phy=phy, data=data, rate.cat=current_rate_category,
      root.p=root.p, pen.type=pen.type, lambda=lambda, node.states=node.states,
      fixed.nodes=fixed.nodes, collapse=collapse, lower.bound=lower.bound,
      upper.bound=upper.bound, lewis.asc.bias=lewis.asc.bias,
      use_RTMB=use_RTMB, opts=opts, tip.fog=tip.fog, fog.ip=fog.ip)
    
    rc_checkpoint_file <- if (!is.null(checkpoint.file)) {
      paste0(tools::file_path_sans_ext(checkpoint.file), 
        "_rc", current_rate_category, ".",
        tools::file_ext(checkpoint.file))
    } else {
      NULL
    }
    
    # Check if checkpoint already has an initial fit we can reuse
    checkpoint_has_state <- !is.null(rc_checkpoint_file) && 
      file.exists(rc_checkpoint_file) && {
        ckpt <- tryCatch(readRDS(rc_checkpoint_file), error = function(e) NULL)
        checkpoint_is_compatible(ckpt, max_index_mat,
          current_rate_category, require.entry, criterion)
      }
    
    if (checkpoint_has_state) {
      if (verbose) cat("Checkpoint found for rate category", current_rate_category, 
        "- skipping initial fit.\n")
      curr_fit <- NULL  # sa_within_rate_category will restore from checkpoint
    } else {
      # Initial fit for this rate category
      curr_fit <- try(corHMMDredgeBase(phy=phy, data=data, 
        rate.cat=current_rate_category, root.p=root.p, pen.type = pen.type, 
        lambda = lambda, rate.mat=curr_index_mat, node.states = node.states, 
        fixed.nodes=fixed.nodes, ip=ip, nstarts=nstarts, tip.fog=tip.fog, 
        n.cores=n.cores, get.tip.states = get.tip.states, fog.ip=fog.ip,
        lewis.asc.bias = lewis.asc.bias, collapse = collapse, 
        lower.bound = lower.bound, upper.bound = upper.bound, 
        opts=opts, p=NULL, use_RTMB=use_RTMB, prep=curr_prep))
      
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
    }
    
    # SIMULATED ANNEALING WITHIN THIS RATE CATEGORY
    sa_result <- sa_within_rate_category(phy, data, curr_fit, 
      curr_index_mat, max_index_mat,
      current_rate_category, root.p, 
      pen.type, lambda, node.states, fixed.nodes, 
      ip, nstarts, n.cores, get.tip.states, 
      lewis.asc.bias, collapse, lower.bound, 
      upper.bound, opts, use_RTMB, criterion,
      drop.threshold, merge.threshold,
      max.iterations, initial.temp, cooling.rate, 
      temp.schedule, verbose, prep=curr_prep, require.entry=require.entry,
      stall.limit=stall.limit, polish.restarts=polish.restarts,
      polish.max.fits=polish.max.fits,
      tip.fog=tip.fog, fog.ip=fog.ip,
      checkpoint.file = rc_checkpoint_file, 
      checkpoint.interval = checkpoint.interval)
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
  
  accepted_models <- list()
  for(i in 1:length(fit_set)){
    accepted_models <- c(accepted_models, fit_set[[i]]$accepted_models)
  }
  accepted_models <- prune_redundant(accepted_models)
  # reattach to final models the big data pieces
  for(i in seq_along(accepted_models)){
    accepted_models[[i]]$phy         <- phy
    accepted_models[[i]]$data        <- data
    accepted_models[[i]]$data.legend <- curr_fit$data.legend
  }
  class(accepted_models) <- "corhmm.dredge"
  attr(accepted_models, "criterion") <- criterion
  attr(accepted_models, "dredge_history") <- fit_set
  attr(accepted_models, "flat_directions") <- flat_directions(
    unlist(lapply(fit_set, function(x) x$every_model), recursive = FALSE))
  return(accepted_models)
}

strip_corhmm <- function(fit) {
  fit$phy         <- NULL
  fit$data        <- NULL
  fit$data.legend <- NULL
  fit
}

information_criteria <- function(loglik, np, ntax) {
  list(
    AIC = -2 * loglik + 2 * np,
    AICc = -2 * loglik + 2 * np * ntax / (ntax - np - 1),
    BIC = -2 * loglik + log(ntax) * np
  )
}

# An exact move: same likelihood, fewer parameters
reindex_fit <- function(fit, index_mat) {
  old_np <- fit$AIC / 2 + fit$loglik
  ntax <- old_np + 1 + 2 * old_np * (old_np + 1) / (fit$AICc - fit$AIC)
  np <- old_np - max(fit$index.mat, na.rm = TRUE) + max(index_mat, na.rm = TRUE)
  fit$index.mat <- index_mat
  fit[names(information_criteria(fit$loglik, np, ntax))] <-
    information_criteria(fit$loglik, np, ntax)
  fit
}

# SA within a single rate category
sa_within_rate_category <- function(phy, data, initial_fit, initial_index_mat, 
  max_index_mat, rate_category, root.p, pen.type, lambda, node.states, 
  fixed.nodes, ip, nstarts, n.cores, get.tip.states, lewis.asc.bias, collapse, 
  lower.bound, upper.bound, opts, use_RTMB, criterion, drop.threshold, 
  merge.threshold, max.iterations, initial.temp, cooling.rate, temp.schedule, 
  verbose, restart.strategy = "fixed_steps", tip.fog=tip.fog, fog.ip=fog.ip,
  restart.interval = 20, restart.threshold = 1.5, restart.probability = 0.05,
  restart.temp.reset = FALSE, prep=NULL, require.entry = FALSE, stall.limit = 500,
  polish.restarts = 2, polish.max.fits = 5,
  checkpoint.file = NULL, checkpoint.interval = 50) {
  
  checkpoint_id <- model_space_id(max_index_mat, rate_category)
  ckpt <- NULL
  if (!is.null(checkpoint.file) && file.exists(checkpoint.file)) {
    ckpt <- tryCatch(readRDS(checkpoint.file), error = function(e) NULL)
    checkpoint_ok <- checkpoint_is_compatible(ckpt, max_index_mat,
      rate_category, require.entry, criterion)
    if(!checkpoint_ok) {
      warning("Ignoring incompatible dredge checkpoint; starting fresh.",
        immediate. = TRUE)
      ckpt <- NULL
    }
  }
  if (!is.null(ckpt)) {
    if (verbose) cat("Resuming from checkpoint:", checkpoint.file, "\n")
    ## unpack checkpoint state and return early if already finished
    if (ckpt$finished) {
      if (verbose) cat("Checkpoint marks run as finished. Returning saved result.\n")
      return(ckpt$result)
    }
    ## restore all SA state
    current_fit        <- ckpt$current_fit
    current_index_mat  <- ckpt$current_index_mat
    current_score      <- ckpt$current_score
    best_fit           <- ckpt$best_fit
    best_index_mat     <- ckpt$best_index_mat
    best_score         <- ckpt$best_score
    accepted_models    <- ckpt$accepted_models
    accepted_index_mats <- ckpt$accepted_index_mats
    accepted_scores    <- ckpt$accepted_scores
    model_ids          <- ckpt$model_ids
    every_model        <- ckpt$every_model
    every_score        <- ckpt$every_score
    every_move         <- ckpt$every_move
    accepted           <- ckpt$accepted
    score_cache        <- list2env(ckpt$score_cache, hash = TRUE)
    rejections         <- ckpt$rejections
    accepted_moves     <- ckpt$accepted_moves
    total_moves        <- ckpt$total_moves
    restart_count      <- ckpt$restart_count
    steps_since_restart <- ckpt$steps_since_restart
    steps_since_accept <- ckpt$steps_since_accept
    steps_since_best   <- ckpt$steps_since_best
    steps_since_new    <- ckpt$steps_since_new
    current_temp       <- ckpt$current_temp
    temperature_epoch_start <- if(is.null(ckpt$temperature_epoch_start)) 1L else ckpt$temperature_epoch_start
    polish_runs        <- if(is.null(ckpt$polish_runs)) 0L else ckpt$polish_runs
    polish_evaluations <- if(is.null(ckpt$polish_evaluations)) 0L else ckpt$polish_evaluations
    polish_failures    <- if(is.null(ckpt$polish_failures)) 0L else ckpt$polish_failures
    loglik_sat         <- ckpt$loglik_sat
    start_iteration    <- ckpt$iteration + 1L
    every_idx          <- ckpt$every_idx
  } else {
    ## --- fresh start ---
    current_fit        <- initial_fit
    score_cache        <- new.env(hash = TRUE)
    rejections         <- character()
    current_index_mat  <- initial_index_mat
    current_score      <- current_fit[[criterion]]
    
    accepted_models     <- list()
    every_model    <- list()
    accepted_index_mats <- list()
    every_score    <- accepted_scores <- numeric()
    model_ids      <- character()
    every_move     <- character()
    accepted       <- numeric()
    
    initial_id <- canonical_id(current_fit$index.mat, rate_category)
    score_cache[[initial_id]] <- strip_corhmm(initial_fit)
    every_model[[1]]    <- accepted_models[[1]] <- initial_fit
    accepted_index_mats[[1]] <- initial_index_mat
    every_score[1]      <- accepted_scores[1] <- current_score
    model_ids[1]        <- initial_id
    every_move[1]       <- "none"
    accepted[1]         <- 1L
    every_idx           <- 1L   
    
    best_fit           <- current_fit
    best_index_mat     <- current_index_mat
    best_score         <- current_score
    ## upper bound on the likelihood of any submodel
    loglik_sat         <- initial_fit$loglik
    current_temp       <- initial.temp
    accepted_moves     <- 0L
    total_moves        <- 0L
    restart_count      <- 0L
    steps_since_restart <- 0L
    steps_since_accept  <- 0L
    steps_since_best    <- 0L
    steps_since_new     <- 0L
    temperature_epoch_start <- 1L
    polish_runs         <- 0L
    polish_evaluations  <- 0L
    polish_failures     <- 0L
    start_iteration    <- 1L
  }
  
  if (verbose) {
    cat("Starting SA with", criterion, "=", round(current_score, 3), "\n")
    cat("Restart strategy:", restart.strategy, "\n")
  }
  
  ## why a proposal died, tallied per move type
  reject <- function(reason) {
    rejections <<- c(rejections, paste(move_result$move_type, reason, sep = ": "))
  }
  
  ## helper to save checkpoint
  save_checkpoint <- function(iteration, finished = FALSE, result = NULL) {
    if (is.null(checkpoint.file)) return(invisible(NULL))
    saveRDS(list(
      version            = 4L,
      model_space_id     = checkpoint_id,
      criterion          = criterion,
      finished           = finished,
      result             = result,
      iteration          = iteration,
      current_fit        = current_fit,
      current_index_mat  = current_index_mat,
      current_score      = current_score,
      best_fit           = best_fit,
      best_index_mat     = best_index_mat,
      best_score         = best_score,
      accepted_models         = accepted_models,
      accepted_index_mats     = accepted_index_mats,
      accepted_scores         = accepted_scores,
      model_ids          = model_ids,
      every_model        = every_model,
      every_score        = every_score,
      every_move         = every_move,
      accepted           = accepted,
      score_cache        = as.list(score_cache),
      rejections         = rejections,
      accepted_moves     = accepted_moves,
      total_moves        = total_moves,
      restart_count      = restart_count,
      steps_since_restart = steps_since_restart,
      steps_since_accept = steps_since_accept,
      steps_since_best   = steps_since_best,
      steps_since_new    = steps_since_new,
      current_temp       = current_temp,
      temperature_epoch_start = temperature_epoch_start,
      polish_runs        = polish_runs,
      polish_evaluations = polish_evaluations,
      polish_failures    = polish_failures,
      loglik_sat         = loglik_sat,
      every_idx          = every_idx
    ), checkpoint.file)
  }
  
  ## number of tip fog parameters, constant across models
  fog_np <- current_fit$AIC / 2 + current_fit$loglik -
    max(current_fit$index.mat, na.rm = TRUE)
  
  ## preps and root.p are rate category specific
  preps <- list()
  preps[[rate_category]] <- prep
  get_root_p <- function(rc) if (is.character(root.p)) root.p else rep(root.p[1], rc)
  get_prep <- function(rc) {
    if (is.null(preps[[rc]])) {
      preps[[rc]] <<- corHMMDredgePrep(phy=phy, data=data, rate.cat=rc,
        root.p=get_root_p(rc), pen.type=pen.type, lambda=lambda,
        node.states=node.states, fixed.nodes=fixed.nodes, collapse=collapse,
        lower.bound=lower.bound, upper.bound=upper.bound,
        lewis.asc.bias=lewis.asc.bias, use_RTMB=use_RTMB, opts=opts,
        tip.fog=tip.fog, fog.ip=fog.ip)
    }
    preps[[rc]]
  }

  polish_neighbor_fn <- function(fit) {
    rc <- fit$rate.cat
    deterministic_neighbors(fit, sub_max(max_index_mat, rc, rate_category),
      rate_category, merge.threshold, require.entry)
  }

  fit_polish_candidate <- function(candidate, fixed_rates = NULL,
    initial_rates = ip) {
    candidate_rc <- candidate$rate_cat
    try(corHMMDredgeBase(phy=phy, data=data, rate.cat=candidate_rc,
      root.p=get_root_p(candidate_rc), pen.type=pen.type, lambda=lambda,
      rate.mat=candidate$new_index_mat, node.states=node.states,
      fixed.nodes=fixed.nodes, ip=initial_rates, nstarts=nstarts,
      n.cores=n.cores, get.tip.states=get.tip.states,
      lewis.asc.bias=lewis.asc.bias, collapse=collapse, fog.ip=fog.ip,
      lower.bound=lower.bound, upper.bound=upper.bound, tip.fog=tip.fog,
      opts=opts, p=fixed_rates, use_RTMB=use_RTMB,
      prep=get_prep(candidate_rc)))
  }

  polish_score_fn <- function(candidate) {
    candidate_rc <- candidate$rate_cat
    id <- canonical_id(candidate$new_index_mat, candidate_rc)
    fit <- score_cache[[id]]
    fitted <- FALSE
    source_score <- candidate$source_scores[[criterion]]
    if(identical(candidate$move_type, "polish_cell_drop") &&
      !is.null(candidate$initial_rates) &&
      !isTRUE(get_prep(candidate_rc)$set.fog) &&
      (is.null(fit) || fit[[criterion]] > source_score + 1e-8)) {
      fit <- NULL
      inherited_fit <- fit_polish_candidate(candidate,
        fixed_rates=candidate$initial_rates)
      if(!inherits(inherited_fit, "try-error") &&
        !is.null(inherited_fit[[criterion]]) &&
        inherited_fit[[criterion]] <= source_score + 1e-8) {
        fit <- inherited_fit
        fitted <- TRUE
      }
    }
    if(is.null(fit)) {
      initial_rates <- if(isTRUE(get_prep(candidate_rc)$set.fog) ||
        is.null(candidate$initial_rates)) ip else
        candidate$initial_rates * get_prep(candidate_rc)$H
      fit <- fit_polish_candidate(candidate, initial_rates=initial_rates)
      fitted <- TRUE
    }
    if(fitted) {
      if(inherits(fit, "try-error")) return(NULL)
      if(is.null(fit$loglik) || length(fit$loglik) == 0) return(NULL)
      if(fit$loglik == -1e+06) return(NULL)
      fit <- strip_corhmm(fit)
      score_cache[[id]] <- fit
      steps_since_new <<- 0L
      polish_evaluations <<- polish_evaluations + 1L
      every_idx <<- every_idx + 1L
      every_model[[every_idx]] <<- fit
      every_score[every_idx] <<- fit[[criterion]]
      every_move[every_idx] <<- candidate$move_type
      accepted[every_idx] <<- 0L
    }
    fit
  }

  run_polish <- function(fit) {
    polish_runs <<- polish_runs + 1L
    if(verbose) cat("Polishing best", criterion, ":",
      round(fit[[criterion]], 3), "- fitting at most", polish.max.fits,
      "new neighbors.\n")
    result <- polish_neighborhood(fit, polish_neighbor_fn, polish_score_fn,
      criterion, max.new.evaluations = polish.max.fits,
      is.cached = function(candidate) {
        id <- canonical_id(candidate$new_index_mat, candidate$rate_cat)
        !is.null(score_cache[[id]])
      })
    if(verbose && result$deferred > 0) cat("Polishing deferred",
      result$deferred, "neighbors after", result$new_evaluations,
      "new fits.\n")
    result
  }

  run_cell_cleanup <- function(fit) {
    polish_neighborhood(fit, function(current) {
      Filter(function(candidate)
        identical(candidate$move_type, "polish_cell_drop"),
        polish_neighbor_fn(current))
    }, polish_score_fn, criterion)
  }

  accept_polish_path <- function(path) {
    if(length(path) == 0) return()
    for(fit in path) {
      id <- canonical_id(fit$index.mat, fit$rate.cat)
      model_count <- length(accepted_models) + 1L
      accepted_models[[model_count]] <<- fit
      accepted_index_mats[[model_count]] <<- fit$index.mat
      accepted_scores[model_count] <<- fit[[criterion]]
      model_ids[model_count] <<- id
      seen <- which(vapply(every_model, function(x)
        !is.null(x) && identical(canonical_id(x$index.mat, x$rate.cat), id),
        logical(1)))
      if(length(seen)) accepted[max(seen)] <<- 1L
    }
  }
  
  stop_reason <- "max.iterations"
  for (iteration in start_iteration:max.iterations) {
    steps_since_restart <- steps_since_restart + 1L
    steps_since_accept  <- steps_since_accept  + 1L
    steps_since_best    <- steps_since_best    + 1L
    steps_since_new     <- steps_since_new     + 1L

    if (polish_due(steps_since_best, steps_since_new, stall.limit,
      restart.interval)) {
      polished <- run_polish(best_fit)
      if(dredge_fit_is_better(polished$best_fit, best_fit, criterion)) {
        accept_polish_path(polished$path)
        best_fit <- current_fit <- polished$best_fit
        best_index_mat <- current_index_mat <- best_fit$index.mat
        best_score <- current_score <- best_fit[[criterion]]
        steps_since_best <- 0L
        steps_since_new <- 0L
        steps_since_accept <- 0L
        steps_since_restart <- 0L
        temperature_epoch_start <- iteration
        current_temp <- initial.temp
        polish_failures <- 0L
        if(verbose) cat("Polishing improved", criterion, "to", round(best_score, 3), "- reheating.\n")
        next
      }
      if(!polished$complete) {
        current_fit <- best_fit
        current_index_mat <- best_index_mat
        current_score <- best_score
        steps_since_best <- 0L
        steps_since_new <- 0L
        steps_since_accept <- 0L
        steps_since_restart <- 0L
        temperature_epoch_start <- iteration
        current_temp <- initial.temp
        if(verbose) cat("Returning to SA before the next polishing batch.\n")
        next
      }
      polish_failures <- polish_failures + 1L
      if(polish_failures > polish.restarts) {
        stop_reason <- "local optimum"
        if(verbose) cat("No improving neighbor after", polish_failures, "polish runs.\n")
        break
      }
      current_fit <- best_fit
      current_index_mat <- best_index_mat
      current_score <- best_score
      steps_since_best <- 0L
      steps_since_new <- 0L
      steps_since_accept <- 0L
      steps_since_restart <- 0L
      temperature_epoch_start <- iteration
      current_temp <- initial.temp
      if(verbose) cat("No improving neighbor - reheating", polish_failures, "of", polish.restarts, ".\n")
      next
    }
    
    
    ## --- periodic checkpoint ---
    if (!is.null(checkpoint.file) && iteration %% checkpoint.interval == 0) {
      save_checkpoint(iteration)
      if (verbose) cat("  [Checkpoint saved at iteration", iteration, "]\n")
    }
    
    ## --- restart logic ---
    should_restart <- FALSE
    restart_reason <- ""
    if (restart.strategy == "fixed_steps" &&
        steps_since_best    >= restart.interval &&
        steps_since_restart >= restart.interval &&
        steps_since_accept  >= restart.interval) {
      should_restart <- TRUE; restart_reason <- "fixed_steps"
    } else if (restart.strategy == "energy_threshold" &&
        current_score > best_score * restart.threshold) {
      should_restart <- TRUE; restart_reason <- "energy_threshold"
    } else if (restart.strategy == "random" && runif(1) < restart.probability) {
      should_restart <- TRUE; restart_reason <- "random"
    }
    
    if (should_restart && restart.strategy != "none") {
      current_id <- canonical_id(current_fit$index.mat, current_fit$rate.cat)
      available_indices <- which(model_ids != current_id)
      if(length(available_indices) == 0) {
        steps_since_new <- restart.interval
        next
      }
      every_idx <- every_idx + 1L
      every_move[every_idx]  <- "restart"
      every_score[every_idx] <- current_score
      every_model[every_idx] <- list(NULL)   ## no new fit on a restart
      accepted[every_idx]    <- NA_integer_
      
      restart_idx       <- sample(available_indices, 1)
      current_fit       <- accepted_models[[restart_idx]]
      current_index_mat <- accepted_index_mats[[restart_idx]]
      current_score     <- accepted_scores[restart_idx]
      if (restart.temp.reset) {
        current_temp <- initial.temp
        temperature_epoch_start <- iteration
      }
      restart_count       <- restart_count + 1L
      steps_since_restart <- 0L
      if (verbose) cat("RESTART", restart_count, "at iteration", iteration,
        "(", restart_reason, ") - Back to a previous model", criterion,
        ":", round(current_score, 3), "\n")
      next
    }
    
    current_temp <- annealing_temperature(iteration, max.iterations,
      initial.temp, cooling.rate, temp.schedule, temperature_epoch_start)
    
    ## --- propose move ---
    current_rc  <- current_fit$rate.cat
    current_max <- sub_max(max_index_mat, current_rc, rate_category)
    move_result <- propose_sa_move_within_rate_cat(current_fit, drop.threshold,
      merge.threshold, current_max)
    proposed_rc <- if (is.null(move_result$rate_cat)) current_rc else move_result$rate_cat

    proposed_max <- sub_max(max_index_mat, proposed_rc, rate_category)
    reason <- validate_index_mat(move_result$new_index_mat, proposed_rc,
      require.entry, proposed_max)
    if (!is.null(reason)) { reject(reason); next }
    
    move_id <- canonical_id(move_result$new_index_mat, proposed_rc)
    total_moves <- total_moves + 1L
    
    ## --- a structure seen before is scored, not refit ---
    proposed_fit <- score_cache[[move_id]]
    cache_hit <- !is.null(proposed_fit)
    
    if (is.null(proposed_fit)) {
      ## saturated bound: no submodel can beat this
      if (lambda == 0 && criterion %in% c("AIC", "BIC") && !is.null(loglik_sat) &&
          proposed_rc == rate_category) {
        k <- max(move_result$new_index_mat, na.rm = TRUE) + fog_np
        optimistic_score <- if(criterion == "AIC") -2 * loglik_sat + 2 * k else
          -2 * loglik_sat + log(length(phy$tip.label)) * k
        if (optimistic_score > best_score) { reject("saturated bound"); next }
      }
      
      ## an exact move keeps the current likelihood, so it needs no fit
      if (isTRUE(move_result$exact)) {
        proposed_fit <- reindex_fit(current_fit, move_result$new_index_mat)
      } else {
        proposed_fit <- try(corHMMDredgeBase(phy=phy, data=data, rate.cat=proposed_rc,
          root.p=get_root_p(proposed_rc), pen.type=pen.type, lambda=lambda,
          rate.mat=move_result$new_index_mat, node.states=node.states,
          fixed.nodes=fixed.nodes, ip=ip, nstarts=nstarts,
          n.cores=n.cores, get.tip.states=get.tip.states,
          lewis.asc.bias=lewis.asc.bias, collapse=collapse, fog.ip=fog.ip,
          lower.bound=lower.bound, upper.bound=upper.bound,tip.fog=tip.fog,
          opts=opts, p=NULL, use_RTMB=use_RTMB, prep=get_prep(proposed_rc)))
      }
      
      if (inherits(proposed_fit, "try-error")) { reject("fit failed"); next }
      if (is.null(proposed_fit$loglik) || length(proposed_fit$loglik) == 0) { reject("fit failed"); next }
      if (proposed_fit$loglik == -1e+06) { reject("no likelihood"); next }
      score_cache[[move_id]] <- strip_corhmm(proposed_fit)
      steps_since_new <- 0L
    }
    
    proposed_score         <- proposed_fit[[criterion]]
    if(!cache_hit) {
      every_idx              <- every_idx + 1L
      every_model[[every_idx]] <- strip_corhmm(proposed_fit)
      every_score[every_idx]   <- proposed_score
      every_move[every_idx]    <- move_result$move_type
    }
    
    ## --- accept/reject ---
    ## models within 2 AIC are indistinguishable, so the chain moves freely
    ## across that band and anneals only on the excess
    delta       <- proposed_score - current_score
    accept_prob <- if (delta <= 2) 1.0 else exp(-(delta - 2) / current_temp)
    
    if (runif(1) < accept_prob) {
      if(!cache_hit) accepted[every_idx] <- 1L
      current_fit          <- proposed_fit
      current_index_mat    <- current_fit$index.mat
      current_score        <- proposed_score
      accepted_moves       <- accepted_moves + 1L
      steps_since_accept   <- 0L
      
      if(!(move_id %in% model_ids)) {
        model_count            <- length(accepted_models) + 1L
        accepted_models[[model_count]]     <- strip_corhmm(current_fit)
        accepted_index_mats[[model_count]] <- current_index_mat
        accepted_scores[model_count]       <- current_score
        model_ids[model_count]        <- move_id
      }
      
      if (dredge_fit_is_better(current_fit, best_fit, criterion)) {
        best_fit        <- current_fit
        best_index_mat  <- current_index_mat
        best_score      <- current_score
        steps_since_best <- 0L
        polish_failures <- 0L
      }
      if (verbose) {
        cat("Iter", iteration, "-", if(cache_hit) "Revisited" else "New",
          paste0(criterion, ":"), round(current_score, 3),
          "- Best", paste0(criterion, ":"), round(best_score, 3), "\n",
          "Move:", move_result$move_type, "- Temp:", round(current_temp, 4),
          "- Steps since restart:", steps_since_restart, "\n",
          "Index Matrix:\n")
        print(current_index_mat)
        cat("\n")
      }
    } else {
      if(!cache_hit) accepted[every_idx] <- 0L
    }
  }

  if(stop_reason == "max.iterations") {
    polished <- run_polish(best_fit)
    if(dredge_fit_is_better(polished$best_fit, best_fit, criterion)) {
      accept_polish_path(polished$path)
      best_fit <- current_fit <- polished$best_fit
      best_index_mat <- current_index_mat <- best_fit$index.mat
      best_score <- current_score <- best_fit[[criterion]]
      steps_since_best <- 0L
    }
    stop_reason <- if(polished$complete) "max.iterations (polished)" else
      "max.iterations (partial polish)"
  }

  cleaned <- run_cell_cleanup(best_fit)
  if(dredge_fit_is_better(cleaned$best_fit, best_fit, criterion)) {
    accept_polish_path(cleaned$path)
    best_fit <- current_fit <- cleaned$best_fit
    best_index_mat <- current_index_mat <- best_fit$index.mat
    best_score <- current_score <- best_fit[[criterion]]
    if(verbose) cat("Final cell cleanup reduced the best model to",
      sum(!is.na(best_index_mat)), "transitions.\n")
  }

  if(!fits_in_model_space(c(accepted_models, every_model), max_index_mat,
    rate_category, require.entry)) {
    stop("Internal error: dredge produced a model outside the permitted transition mask.",
      call. = FALSE)
  }
  
  acceptance_rate <- if (total_moves > 0) accepted_moves / total_moves else 0
  
  model_summary <- data.frame(
    model_id = model_ids,
    score    = accepted_scores,
    stringsAsFactors = FALSE
  )
  model_summary <- model_summary[order(model_summary$score), ]
  
  result <- list(
    best_fit             = best_fit,
    best_index_mat       = best_index_mat,
    best_score           = best_score,
    accepted_models      = accepted_models,
    accepted_index_mats  = accepted_index_mats,
    accepted_scores      = accepted_scores,
    model_ids            = model_ids,
    model_summary        = model_summary,
    iterations           = iteration,
    acceptance_rate      = acceptance_rate,
    total_unique_models  = length(accepted_models),
    rate_category        = rate_category,
    criterion            = criterion,
    stop_reason          = stop_reason,
    total_moves          = total_moves,
    rejections           = table(rejections),
    unique_structures    = length(score_cache),
    restart_count        = restart_count,
    polish_runs           = polish_runs,
    polish_evaluations    = polish_evaluations,
    polish_failures       = polish_failures,
    final_temperature     = current_temp,
    final_steps_since_best = steps_since_best,
    every_model          = every_model,
    every_score          = every_score,
    every_move           = every_move,
    accepted             = accepted
  )
  
  ## save final checkpoint marking run as complete
  save_checkpoint(iteration, finished = TRUE, result = result)
  
  return(result)
}

# Propose stochastic moves within rate category
propose_sa_move_within_rate_cat <- function(current_fit, drop.threshold,
  merge.threshold, max_index_mat) {
  
  moves <- c("drop", "merge", "free", "lump")
  probs <- c(0.3, 0.3, 0.2, 0.2)
  if (current_fit$rate.cat > 1) {
    moves <- c(moves, "collapse")
    probs <- c(0.275, 0.275, 0.2, 0.2, 0.05)
  }
  move_type <- sample(moves, 1, prob = probs)
  
  if (move_type == "drop") {
    return(c(list(move_type = move_type),
      propose_stochastic_drop(current_fit, drop.threshold)))
  }
  if (move_type == "merge") {
    return(list(move_type = move_type,
      new_index_mat = propose_stochastic_merge(current_fit, merge.threshold)))
  }
  if (move_type == "free") {
    return(list(move_type = move_type,
      new_index_mat = propose_stochastic_free(current_fit, max_index_mat)))
  }
  if (move_type == "lump") {
    # exact: the likelihood is unchanged, so no refit is needed
    return(c(list(move_type = move_type), propose_lumpable_merge(current_fit)))
  }
  return(c(list(move_type = move_type), propose_collapse_class(current_fit)))
}

# Stochastic parameter dropping
propose_stochastic_drop <- function(current_fit, drop.threshold) {
  Q <- current_fit$solution
  Q[is.na(Q)] <- Inf
  index_mat <- current_fit$index.mat
  
  # parameter values, one per free parameter
  pars <- sort(unique(na.omit(as.vector(index_mat))))
  par_values <- sapply(pars, function(i) min(Q[which(index_mat == i)]))
  small_pars <- pars[par_values < drop.threshold]
  
  # nothing tiny, so consider the smallest quartile
  if(length(small_pars) == 0) {
    if(length(pars) == 0) return(NULL)
    small_pars <- pars[par_values <= quantile(par_values, 0.25)]
  }
  if(length(small_pars) == 0) return(NULL)
  
  drop_probs <- 1 / (par_values[match(small_pars, pars)] + 1e-10)
  n_to_drop <- sample(1:min(3, length(small_pars)), 1)
  if(length(small_pars) == 1){
    to_drop <- small_pars
  }else{
    to_drop <- sample(small_pars, n_to_drop, prob = drop_probs / sum(drop_probs))
  }
  new_index_mat <- dropStateMatPars(index_mat, to_drop)
  return(prune_isolated_states(new_index_mat, current_fit$rate.cat))
}

# A drop that strands a whole hidden class removes it and lowers the rate
# category
prune_isolated_states <- function(index_mat, rate_cat) {
  isolated <- vapply(seq_len(nrow(index_mat)), function(i)
    all(is.na(index_mat[i, ])) && all(is.na(index_mat[, i])), logical(1))
  if(!any(isolated)) return(list(new_index_mat = renumber_index_mat(index_mat)))
  
  n <- nrow(index_mat) / rate_cat
  class_of <- rep(seq_len(rate_cat), each = n)
  gone <- which(vapply(split(isolated, class_of), all, logical(1)))
  # nothing to collapse: hand it on and let the validator name the problem
  if(length(gone) == 0 || length(gone) == rate_cat) {
    return(list(new_index_mat = renumber_index_mat(index_mat)))
  }
  
  rc <- rate_cat - length(gone)
  new_index_mat <- renumber_index_mat(index_mat[!class_of %in% gone, !class_of %in% gone])
  rownames(new_index_mat) <- colnames(new_index_mat) <-
    paste0("(", rep(seq_len(n), rc), ",R", rep(seq_len(rc), each = n), ")")
  list(new_index_mat = new_index_mat, rate_cat = rc)
}

# Delete a whole hidden class and drop the rate category by one
propose_collapse_class <- function(current_fit) {
  rc <- current_fit$rate.cat
  index_mat <- current_fit$index.mat
  n <- nrow(index_mat) / rc
  Q <- current_fit$solution
  Q[is.na(Q)] <- 0
  class_of <- rep(seq_len(rc), each = n)
  
  dead <- which(sapply(seq_len(rc), function(k)
    all(Q[class_of != k, class_of == k] == 0)))
  if (length(dead) == 0) return(list(new_index_mat = NULL))
  gone <- dead[sample(length(dead), 1)]
  keep <- setdiff(seq_len(nrow(index_mat)), ((gone - 1) * n + 1):(gone * n))
  new_index_mat <- renumber_index_mat(index_mat[keep, keep])
  rownames(new_index_mat) <- colnames(new_index_mat) <-
    paste0("(", rep(seq_len(n), rc - 1), ",R", rep(seq_len(rc - 1), each = n), ")")
  list(new_index_mat = new_index_mat, rate_cat = rc - 1)
}

# Kemeny-Snell lumpability: two states lump exactly when their rates into
# every other state agree. Equating those rates costs no likelihood.
propose_lumpable_merge <- function(current_fit) {
  Q <- current_fit$solution
  Q[is.na(Q)] <- 0
  index_mat <- current_fit$index.mat
  n <- nrow(Q)
  if (n < 3 || max(index_mat, na.rm = TRUE) < 2) return(list(new_index_mat = NULL))
  tol <- 1e-8 * max(Q)
  
  candidates <- combn(n, 2)
  candidates <- candidates[, sample(ncol(candidates)), drop = FALSE]
  for (k in seq_len(ncol(candidates))) {
    i <- candidates[1, k]; j <- candidates[2, k]
    others <- setdiff(seq_len(n), c(i, j))
    if (any(abs(Q[i, others] - Q[j, others]) > tol)) next
    # only structurally equatable when both cells are free
    pairs <- others[!is.na(index_mat[i, others]) & !is.na(index_mat[j, others])]
    pairs <- pairs[index_mat[i, pairs] != index_mat[j, pairs]]
    if (length(pairs) == 0) next
    if (any(xor(is.na(index_mat[i, others]), is.na(index_mat[j, others])))) next
    new_index_mat <- equateStateMatPars(index_mat,
      lapply(pairs, function(m) c(index_mat[i, m], index_mat[j, m])))
    return(list(new_index_mat = renumber_index_mat(new_index_mat), exact = TRUE))
  }
  return(list(new_index_mat = NULL))
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
  return(renumber_index_mat(new_index_mat))
}

propose_stochastic_free <- function(current_fit, max_index_mat) {
  if(!identical(dim(current_fit$index.mat), dim(max_index_mat))) return(NULL)
  max_index_mat[max_index_mat == 0] <- NA
  diag(max_index_mat) <- NA
  duplicates <- !is.na(current_fit$index.mat) & 
    !is.na(max_index_mat) &
    duplicated(current_fit$index.mat, MARGIN = 0)
  dropped <- is.na(current_fit$index.mat) & !is.na(max_index_mat)
  if(sum(dropped | duplicates) == 0) return(NULL)
  candidates <- which(dropped | duplicates)
  n_free <- sample.int(min(3, length(candidates)), 1)
  focal_free <- candidates[sample.int(length(candidates), n_free)]
  new_index_mat <- current_fit$index.mat
  new_index_mat[focal_free] <- max(current_fit$index.mat, na.rm = TRUE)+1:n_free
  return(renumber_index_mat(new_index_mat))
}

# Stochastic version of merge_current_pars
stochastic_merge_pars <- function(current_pars, merge.threshold) {
  if(length(current_pars) <= 1) return(NULL)

  estimated <- which(current_pars > 1e-8)
  if(length(estimated) < length(current_pars)){
    merger <- stochastic_merge_pars(current_pars[estimated], merge.threshold)
    return(if(is.null(merger)) NULL else estimated[merger])
  }
  
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

# Models sharing a likelihood but not a parameter count mark a flat direction
flat_directions <- function(model_list){
  model_list <- model_list[!sapply(model_list, is.null)]
  if(length(model_list) < 2) return(NULL)
  lnLik <- round(sapply(model_list, "[[", "loglik"), 6)
  np <- sapply(model_list, function(x) max(x$index.mat, na.rm = TRUE))
  groups <- split(np, lnLik)
  groups <- groups[sapply(groups, function(x) length(unique(x)) > 1)]
  if(length(groups) == 0) return(NULL)
  data.frame(lnLik = as.numeric(names(groups)),
    n_models = sapply(groups, length),
    np = sapply(groups, function(x) paste(range(x), collapse = "-")),
    row.names = NULL)
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

# Relabel parameters by first appearance in a row-major scan
renumber_index_mat <- function(index_mat){
  labs <- unique(na.omit(c(t(index_mat))))
  out <- index_mat
  out[] <- match(index_mat, labs)
  out
}

perms <- function(x){
  if(length(x) == 1) return(matrix(x))
  do.call(rbind, lapply(seq_along(x), function(i) cbind(x[i], perms(x[-i]))))
}

# Hash invariant to hidden class relabelling
canonical_id <- function(index_mat, rate_cat){
  if(rate_cat < 2){
    return(paste0(c(renumber_index_mat(index_mat)), collapse = "_"))
  }
  n <- nrow(index_mat) / rate_cat
  ids <- apply(perms(seq_len(rate_cat)), 1, function(perm){
    ord <- as.vector(sapply(perm, function(k) ((k - 1) * n + 1):(k * n)))
    paste0(c(renumber_index_mat(index_mat[ord, ord])), collapse = "_")
  })
  min(ids)
}

# The upper bound matrix restricted to the first rc rate classes
sub_max <- function(max_index_mat, rc, full_rc){
  if(rc >= full_rc) return(max_index_mat)
  n <- nrow(max_index_mat) / full_rc
  keep <- seq_len(rc * n)
  renumber_index_mat(max_index_mat[keep, keep])
}

model_space_id <- function(max_index_mat, rate_cat) {
  paste(c(rate_cat, dim(max_index_mat), as.integer(!is.na(max_index_mat))),
    collapse = "_")
}

fits_in_model_space <- function(fits, max_index_mat, full_rate_category,
  require.entry = FALSE) {
  if(length(fits) == 0) return(TRUE)
  all(vapply(fits, function(fit) {
    if(is.null(fit)) return(TRUE)
    if(is.null(fit$index.mat) || length(fit$rate.cat) != 1 ||
      !is.finite(fit$rate.cat) || fit$rate.cat < 1 ||
      fit$rate.cat > full_rate_category) return(FALSE)
    allowed <- sub_max(max_index_mat, fit$rate.cat, full_rate_category)
    is.null(validate_index_mat(fit$index.mat, fit$rate.cat, require.entry,
      allowed))
  }, logical(1)))
}

checkpoint_is_compatible <- function(ckpt, max_index_mat, rate_category,
  require.entry = FALSE, criterion = "AIC") {
  if(is.null(ckpt) || !identical(ckpt$version, 4L) ||
    !identical(ckpt$criterion, criterion) ||
    !identical(ckpt$model_space_id,
      model_space_id(max_index_mat, rate_category)) ||
    (is.null(ckpt$current_fit) && !isTRUE(ckpt$finished))) return(FALSE)
  fits <- c(list(ckpt$current_fit, ckpt$best_fit), ckpt$accepted_models,
    ckpt$every_model, ckpt$score_cache,
    if(isTRUE(ckpt$finished)) ckpt$result$accepted_models else list())
  fits_in_model_space(fits, max_index_mat, rate_category, require.entry)
}

polish_due <- function(steps_since_best, steps_since_new, stall.limit,
  restart.interval) {
  steps_since_best >= stall.limit || steps_since_new >= restart.interval
}

annealing_temperature <- function(iteration, max.iterations, initial.temp,
  cooling.rate, temp.schedule, epoch.start = 1L) {
  final.temp <- max(0.001, initial.temp * cooling.rate ^ max.iterations)
  span <- max(1, max.iterations - epoch.start)
  progress <- min(1, max(0, (iteration - epoch.start) / span))
  if(temp.schedule == "exponential") {
    return(initial.temp * exp(log(final.temp / initial.temp) * progress))
  }
  if(temp.schedule == "linear") {
    return(initial.temp + (final.temp - initial.temp) * progress)
  }
  if(temp.schedule == "logarithmic") {
    scaled <- log1p(9 * progress) / log(10)
    return(initial.temp * exp(log(final.temp / initial.temp) * scaled))
  }
  stop("Unknown temperature schedule: ", temp.schedule, call. = FALSE)
}

polish_neighborhood <- function(initial_fit, neighbor_fn, score_fn, criterion,
  max.new.evaluations = Inf, is.cached = function(candidate) FALSE) {
  current <- initial_fit
  path <- list()
  evaluations <- 0L
  new_evaluations <- 0L
  deferred <- 0L
  repeat {
    proposals <- neighbor_fn(current)
    if(length(proposals) == 0) break
    fits <- lapply(proposals, function(candidate) {
      cached <- isTRUE(is.cached(candidate))
      if(!cached && new_evaluations >= max.new.evaluations) {
        deferred <<- deferred + 1L
        return(NULL)
      }
      if(!cached) new_evaluations <<- new_evaluations + 1L
      evaluations <<- evaluations + 1L
      score_fn(candidate)
    })
    keep <- !vapply(fits, is.null, logical(1))
    if(!any(keep)) break
    fits <- fits[keep]
    scores <- vapply(fits, function(x) x[[criterion]], numeric(1))
    better <- vapply(fits, dredge_fit_is_better, logical(1), current = current,
      criterion = criterion)
    if(!any(better)) break
    best_score <- min(scores[better])
    eligible <- which(better & scores <= best_score + 1e-8)
    complexity <- t(vapply(fits[eligible], dredge_fit_complexity, numeric(2)))
    best <- eligible[order(complexity[, 1], complexity[, 2], scores[eligible])[1]]
    current <- fits[[best]]
    path[[length(path) + 1L]] <- current
  }
  list(best_fit = current, path = path, evaluations = evaluations,
    new_evaluations = new_evaluations, deferred = deferred,
    complete = deferred == 0L)
}

dredge_fit_complexity <- function(fit) {
  c(max(fit$index.mat, na.rm = TRUE), sum(!is.na(fit$index.mat)))
}

dredge_inherited_rates <- function(current_fit, candidate_mat) {
  pars <- sort(unique(na.omit(as.vector(candidate_mat))))
  parent_rates <- unname(MatrixToPars(current_fit))
  parent_rates <- parent_rates[is.finite(parent_rates) & parent_rates > 0]
  if(length(parent_rates) == 0) return(NULL)
  fallback <- median(parent_rates)
  if(!identical(dim(current_fit$index.mat), dim(candidate_mat))) {
    return(rep(fallback, length(pars)))
  }
  rates <- vapply(pars, function(par) {
    values <- current_fit$solution[which(candidate_mat == par)]
    values <- unique(values[is.finite(values) & values > 0])
    if(length(values) == 0) fallback else mean(values)
  }, numeric(1))
  if(any(!is.finite(rates))) return(NULL)
  rates
}

dredge_fit_is_better <- function(candidate, current, criterion, tolerance = 1e-8) {
  candidate_score <- candidate[[criterion]]
  current_score <- current[[criterion]]
  if(!is.finite(candidate_score)) return(FALSE)
  if(!is.finite(current_score)) return(TRUE)
  if(candidate_score < current_score - tolerance) return(TRUE)
  if(abs(candidate_score - current_score) > tolerance) return(FALSE)
  candidate_complexity <- dredge_fit_complexity(candidate)
  current_complexity <- dredge_fit_complexity(current)
  candidate_complexity[1] < current_complexity[1] ||
    (candidate_complexity[1] == current_complexity[1] &&
      candidate_complexity[2] < current_complexity[2])
}

deterministic_neighbors <- function(current_fit, max_index_mat, full_rate_category,
  merge.threshold = 0, require.entry = FALSE) {
  index_mat <- current_fit$index.mat
  rc <- current_fit$rate.cat
  out <- list()
  ids <- canonical_id(index_mat, rc)
  add <- function(mat, candidate_rc, move_type, inherit = TRUE) {
    candidate_max <- sub_max(max_index_mat, candidate_rc, rc)
    if(!is.null(validate_index_mat(mat, candidate_rc, require.entry,
      candidate_max))) return()
    id <- canonical_id(mat, candidate_rc)
    if(id %in% ids) return()
    ids <<- c(ids, id)
    candidate <- list(new_index_mat = mat, rate_cat = candidate_rc,
      move_type = move_type)
    if(inherit) {
      candidate$initial_rates <- dredge_inherited_rates(current_fit, mat)
      candidate$source_scores <- current_fit[c("AIC", "AICc", "BIC")]
    }
    out[[length(out) + 1L]] <<- candidate
  }

  pars <- sort(unique(na.omit(as.vector(index_mat))))
  for(par in pars) {
    dropped <- prune_isolated_states(dropStateMatPars(index_mat, par), rc)
    candidate_rc <- if(is.null(dropped$rate_cat)) rc else dropped$rate_cat
    add(dropped$new_index_mat, candidate_rc, "polish_drop")
  }

  frequencies <- tabulate(index_mat, nbins = if(length(pars)) max(pars) else 0)
  tied <- !is.na(index_mat)
  tied[tied] <- frequencies[index_mat[tied]] > 1
  for(position in which(tied)) {
    reduced <- index_mat
    reduced[position] <- NA
    reduced <- prune_isolated_states(reduced, rc)
    candidate_rc <- if(is.null(reduced$rate_cat)) rc else reduced$rate_cat
    add(reduced$new_index_mat, candidate_rc, "polish_cell_drop", inherit=TRUE)
  }

  par_values <- vapply(pars, function(par)
    min(current_fit$solution[which(index_mat == par)]), numeric(1))
  floor_pars <- pars[par_values <= 1e-8]
  for(par in floor_pars) {
    dropped <- dropStateMatPars(index_mat, par)
    targets <- sort(unique(na.omit(as.vector(dropped))))
    positions <- which(is.na(dropped) & !is.na(max_index_mat))
    for(position in positions) {
      for(target in targets) {
        relocated <- dropped
        relocated[position] <- target
        add(renumber_index_mat(relocated), rc, "polish_relocate")
      }
    }
  }

  current_pars <- MatrixToPars(current_fit)
  eligible <- which(current_pars > 1e-8)
  groups <- list(seq_along(current_pars))
  if(rc > 1) {
    within <- lapply(seq_len(rc), function(i)
      grep(paste0("R", i, " .* -> R", i), names(current_pars)))
    groups <- c(within, list(setdiff(seq_along(current_pars), unlist(within))))
  }
  for(group in groups) {
    group <- intersect(group, eligible)
    if(length(group) < 2) next
    pairs <- combn(group, 2)
    distances <- abs(current_pars[pairs[1, ]] - current_pars[pairs[2, ]])
    chosen <- which(distances <= merge.threshold)
    if(length(chosen) == 0) chosen <- which(distances == min(distances))
    for(i in chosen) {
      merged <- renumber_index_mat(equateStateMatPars(index_mat,
        c(pairs[1, i], pairs[2, i])))
      add(merged, rc, "polish_merge")
    }
  }

  dropped <- is.na(index_mat) & !is.na(max_index_mat)
  free_positions <- which((tied | dropped) & !is.na(max_index_mat))
  for(position in free_positions) {
    freed <- index_mat
    freed[position] <- max(index_mat, na.rm = TRUE) + 1
    add(renumber_index_mat(freed), rc, "polish_free")
  }

  if(rc > 1) {
    n <- nrow(index_mat) / rc
    Q <- current_fit$solution
    Q[is.na(Q)] <- 0
    class_of <- rep(seq_len(rc), each = n)
    dead <- which(vapply(seq_len(rc), function(k)
      all(Q[class_of != k, class_of == k] == 0), logical(1)))
    for(gone in dead) {
      keep <- which(class_of != gone)
      collapsed <- renumber_index_mat(index_mat[keep, keep])
      rownames(collapsed) <- colnames(collapsed) <-
        paste0("(", rep(seq_len(n), rc - 1), ",R",
          rep(seq_len(rc - 1), each = n), ")")
      add(collapsed, rc - 1, "polish_collapse")
    }
  }

  out
}

root_connected <- function(index_mat) {
  reachable <- !is.na(index_mat)
  diag(reachable) <- TRUE
  for(k in seq_len(nrow(reachable))) {
    reachable <- reachable | outer(reachable[, k], reachable[k, ], `&`)
  }
  any(apply(reachable, 1, all))
}

# The one coherence check for a proposed model. Returns NULL when the matrix
# is a model, otherwise the reason it is not. A state with rates in only one
# direction is legitimate (e.g. an irreversible model with an absorbing
# state); set require.entry to forbid unenterable states.
validate_index_mat <- function(index_mat, rate_cat = 1, require.entry = FALSE,
  allowed_mat = NULL){
  if(is.null(index_mat)) return("no proposal")
  index_mat[index_mat == 0] <- NA
  if(length(dim(index_mat)) != 2 || nrow(index_mat) != ncol(index_mat)) return("not square")
  if(any(!is.na(diag(index_mat)))) return("diagonal parameter")
  if(!is.null(allowed_mat)) {
    allowed_mat[allowed_mat == 0] <- NA
    if(!identical(dim(index_mat), dim(allowed_mat))) return("dimensions do not match model space")
    if(any(!is.na(index_mat) & is.na(allowed_mat))) return("transition outside model space")
  }
  if(all(is.na(index_mat))) return("no free parameters")
  
  entered <- apply(!is.na(index_mat), 2, any)
  left <- apply(!is.na(index_mat), 1, any)
  if(any(!entered & !left)) return("isolated state")
  if(!root_connected(index_mat)) return("multiple source components")
  if(require.entry && any(!entered)) return("unenterable state")
  
  if(rate_cat > 1){
    if(nrow(index_mat) %% rate_cat != 0) return("dimensions do not match rate.cat")
    class_of <- rep(seq_len(rate_cat), each = nrow(index_mat) / rate_cat)
    if(any(!tapply(entered | left, class_of, any))) return("empty rate class")
  }
  NULL
}

is_valid_index_mat <- function(index_mat, require.entry = FALSE, allowed_mat = NULL){
  is.null(validate_index_mat(index_mat, require.entry = require.entry,
    allowed_mat = allowed_mat))
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
corHMMDredgeBase <- function(phy, data, rate.cat, root.p="maddfitz", tip.fog=NULL, fog.ip=0.01,
  pen.type="l1", lambda=1, rate.mat=NULL, node.states="marginal", fixed.nodes=FALSE, 
  ip=NULL, nstarts=0, n.cores=1,get.tip.states=FALSE, lewis.asc.bias=FALSE, collapse=FALSE, 
  lower.bound=1e-10, upper.bound=100, opts=NULL, p=NULL, use_RTMB=FALSE, prep=NULL) {

  if (!is.null(rate.mat)) {
    rate_reason <- validate_index_mat(rate.mat, rate.cat)
    if(!is.null(rate_reason)) {
      stop("rate.mat is not valid: ", rate_reason, call. = FALSE)
    }
  }

  ## --- use precomputed prep if available, otherwise compute it ---
  if (is.null(prep)) {
    prep <- corHMMDredgePrep(phy=phy, data=data, rate.cat=rate.cat, root.p=root.p,
      pen.type=pen.type, lambda=lambda, node.states=node.states,
      fixed.nodes=fixed.nodes, collapse=collapse, lower.bound=lower.bound,
      upper.bound=upper.bound, lewis.asc.bias=lewis.asc.bias,
      use_RTMB=use_RTMB, opts=opts, tip.fog=tip.fog, fog.ip=fog.ip)
  }
  
  ## unpack prep
  phy          <- prep$phy
  phy_original <- prep$phy_original
  data         <- prep$data
  input.data   <- prep$input.data
  data.legend  <- prep$data.legend
  CorData      <- prep$CorData
  nObs         <- prep$nObs
  H            <- prep$H
  lb           <- prep$lb
  ub           <- prep$ub
  lower.bound  <- prep$lower.bound
  upper.bound  <- prep$upper.bound
  mean.change  <- prep$mean.change
  StateNames   <- prep$StateNames
  root.p       <- prep$root.p
  node.states  <- prep$node.states
  lambda       <- prep$lambda
  opts         <- prep$opts
  nb.tip       <- prep$nb.tip
  nb.node      <- prep$nb.node
  use_RTMB     <- prep$use_RTMB
  collapse     <- prep$collapse
  pen.type     <- prep$pen.type
  lewis.asc.bias <- prep$lewis.asc.bias
  levels       <- prep$levels

  order.test <- FALSE
  model <- "ARD"
  
  ## --- fast rate.cat.set assembly using precomputed prep components ---
  if (!is.null(prep)) {
    nTraits <- prep$nTraits
    tmp     <- prep$tip_liks_template
    if (is.null(rate.mat)) {
      base_r <- prep$base_rate_mat
      if (rate.cat > 1) {
        StateMats <- replicate(rate.cat, base_r, simplify=FALSE)
        index.matrix <- getFullMat(StateMats)
      } else {
        index.matrix <- base_r
      }
      index.matrix[index.matrix == 0] <- NA
      index.matrix[!is.na(index.matrix)] <- seq_len(sum(!is.na(index.matrix)))
    } else {
      index.matrix <- rate.mat
      index.matrix[index.matrix == 0] <- NA
    }
    rate <- index.matrix
    rate[is.na(rate)] <- max(rate, na.rm=TRUE) + 1
    ## liks: just replicate tmp across rate.cat -- O(n) not O(n * full_setup)
    liks <- matrix(rep(tmp, rate.cat), nrow(tmp), nTraits * rate.cat)
    Q <- matrix(0, nrow(rate), ncol(rate))
    model.set.final <- list(
      rate.cat     = rate.cat,
      np           = max(rate) - 1,
      rate         = rate,
      index.matrix = index.matrix,
      liks         = liks,
      Q            = Q
    )
  } else {
    ## fallback: full computation when no prep available
    model.set.final <- rate.cat.set.corHMM.JDB(phy=phy, data=input.data,
      rate.cat=rate.cat, ntraits=nObs, model=model,
      rate.mat=rate.mat, collapse=collapse)
    if (is.null(rate.mat)) {
      model.set.final$index.matrix[!is.na(model.set.final$index.matrix)] <-
        seq_len(sum(!is.na(model.set.final$index.matrix)))
      model.set.final$rate <- model.set.final$index.matrix
      model.set.final$rate[is.na(model.set.final$rate)] <- max(model.set.final$rate, na.rm=TRUE) + 1
      model.set.final$np <- max(model.set.final$index.matrix, na.rm=TRUE)
    }
  }
  
  ## --- apply tip fog (after liks is built, since liks depends on rate.cat) ---
  set.fog <- prep$set.fog
  fog.vec <- prep$fog.vec
  fog.est <- prep$fog.est   ## non-NULL only in fixed mode
  
  if (!is.null(fog.est)) {
    ## fixed fog: bake into liks once, no extra parameters
    tip.fog_expanded <- if (rate.cat > 1) rep(fog.est, rate.cat) else fog.est
    for (tip.index in seq_len(nb.tip)) {
      num.zeros <- sum(liks[tip.index, ] == 0)
      if (num.zeros > 0) {
        if (rate.cat > 1) {
          liks[tip.index, which(liks[tip.index,] == 1)] <-
            1 - (sum(tip.fog_expanded[which(liks[tip.index,] != 1)]) / rate.cat)
        } else {
          liks[tip.index, which(liks[tip.index,] == 1)] <-
            1 - sum(tip.fog_expanded[which(liks[tip.index,] != 1)])
        }
        liks[tip.index, which(liks[tip.index,] == 0)] <-
          tip.fog_expanded[which(liks[tip.index,] == 0)]
      }
    }
    model.set.final$liks <- liks
    fog.est_out <- fog.est   ## pass through to output
  } else if (set.fog) {
    ## estimated fog: expand fog.vec for rate.cat and attach to model.set.final
    model.set.final$fog.vec <- if (rate.cat > 1) rep(fog.vec, rate.cat) else fog.vec
    fog.est_out <- NULL   ## estimated, not known until after optimization
  } else {
    fog.est_out <- NULL
  }
  
  ## --- np_tot accounts for fog parameters if estimated ---
  n_fog_pars <- if (set.fog) length(unique(model.set.final$fog.vec)) else 0L
  np_tot <- model.set.final$np + n_fog_pars
  
  
  lower <- rep(lb, model.set.final$np)
  upper <- rep(ub, model.set.final$np)
  
  ## --- StateNames with rate categories ---
  if (collapse) {
    StateNames_rc <- rep(gsub("_", "|", CorData$ObservedTraits), rate.cat)
    RCNames <- rep(paste("R", 1:rate.cat, sep=""), each=length(CorData$ObservedTraits))
  } else {
    StateNames_rc <- rep(gsub("_", "|", CorData$PossibleTraits), rate.cat)
    RCNames <- rep(paste("R", 1:rate.cat, sep=""), each=length(CorData$PossibleTraits))
  }
  if (rate.cat > 1) StateNames_rc <- paste(RCNames, StateNames_rc)
  
  ## --- optimizer setup ---
  if (use_RTMB) {
    RTMB_obj <- mkdev.corhmm_rtmb(rep(0, np_tot), phy,
      liks=model.set.final$liks, Q=model.set.final$Q,
      rate=model.set.final$rate, root.p=root.p,
      rate.cat = rate.cat,
      ## these two are intentionally disabled
      order.test = FALSE,
      lewis.asc.bias = FALSE,
      set.fog = set.fog,
      fog.vec = model.set.final$fog.vec,
      pen.type=pen.type, lambda=lambda)
    
    devfun      <- function(p, ...) RTMB_obj$fn(p)
    devfun_grad <- function(p, ...) RTMB_obj$gr(p)
  }
  
  ## --- optimization ---
  if (!is.null(p)) {
    cat("Calculating likelihood from a set of fixed parameters\n")
    out <- NULL
    p <- p * H
    est.pars <- log(p)
    out$objective <- dev.corhmm.dredge(est.pars, phy=phy,
      liks=model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate,
      root.p=root.p, rate.cat=rate.cat, order.test=order.test,
      lewis.asc.bias=lewis.asc.bias, pen.type=pen.type, lambda=lambda)
    loglik   <- -out$objective
    est.pars <- exp(est.pars)
    
  } else if (is.null(ip)) {
    
    random.restart <- function(nstarts) {
      starts <- if (mean.change == 0) {
        rep(0.01 + exp(lb), model.set.final$np)
      } else {
        sort(rexp(model.set.final$np, 1/mean.change), decreasing=TRUE)
      }
      starts[starts < exp(lb)] <- exp(lb)
      starts[starts > exp(ub)] <- exp(lb)
      
      lower_r <- lower
      upper_r <- upper
      if (set.fog) {
        starts  <- c(rep(prep$fog.ip, n_fog_pars), starts)
        lower_r <- c(rep(lb, n_fog_pars), lower_r)
        upper_r <- c(rep(log(0.50), n_fog_pars), upper_r)
        tmp <- matrix(, 1, ncol=(1 + np_tot))
      } else {
        tmp <- matrix(, 1, ncol=(1 + model.set.final$np))
      }
      
      if (use_RTMB) {
        out <- tryCatch(
          nlminb(log(starts), objective=devfun, gradient=devfun_grad,
            lower=lower_r, upper=upper_r),
          error   = function(e) NULL,
          warning = function(w) suppressWarnings(
            nlminb(log(starts), objective=devfun, lower=lower_r, upper=upper_r))
        )
        if (is.null(out)) return(NULL)
        tmp[, 1] <- out$objective
        tmp[, 2:(np_tot + 1)] <- out$par
      } else {
        out <- nloptr(x0=log(starts), eval_f=dev.corhmm.dredge,
          lb=lower_r, ub=upper_r, opts=opts, phy=phy,
          liks=model.set.final$liks, Q=model.set.final$Q,
          rate=model.set.final$rate, root.p=root.p, rate.cat=rate.cat,
          order.test=order.test, lewis.asc.bias=lewis.asc.bias,
          pen.type=pen.type, lambda=lambda)
        tmp[, 1] <- out$objective
        tmp[, 2:(model.set.final$np + 1)] <- out$solution
      }
      tmp
    }    
    restart.set <- if (n.cores > 1) {
      mclapply(1:nstarts, random.restart, mc.cores=n.cores)
    } else {
      lapply(1:nstarts, random.restart)
    }
    restart.set <- Filter(Negate(is.null), restart.set)  ## drop failed restarts
    if (length(restart.set) == 0) return(NULL)           ## all failed, bail out
    best.fit <- which.min(unlist(lapply(restart.set, function(x) x[1])))
    out <- list(
      objective = unlist(restart.set[[best.fit]][, 1]),
      solution  = unlist(restart.set[[best.fit]][, 2:(model.set.final$np + 1)])
    )
    loglik   <- -out$objective
    est.pars <- exp(out$solution)
    if (set.fog) {
      fog.est_out <- est.pars[seq_len(n_fog_pars)]
      est.pars    <- est.pars[-seq_len(n_fog_pars)]
    }
    
  } else {
    
    if (use_RTMB) {
      out <- nlminb(rep(log(ip), length.out=model.set.final$np),
        objective=devfun, gradient=devfun_grad,
        lower=lower, upper=upper)
      loglik   <- -out$objective
      est.pars <- exp(if (use_RTMB) out$par else out$solution)
      if (set.fog) {
        fog.est_out <- est.pars[seq_len(n_fog_pars)]
        est.pars    <- est.pars[-seq_len(n_fog_pars)]
      }
    } else {
      out <- nloptr(x0=rep(log(ip), length.out=model.set.final$np),
        eval_f=dev.corhmm.dredge, lb=lower, ub=upper, opts=opts,
        phy=phy, liks=model.set.final$liks, Q=model.set.final$Q,
        rate=model.set.final$rate, root.p=root.p, rate.cat=rate.cat,
        order.test=order.test, lewis.asc.bias=lewis.asc.bias,
        pen.type=pen.type, lambda=lambda)
      loglik   <- -out$objective
      est.pars <- exp(if (use_RTMB) out$par else out$solution)
      if (set.fog) {
        fog.est_out <- est.pars[seq_len(n_fog_pars)]
        est.pars    <- est.pars[-seq_len(n_fog_pars)]
      }
    }
  }
  
  ## --- ancestral state reconstruction ---
  TIPS <- 1:nb.tip
  if (node.states %in% c("marginal", "scaled")) {
    lik.anc <- ancRECON(phy, input.data, est.pars, rate.cat,
      rate.mat=model.set.final$rate, method=node.states, ntraits=NULL,
      root.p=root.p, model=model, get.tip.states=get.tip.states,
      tip.fog=fog.est_out, collapse=collapse)
    pr <- apply(lik.anc$lik.anc.states, 1, which.max)
    phy$node.label <- pr
    tip.states <- lik.anc$lik.tip.states
    row.names(tip.states) <- phy$tip.label
  } else if (node.states == "joint") {
    lik.anc <- ancRECON(phy, input.data, est.pars, rate.cat,
      rate.mat=model.set.final$rate, method=node.states, ntraits=NULL,
      root.p=root.p, model=model, get.tip.states=get.tip.states, collapse=collapse)
    phy$node.label <- lik.anc$lik.anc.states
    tip.states <- lik.anc$lik.tip.states
  } else {
    lik.anc <- list(lik.tip.states=NA, lik.anc.states=NA, info.anc.states=NA)
    phy$node.label <- NA
    tip.states <- NA
  }
  
  ## --- finalize output ---
  solution <- matrix(est.pars[model.set.final$index.matrix], dim(model.set.final$index.matrix))
  np <- model.set.final$np
  index.matrix <- model.set.final$index.matrix
  try({rownames(solution) <- colnames(solution) <- StateNames_rc})
  try({rownames(index.matrix) <- colnames(index.matrix) <- StateNames_rc})
  solution <- solution / H
  solution[solution < lower.bound] <- lower.bound
  
  np_for_aic <- model.set.final$np + n_fog_pars
  fit_criteria <- information_criteria(loglik, np_for_aic, nb.tip)
  AIC <- fit_criteria$AIC
  AICc <- fit_criteria$AICc
  BIC <- fit_criteria$BIC
  
  if (is.character(node.states) && node.states %in% c("marginal", "scaled")) {
    colnames(lik.anc$lik.anc.states) <- StateNames_rc
  }
  
  if (loglik == -1e+06) {
    warning("corHMM may have failed to optimize correctly, consider checking inputs and running again.",
      immediate. = TRUE)
  }
  
  tip.fog.probs = if (set.fog && !is.null(fog.est_out)) {
    probs <- numeric(length(model.set.final$fog.vec))
    probs[] <- c(fog.est_out, 0)[model.set.final$fog.vec]
    names(probs) <- StateNames_rc
    probs
  } else {
    NULL
  }
  
  if(use_RTMB){
      rtmb_obj <- RTMB_obj
    }else{
      rtmb_obj <- NULL
      }
  
  obj <- list(
    loglik      = loglik,
    AIC         = AIC,
    AICc        = AICc,
    BIC         = BIC,
    rate.cat    = rate.cat,
    solution    = solution,
    index.mat   = index.matrix,
    data        = input.data,
    data.legend = data.legend,
    phy         = phy_original,
    states      = lik.anc$lik.anc.states,
    tip.states  = tip.states,
    states.info = lik.anc$info.anc.states,
    iterations  = out$iterations,
    collapse    = collapse,
    root.p      = root.p,
    pen.type    = pen.type,
    lambda      = lambda,
    lower.bound = lower.bound,
    use_RTMB    = use_RTMB,
    tip.fog     = tip.fog,
    fog.ip      = fog.ip,
    tip.fog.p   = tip.fog.probs,
    rtmb_obj    = rtmb_obj
  )
  class(obj) <- "corhmm"
  return(obj)
}

## Compute everything that doesn't depend on rate.mat
corHMMDredgePrep <- function(phy, data, rate.cat, root.p="maddfitz", tip.fog=NULL, fog.ip=0.01,
  pen.type="l1", lambda=1, node.states="marginal", fixed.nodes=FALSE, collapse=FALSE,
  lower.bound=1e-10, upper.bound=100, lewis.asc.bias=FALSE, use_RTMB=FALSE, opts=NULL) {
  
  ## --- node.states validation ---
  if (is.null(node.states)) {
    stop("No model for ancestral states selected. Please pass one of: joint, marginal, scaled, none.")
  }
  valid.models <- c("joint", "marginal", "scaled", "none")
  if (!any(valid.models == node.states)) {
    stop(paste0("'", node.states, "' is not a valid node.states option."))
  }
  if (length(node.states) > 1) {
    node.states <- "marginal"
    cat("No model selected for 'node.states'. Will perform marginal ancestral state estimation.\n")
  }
  
  ## --- fixed.nodes check ---
  if (fixed.nodes == FALSE && !is.null(phy$node.label)) {
    phy$node.label <- NULL
    cat("You specified 'fixed.nodes=FALSE' but included node labels. These have been removed.\n")
  }
  
  ## --- root.p normalization ---
  if (!is.null(root.p) && !is.character(root.p)) {
    root.p <- root.p / sum(root.p)
  }
  
  ## --- pen.type ---
  if (pen.type == "unreg") lambda <- 0
  
  ## --- rescale phy ---
  phy_original <- phy
  H <- max(node.depth.edgelength(phy))
  phy$edge.length <- phy$edge.length / H
  upper.bound <- upper.bound * H
  lower.bound <- lower.bound / H
  lb <- log(lower.bound)
  ub <- log(upper.bound)
  
  ## --- process data ---
  input.data <- data
  CorData <- corProcessData(data, collapse = collapse)
  data.legend <- data <- CorData$corData
  nObs <- length(CorData$ObservedTraits)
  
  ## --- match tree and data ---
  matching <- match.tree.data(phy, data)
  data <- matching$data
  phy  <- matching$phy
  
  ## --- invariant character checks ---
  if (nlevels(as.factor(data[, 1])) <= 1) {
    stop("Character is invariant. Analysis stopped.")
  }
  lvls <- as.factor(data[, 1])
  if (nlevels(lvls) == 2 && length(which(lvls == "?"))) {
    stop("Character is invariant. Analysis stopped.")
  }
  
  ## --- zero branch lengths ---
  if (any(phy$edge.length <= .Machine$double.eps)) {
    warning(paste0("Branch lengths of 0 detected. Adding ", sqrt(.Machine$double.eps)),
      immediate. = TRUE)
    phy$edge.length <- phy$edge.length + sqrt(.Machine$double.eps)
  }
  
  ## --- data.sort for parsimony starting values ---
  data.sort <- data.frame(data[, 2], data[, 2], row.names = data[, 1])
  data.sort <- data.sort[phy$tip.label, ]
  counts <- table(data.sort[, 1])
  levels <- levels(as.factor(data.sort[, 1]))
  
  ## --- tree reordering ---
  phy <- reorder(phy, "pruningwise")
  
  ## --- state names ---
  if (collapse) {
    StateNames <- gsub("_", "|", CorData$ObservedTraits)
  } else {
    StateNames <- gsub("_", "|", CorData$PossibleTraits)
  }
  
  ## --- parsimony-based mean.change for starting values ---
  taxa.missing.data.drop <- which(is.na(data.sort[, 1]))
  if (length(taxa.missing.data.drop) != 0) {
    dat <- as.matrix(data.sort)
    dat.red <- dat[-taxa.missing.data.drop, ]
    phy.red <- drop.tip(phy, taxa.missing.data.drop)
    dat.red <- phyDat(dat.red, type = "USER", levels = levels)
    phy.tmp <- multi2di(phy.red)
    par.score <- parsimony(phy.tmp, dat.red, method = "fitch") / 2
  } else {
    dat <- as.matrix(data.sort)
    dat <- phyDat(dat, type = "USER", levels = levels)
    phy.tmp <- multi2di(phy)
    par.score <- parsimony(phy.tmp, dat, method = "fitch") / 2
  }
  tl <- sum(phy$edge.length)
  mean.change <- par.score / tl
  
  ## --- opts default ---
  if (is.null(opts)) {
    if (use_RTMB) {
      opts <- list("algorithm" = "NLOPT_LD_MMA", "maxeval" = "1000000",
        "ftol_rel" = .Machine$double.eps^0.5)
    } else {
      opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = "1000000",
        "ftol_rel" = .Machine$double.eps^0.5)
    }
  }
  
  ## --- precompute the parts of rate.cat.set that don't depend on rate.mat ---
  ## base single-rate-cat state matrix (used to build liks and as fallback)
  base_rate_mat <- getStateMat4Dat(input.data, model="ARD", collapse=collapse)$rate.mat
  nTraits <- nrow(base_rate_mat)
  
  ## tip state likelihood template (nb.tip + nb.node) x nTraits
  ## only depends on matched data and nTraits -- rate.mat invariant
  matching_for_liks <- match.tree.data(phy, data)  ## phy already rescaled+reordered here
  nb.tip_prep  <- length(phy$tip.label)
  nb.node_prep <- phy$Nnode
  tmp <- matrix(0, nb.tip_prep + nb.node_prep, nTraits)
  for (i in seq_len(nb.tip_prep)) {
    focal_state <- matching_for_liks$data[i, 2]
    if (focal_state == "?") {
      tmp[i, ] <- 1
    } else {
      state_index <- as.numeric(unlist(strsplit(as.character(focal_state), "&")))
      tmp[i, state_index] <- 1
    }
  }
  
  ## --- tip fog preprocessing ---
  ## We can determine fog mode here; fixed fog modifies liks (static),
  ## estimated fog just sets up fog.vec for use in corHMMDredgeBase
  set.fog  <- FALSE
  fog.vec  <- NULL
  fog.est  <- NULL   ## will hold fixed fog values if fixed mode
  
  if (!is.null(tip.fog)) {
    if (sum(tip.fog) < 1) {
      ## fixed fog mode: expand to cover all states if scalar
      if (length(tip.fog) == 1) {
        tip.fog <- rep(tip.fog, nTraits)
      }
      ## store for application in corHMMDredgeBase after liks is built
      ## (liks depends on rate.cat which isn't known in prep)
      fog.est <- tip.fog
      set.fog <- FALSE
    } else {
      ## estimated fog mode: tip.fog is used as fog.vec index
      fog.vec <- tip.fog
      set.fog <- TRUE
    }
  }
  
  list(
    phy          = phy,
    phy_original = phy_original,
    data         = data,
    input.data   = input.data,
    data.legend  = data.legend,
    data.sort    = data.sort,
    CorData      = CorData,
    nObs         = nObs,
    counts       = counts,
    levels       = levels,
    H            = H,
    lb           = lb,
    ub           = ub,
    lower.bound  = lower.bound,
    upper.bound  = upper.bound,
    mean.change  = mean.change,
    StateNames   = StateNames,
    root.p       = root.p,
    node.states  = node.states,
    lambda       = lambda,
    opts         = opts,
    nb.tip       = length(phy$tip.label),
    nb.node      = phy$Nnode,
    use_RTMB     = use_RTMB,
    collapse     = collapse,
    pen.type     = pen.type,
    lewis.asc.bias = lewis.asc.bias,
    base_rate_mat = base_rate_mat,
    nTraits       = nTraits,
    tip_liks_template = tmp,   ## the nb.tip+nb.node x nTraits matrix
    tip.fog  = tip.fog,
    fog.ip   = fog.ip,
    fog.est  = fog.est,   ## fixed fog values (NULL if estimated or no fog)
    fog.vec  = fog.vec,   ## fog index vector (NULL if fixed or no fog)
    set.fog  = set.fog
    
  )
}

### The function used to optimize parameters:

dev.corhmm.dredge <- function(p,phy,liks,Q,rate,root.p,rate.cat,order.test,lewis.asc.bias,pen.type="l1",lambda=1){
  p = exp(p)
  cp_root.p <- root.p
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  TIPS <- 1:nb.tip
  comp <- numeric(nb.tip + nb.node)
  #Obtain an object of all the unique ancestors
  edge1 <- phy$edge[,1]
  edge2 <- phy$edge[,2]
  edge.length <- phy$edge.length
  anc <- unique(edge1)
  #Map every edge to its ancestor in one pass instead of rescanning phy$edge per node
  desRowsList <- getDesRows(edge1, anc)
  k.rates <- dim(Q)[2] / 2
  if (any(is.nan(p)) || any(is.infinite(p))) return(1000000)
  
  Q[] <- c(p, 0)[rate]
  diag(Q) <- -rowSums(Q)
  #Q is fixed for the rest of this call, so decompose it once and reuse it for
  #every branch rather than exponentiating per edge:
  Pv <- makeExpmFuns(Q)$Pv
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
    desRows <- desRowsList[[i]]
    desNodes <- edge2[desRows]
    v <- 1
    #Loops through all descendants of focal (how we deal with polytomies):
    for (desIndex in seq_along(desRows)){
      v <- v*Pv(edge.length[desRows[desIndex]], liks[desNodes[desIndex],])
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
  node.loglik <- sum(log(comp[-TIPS]))
  if (is.na(node.loglik)){return(1000000)}
  liks.root <- liks[root,]

  if (is.null(root.p)){
    equil.root <- getEquilRoot(Q)
    flat.root = equil.root
    k.rates <- 1/length(which(!is.na(equil.root)))
    flat.root[!is.na(flat.root)] = k.rates
    flat.root[is.na(flat.root)] = 0
    loglik<- -(node.loglik + log(sum(flat.root * liks.root)))
  }
  if(is.character(root.p)){
    # root.p==yang will fix root probabilities based on the inferred rates: q10/(q01+q10)
    if(root.p == "yang"){
      #root.test was already computed above from the same Q
      root.p <- c(root.test/sum(root.test))
      loglik <- -(node.loglik + log(sum(root.p * liks.root)))
      if(is.infinite(loglik)){
        return(1000000)
      }
    }else{
      # root.p==maddfitz will fix root probabilities according to FitzJohn et al 2009 Eq. 10:
      root.p = liks.root / sum(liks.root)
      loglik <- -(node.loglik + log(sum(exp(log(root.p)+log(liks.root)))))
    }
  }else{
    if(is.numeric(root.p[1])){
      loglik <- -(node.loglik + log(sum(exp(log(root.p)+log(liks.root)))))
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
  return(loglik + (pen_score * lambda))
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
                                         use_RTMB = corhmm_obj$use_RTMB,
                                         fog.ip = corhmm_obj$fog.ip,
                                         tip.fog = corhmm_obj$tip.fog)
      
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

#' @export
#' @method print corhmm.dredge
print.corhmm.dredge <- function(x, ...) {
  sa_fits <- attr(x, "dredge_history")
  criterion <- attr(x, "criterion")
  if(is.null(criterion)) criterion <- if(!is.null(sa_fits[[1]]$criterion))
    sa_fits[[1]]$criterion else "AIC"
  
  # Header
  n_rc <- length(sa_fits)
  rate_cats <- sapply(sa_fits, "[[", "rate_category")
  cat("corhmm.dredge object\n")
  cat("  Models retained :", length(x), "\n")
  cat("  Rate categories :", paste(rate_cats, collapse = ", "), "\n")
  
  # SA diagnostics if available
  if (!is.null(sa_fits)) {
    for (i in seq_along(sa_fits)) {
      s <- sa_fits[[i]]
      cat(sprintf("  [Rate cat %d] Iterations: %d | Acceptance rate: %.1f%% | Restarts: %d | Polish: %d | Best %s: %.3f\n",
        rate_cats[i],
        s$iterations,
        s$acceptance_rate * 100,
        s$restart_count,
        if(is.null(s$polish_runs)) 0L else s$polish_runs,
        if(is.null(s$criterion)) criterion else s$criterion,
        s$best_score))
      cat(sprintf("               Fits: %d | Stopped: %s\n",
        s$unique_structures, s$stop_reason))
    }
  }
  
  cat("\n")
  
  # Model table
  tbl <- getModelTable(x, type = criterion)
  print(tbl)
  
  # Why proposals died
  rej <- unlist(lapply(sa_fits, "[[", "rejections"))
  if (length(rej) > 0) {
    cat("\nRejected proposals:\n")
    print(sort(tapply(rej, names(rej), sum), decreasing = TRUE))
  }
  
  # Flat directions
  flat <- attr(x, "flat_directions")
  if (!is.null(flat)) {
    cat("\nFlat directions (same lnLik, different np):\n")
    print(flat)
  }
  
  best_idx <- which.min(tbl[[paste0("d", criterion)]])
  best <- x[[best_idx]]
  floor_bound <- if (!is.null(best$lower.bound)) best$lower.bound else NA_real_
  floored <- !is.na(best$solution) & !is.na(floor_bound) &
    best$solution <= floor_bound
  cat("\n--- Best model ---\n")
  print(best)
  if (any(floored)) {
    cat(sum(floored), "free rate(s) are at the optimizer lower bound; they are not structural zeros\n")
  }
  
  # Footer hint
  cat("\nAccess SA trace: attr(x, 'dredge_history')\n")
  
  invisible(x)
}

plotDredgeTrace <- function(dredge_fits,
  break_size = 5,
  palette = c("drop" = "#A23B72",
    "merge" = "#2E86AB",
    "free" = "#F18F01",
    "restart" = "#7209B7",
    "lump" = "#C73E1D",
    "collapse" = "#3B7A57",
    "none" = "grey60"),
  legend = TRUE,
  legend.pos = "topright",
  ...) {
  
  # 1. Input Validation
  if (!is.list(dredge_fits)) {
    stop("Input must be a corHMMDredge object.")
  }
  sa_fits <- attr(dredge_fits, "dredge_history")
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
    
    # Filter out NAs which can occur if the loop ends prematurely and can happen when the model is not possible to be fit
    valid_indices <- !is.na(fit_data$every_move) & 
      !is.na(fit_data$accepted) &
      !is.infinite(fit_data$every_score)
    
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
