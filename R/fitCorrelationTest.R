# a function for silencing models
silence <- function(x){
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

fitCorrelationTest <- function(phy, data, simplified_models=FALSE){
  cat("Begining test of correlation...\n")
  rate_cat_mat <- getRateCatMat(2)
  indep_model_1 <- getStateMat4Dat(data, "ARD", collapse = FALSE, indep = TRUE)$rate.mat
  indep_model_2 <- getFullMat(list(indep_model_1, indep_model_1), RateClassMat = rate_cat_mat)
  corr_model_1 <- getStateMat4Dat(data, "ARD", collapse = FALSE, indep = FALSE)$rate.mat
  corr_model_2 <- getFullMat(list(corr_model_1, corr_model_1), RateClassMat = rate_cat_mat)
  cat("\nFitting an independent model...\n")
  independent_model_fit <- silence(corHMM(phy = phy, data = data, rate.cat = 1, rate.mat = indep_model_1))
  cat("Fitting a hidden Markov independent model...\n")
  hidden_independent_model_fit <- silence(corHMM(phy = phy, data = data, rate.cat = 2, rate.mat = indep_model_2))
  cat("Fitting a correlated model...\n")
  correlated_model_fit <- silence(corHMM(phy = phy, data = data, rate.cat = 1, rate.mat = corr_model_1))
  cat("Fitting a hidden Markov correlated model...\n")
  hidden_correlated_model_fit <- silence(corHMM(phy = phy, data = data, rate.cat = 2, rate.mat = corr_model_2))
  model_list <- list(independent_model_fit = independent_model_fit, hidden_Markov_independent_model_fit = hidden_independent_model_fit, correlated_model_fit = correlated_model_fit, hidden_Markov_correlated_model_fit = hidden_correlated_model_fit)
  if(simplified_models){
    simp_indep_model_1 <- equateStateMatPars(indep_model_1, list(c(1,3), c(2,4)))
    simp_indep_model_2 <- getFullMat(list(simp_indep_model_1, simp_indep_model_1), rate_cat_mat)
    simp_corr_model_1 <- equateStateMatPars(corr_model_1, list(c(1,3),c(2,5), c(4,7), c(6,8)))
    simp_corr_model_2 <- getFullMat(list(simp_corr_model_1, simp_corr_model_1), rate_cat_mat)
    cat("\nFitting an simplified independent model...\n")
    simplified_independent_model_fit <- silence(corHMM(phy = phy, data = data, rate.cat = 1, rate.mat = simp_indep_model_1))
    cat("Fitting a simplified hidden Markov independent model...\n")
    simplified_hidden_independent_model_fit <- silence(corHMM(phy = phy, data = data, rate.cat = 2, rate.mat = simp_indep_model_2))
    cat("Fitting a simplified correlated model...\n")
    simplified_correlated_model_fit <- silence(corHMM(phy = phy, data = data, rate.cat = 1, rate.mat = simp_corr_model_1))
    cat("Fitting a simplified hidden Markov correlated model...\n")
    simplified_hidden_correlated_model_fit <- silence(corHMM(phy = phy, data = data, rate.cat = 2, rate.mat = simp_corr_model_2))
    simplified_model_list <- list(simplified_independent_model_fit = simplified_independent_model_fit, simplified_hidden_Markov_independent_model_fit = simplified_hidden_independent_model_fit, simplified_correlated_model_fit = simplified_correlated_model_fit, simplified_hidden_Markov_correlated_model_fit = simplified_hidden_correlated_model_fit)
    model_list <- c(simplified_model_list, model_list)
  }
  cat("Done.\n")
  class(model_list) <- "corhmm_list"
  return(model_list)
}

getModelTable <- function(model_list, type="AIC"){
  # checks
  ParCount <- unlist(lapply(model_list, function(x) max(x$index.mat, na.rm = TRUE)))
  nTip <- length(model_list[[1]]$phy$tip.label)
  AIC <- simplify2array(lapply(model_list, "[[", type))
  dAIC <- AIC - min(AIC)
  AICwt <- exp(-0.5 * dAIC)/sum(exp(-0.5 * dAIC))
  LogLik <- simplify2array(lapply(model_list, "[[", "loglik"))
  out <- data.frame(np = ParCount, lnLik = LogLik, AIC = AIC, dAIC = dAIC, AICwt = AICwt)
  colnames(out) <- gsub("AIC", type, colnames(out))
  return(out)
}


print.corhmm_list <- function(x, ...){
  model_table <- getModelTable(x)
  model_names <- gsub("_", " ", gsub("_fit", "", names(x)))
  if(length(which(model_table$dAIC < 2)) > 1){
    my_message <- paste0("Multiple models fit well to the data (dAIC below 2):\n", paste(model_names[which(model_table$dAIC < 2)], collapse = "\n"))
  }else{
    my_message <- paste0("The ", model_names[which.min(model_table$dAIC)], " was the best fit to the data.\n")
  }
  cat(my_message)
  cat("\n")
  print(round(model_table, 2))
  cat("\n")
}


### testing 
# require(corHMM)
# 
# data(primates)
# phy <- multi2di(primates[[1]])
# data <- primates[[2]]
# 
# test <- fitCorrelationTest(phy, data, TRUE)
# test

