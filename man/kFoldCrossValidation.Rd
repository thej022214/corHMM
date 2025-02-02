\name{kFoldCrossValidation}
\alias{kFoldCrossValidation}
\title{Perform K-Fold Cross-Validation for corHMM Models}
\description{
  This function performs k-fold cross-validation on a given \code{corHMM} model by dividing the data into k equally sized subsets. The function evaluates model performance across multiple lambda regularization values, if provided. Optionally, it can save the trained models for each fold and return the cross-validation results.
}
\usage{
kFoldCrossValidation(corhmm_obj, k, lambdas = NULL, return_model = TRUE, 
save_model_dir = NULL, model_name = NULL)
}
\arguments{
  \item{corhmm_obj}{A \code{corHMM} object that contains a fitted model.}
  \item{k}{An integer specifying the number of folds to divide the data into for cross-validation.}
  \item{lambdas}{A numeric vector of lambda regularization values to evaluate during cross-validation. If \code{NULL}, the lambda value from \code{corhmm_obj} will be used. Defaults to \code{NULL}.}
  \item{return_model}{A logical value indicating whether to return the trained models for each fold. Defaults to \code{TRUE}.}
  \item{save_model_dir}{A character string specifying the directory to save the trained models for each fold. If \code{NULL}, models will not be saved. Defaults to \code{NULL}.}
  \item{model_name}{A character string specifying the base name for saved model files. If \code{NULL}, a default name \code{"corhmm.obj"} is used. Defaults to \code{NULL}.}
}
\details{
  The function splits the data into \code{k} folds and trains a separate \code{corHMM} model for each fold by leaving one fold out as the test set. The remaining folds are used for training the model. The performance of the model is evaluated on the test set using a divergence-based (Jensen-Shannon Divergence) scoring method. Evaluations are based on estimating the tips which were removed for that particular fold given the newly fitted model.

  The function supports evaluating models across different lambda regularization values. If \code{lambdas} are provided, models are trained and evaluated for each lambda value. The results, including the models (if \code{return_model = TRUE}) and cross-validation scores, are returned as a list.
}
\value{
  A list of cross-validation results, including the following components:
  \item{models}{A list of the trained models for each fold (if \code{return_model = TRUE}).}
  \item{scores}{A numeric vector of the cross-validation scores for each fold.}
  \item{averageScore}{The average cross-validation score across all folds.}
}
\examples{
\donttest{
data(primates)
phy <- multi2di(primates[[1]])
data <- primates[[2]]
dredge_fits <- corHMMDredge(phy = phy, data = data, max.rate.cat = 1, pen.type = "l1", 
	root.p = "maddfitz", lambda = 1, nstarts = 10, n.cores = 10)
model_table <- getModelTable(dredge_fits)
dredge_model <- dredge_fits[[which.min(model_table$dAIC)]]
k_fold_res <- kFoldCrossValidation(dredge_model, k = 5, lambdas = c(0,0.25,0.5,0.75,1))
getCVTable(k_fold_res)
}
}
\author{
  James D. Boyko
}
