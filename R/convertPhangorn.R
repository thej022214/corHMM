#' @title Convert phangorn reconstruction to a vector
#'
#' @description
#' Converts a character reconstruction from phangorn into a vector
#'
#' @param x The phyDat object that contains a character reconstruction from phangorn
#' @param site The character number to convert into a vector
#' @param best A logical indicating whether the state that maximizes some function (likelihood, parsimony, etc.) is to be returned.
#'
#' @details
#' Creates a vector that contains the best tips and node state from a phangorn reconstruction.
ConvertPhangornReconstructions <- function(x, site=1, best=TRUE) {
   x.local <- subset(x, , site) 
   nc = attr(x.local, "nc")
   y = matrix(unlist(x.local[]), ncol = nc, byrow = TRUE)
   rownames(y) <- names(x.local[])
   colnames(y) <- attr(x.local, "levels")
   result <- y
   if(best) {
       best.vector <- rep(NA,nrow(result))
       for (i in sequence(nrow(result))) {
          best.vector[i] <- sample(which.max(result[i,]), 1) #so resolve ties randomly
       }
       names(best.vector) <-rownames(y)
       result<-best.vector
   }
   return(result)
}


