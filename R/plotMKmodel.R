getRateMatCoord <- function(RateMat){
  x.vec <- seq(from = 2, to = 4, length.out = dim(RateMat)[1]+1)
  y.vec <- seq(from = 1, to = -1, length.out = dim(RateMat)[1]+1)
  x.mat <- matrix(x.vec, length(x.vec), length(x.vec), byrow = TRUE)
  y.mat <- matrix(y.vec, length(y.vec), length(y.vec))
  diag.coord <- cbind(diag(x.mat), diag(y.mat))
  diag(x.mat) <- NA
  diag(y.mat) <- NA
  mat.coord <- cbind(as.vector(na.exclude(as.vector(x.mat))), as.vector(na.exclude(as.vector(y.mat))))
  Ind <- mat.coord[,1] == 2 | mat.coord[,2] == 1
  name.coord <- mat.coord[Ind,]
  mat.coord <- mat.coord[!Ind,]
  
  return(list(MatrixCoord = mat.coord, NameCoord = name.coord, DiagCoord = diag.coord))
}

getRateMats <- function(pp=NULL, rate.mat=NULL, rate.cat=NULL){
  if(!is.null(pp)){
    rate.cat <- pp$rate.cat
    rate.mat <- pp$solution
  }
  # find internal rate.cat dynamics
  RateCats <- paste("R", 1:rate.cat, sep = "")
  edge.rates <- RateMats <- vector("list", rate.cat)
  ind <- matrix(1:dim(rate.mat)[1], dim(rate.mat)[1]/rate.cat, rate.cat)
  # subset the larger matrix into the internal dynamics
  for(i in 1:length(RateCats)){
    RateMats[[i]] <- rate.mat[ind[,i], ind[,i]]
  }
  # if there are hidden rate categories then we need a special matrix to describe them
  if(rate.cat > 1){
    RateCoord <- ind[1,]
    tmpMat1 <- matrix(RateCoord, length(RateCats) , length(RateCats))
    tmpMat2 <- matrix(RateCoord, length(RateCats) , length(RateCats), byrow = TRUE)
    diag(tmpMat1) <- diag(tmpMat2) <- NA
    ClassCoord <- cbind((as.vector(tmpMat1)), (as.vector(tmpMat2)))
    RateClassMat <- matrix(apply(ClassCoord, 1, function(x) rate.mat[x[1], x[2]]), length(RateCats), length(RateCats), dimnames = list(RateCats, RateCats))
    RateMats$RateClassMat <- RateClassMat
  }
  return(RateMats)
}

plotMKmodel <- function(pp, rate.cat = NULL, col.func = plasma, arrow.scale = 2, text.scale = 3){
  
  if(class(pp) == "matrix" & is.null(rate.cat)){
    return(cat("Error: user provided a rate matrix without providing the number of rate categories."))
  }

  # decompose the matrix solution
  if(class(pp) == "corhmm"){
    RateMats <- getRateMats(pp = pp)
    rate.cat <- pp$rate.cat
  }
  if(class(pp) == "matrix"){
    RateMats <- getRateMats(rate.mat = pp, rate.cat = rate.cat)
  }
  ## getting the colours set up and changing mats into vectors
  # the arrows need to be the same colours as the matrix, but i want the arrows to be transparent
  # get the rates from the matrix and ensure it's in the same order as the edges
  edge.rates <- lapply(RateMats, function(x) as.vector(na.exclude(as.vector(t(x)))))
  max.rate <- max(do.call(rbind, edge.rates))
  scaled.edge.rates <- lapply(edge.rates, function(x) x/max.rate)
  cols <- col.func(101, alpha = 0.5)
  edge.cols <- lapply(scaled.edge.rates, function(x) cols[round(x*100)+1])

  # if there is one rate cat there is only one matrix (any more and we need an additional one)
  if(rate.cat == 1){
    xlim <- c(-1,5)
  }else{
    par(mfrow=c(rate.cat+1,1))
    xlim <- c(0,3)
  }
  
  # loop over all the plots
  for(i in 1:length(RateMats)){
    nCol <- dim(RateMats[[i]])[2]
    Coordinates <- getRateMatCoord(RateMats[[i]])
    IndMat <- matrix(1:prod(dim(RateMats[[i]])), dim(RateMats[[i]])[1], dim(RateMats[[i]])[2])
    diag(IndMat) <- NA
    arrow.scalar <- arrow.scale/nCol
    text.scalar <- text.scale/nCol
    
    # setup the matrix with non-diagonal elements
    IndVec <- as.vector(na.omit(as.vector(IndMat)))
    edge.rates <- RateMats[[i]][IndVec]
    scaled.edge.rates <- edge.rates/max.rate
    cols <- col.func(101, alpha = 1)
    mat.cols <- cols[round(scaled.edge.rates*100)+1]
    mat.cols[is.na(mat.cols)] <- "black"
        
    # define the main title of the plots
    main.title <- paste("Rate Category ", i, " (R",i, ")", sep = "")
    if(i == length(RateMats) & rate.cat == 1){
      main.title <- paste("Rate Category ", i, " (R",i, ")", sep = "")
    }
    if(i == length(RateMats) & rate.cat > 1){
      main.title <- "Rate Category Transitions"
    }
    
    # make an igraph object
    tmp <- RateMats[[i]]
    tmp[is.na(tmp)] <- 0
    g1 <- graph_from_adjacency_matrix(adjmatrix = tmp, 
                                      mode = "directed", 
                                      weighted = TRUE, 
                                      diag = FALSE)
    # plot the igraph object
    plot.igraph(g1, xlim = xlim,
                layout = layout_in_circle(g1),
                vertex.label.font = 2,
                vertex.shape = "none",
                vertex.label.cex = 1,
                vertex.label.color = "black",
                vertex.label.family="Helvetica",
                edge.color = edge.cols[[i]],
                edge.width = 2,
                edge.arrow.size = arrow.scalar,
                main = main.title)
    # plot the associated markov matrix
    label.text <- as.character(round(edge.rates, 2))
    label.text[label.text == 0] <- "<0.01"
    label.text[is.na(label.text)] <- "--"
    label.names <- vertex.attributes(g1)$name
    rect(xleft = 1.6, ybottom = -1.3, xright = 4.4, ytop = 1.3, col = rgb(0,0,0,0.2))
    text(Coordinates[[1]], labels = label.text, col = mat.cols, font = 2, cex = text.scalar)
    text(Coordinates[[2]], labels = label.names, col = "black", font = 2, cex = text.scalar)
    text(Coordinates[[3]], labels = c(" ", rep("--", dim(Coordinates[[3]])[1]-1)), font = 2, cex = text.scalar)
  }
}
