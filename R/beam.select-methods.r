###################################
## Methods for "beam.select-class"
###################################

#' @rdname beam.select-class
#' @aliases print
#' @param x An object of class \code{beam.select-class}
#' @param ... further arguments passed to or from other methods.
setMethod(
  f = "print",
  signature = "beam.select",
  definition = function(x,...){
    cat("Class \"beam.select\"\n")
  }
)

#' @rdname beam.select-class
#' @aliases show
#' @param object An object of class \code{beam.select-class}
setMethod(
  f = "show",
  signature = "beam.select",
  definition = function(object) print(object)
)

#' @rdname beam.select-class
#' @aliases summary
setMethod(
  f = "summary",
  signature = "beam.select",
  definition = function(object, ...){
    cat("\n")
    cat("Object of class \"beam.select\"\n\n")
    cat("+ INPUT DATA \n")
    cat("number of observations:", object@dimX[1], "\n")
    cat("number of variables:", object@dimX[2], "\n\n")
    
    p <- object@dimX[2]
    
    cat("+ SELECTION METHOD \n")
    cat("Adjustment method: ")
    if(object@method=="BH"){
      cat("Benjamini-Hochberg\n")
    }
    if(object@method=="holm"){
      cat("Holm\n")
    }
    if(object@method=="bonferroni"){
      cat("Bonferroni\n")
    }
    if(object@method=="BY"){
      cat("Benjamini-Yekutieli\n")
    }
    if(object@method=="HC"){
      cat("Higher Criticism\n")
    }
    cat("Significance threshold:", object@thres,"\n\n")
    
    if(!(object@type == 'conditional')){
      cat("+ MARGINAL INDEPENDENCE GRAPH\n")
      marg <- object@marginal
      
      tot.edges = p*(p-1)/2
      
      edges <- .Idx2RowCol(as.numeric(rownames(marg))) # dataframe with all edges
      vertices <- unique(union(edges[,1], edges[,2]))
      
      sel.edges <- nrow(edges)
      sel.vertices <- length(vertices)
      tot.vertices <- p
      
      edge.perc <- sel.edges/tot.edges*100
      vert.perc <- sel.vertices/tot.vertices*100
      
      cat(format(sel.edges, big.mark=","), " selected edges out of ", format(tot.edges, big.mark=","), " (", round(edge.perc, digits=2),"%)\n", sep="")
      cat(format(sel.vertices, big.mark=","), " nodes with nonzero degree out of ", format(tot.vertices, big.mark=","), " (", round(vert.perc, digits=2), "%)\n\n", sep="")
      
    }
    if(!(object@type == 'marginal')){
      cat("+ CONDITIONAL INDEPENDENCE GRAPH\n")
      cond <- object@conditional
      
      tot.edges = p*(p-1)/2
      
      edges <- .Idx2RowCol(as.numeric(rownames(cond))) # dataframe with all edges
      vertices <- unique(union(edges[,1], edges[,2]))
      
      sel.edges <- nrow(edges)
      sel.vertices <- length(vertices)
      tot.vertices <- p
      
      edge.perc <- sel.edges/tot.edges*100
      vert.perc <- sel.vertices/tot.vertices*100
      
      cat(format(sel.edges, big.mark=","), " selected edges out of ", format(tot.edges, big.mark=","), " (", round(edge.perc, digits=2),"%)\n", sep="")
      cat(format(sel.vertices, big.mark=","), " nodes with nonzero degree out of ", format(tot.vertices, big.mark=","), " (", round(vert.perc, digits=2), "%)\n", sep="")
      
    }
    cat("\n")
  }
)

#' @rdname beam.select-class
#' @aliases marg
#' @param object An object of class \code{beam.select-class}
setMethod(
  f = "marg",
  signature = "beam.select",
  definition = function(object){
    if(object@type == 'conditional'){
      warning('No information about marginal dependence structure available. Check "type" argument in beam call')
      return(data.frame())
    }else{
      sel.idxs <- as.numeric(rownames(object@marginal))
      rowcol <- .Idx2RowCol(sel.idxs)
      return(cbind(rowcol, object@marginal))
    }

  }
)

#' @rdname beam.select-class
#' @aliases cond
#' @param object An object of class \code{beam.select-class}
setMethod(
  f = "cond",
  signature = "beam.select",
  definition = function(object){
    if(object@type == 'marginal'){
      warning('No information about condtional dependence structure available. Check "type" argument in beam call')
      return(data.frame())
    }else{
      sel.idxs <- as.numeric(rownames(object@conditional))
      rowcol <- .Idx2RowCol(sel.idxs)
      return(cbind(rowcol, object@conditional))
    }
  }
)


#' @rdname beam.select-class
#' @aliases mcor
#' @param object An object of class \code{beam.select-class}
setMethod(
  f = "mcor",
  signature = "beam.select",
  definition = function(object){
    df <- object@marginal
    p <- object@dimX[2]
    idxs <- .Idx2RowCol(as.numeric(rownames(df)))

    if('m_cor' %in% colnames(df)){
      matcor <- as.matrix(Matrix::sparseMatrix(i=c(idxs$row,idxs$col), j=c(idxs$col,idxs$row), x=rep(df$m_cor,2), dims=c(p,p)))
      diag(matcor) <- 1
      return(matcor)
    }else{
      stop('Method not available: either "cor" is not included in return.only or type = "marginal" in calling function beam')
    }
  }
)

#' @rdname beam-class
#' @aliases pcor
#' @param object An object of class \code{beam-class}
setMethod(
  f = "pcor",
  signature = "beam.select",
  definition = function(object){
    df <- object@conditional
    p <- object@dimX[2]
    idxs <- .Idx2RowCol(as.numeric(rownames(object@conditional)))

    if('p_cor' %in% colnames(df)){
      matcor <- as.matrix(Matrix::sparseMatrix(i=c(idxs$row,idxs$col), j=c(idxs$col,idxs$row), x=rep(df$p_cor,2), dims=c(p,p)))
      diag(matcor) <- 1
      return(matcor)
    }else{
      stop('Method not available: either "cor" is not included in return.only or type = "marginal" in calling function beam')
    }

  }
)

#' @rdname beam.select-class
#' @aliases plotML
setMethod(
  f = "plotML",
  signature = "beam.select",
  definition = function(object, ...){
    plot(object@gridAlpha[,2], object@gridAlpha[,3], type="l", xlab=expression(alpha), ylab="log-marginal likelihood")
    abline(v=object@alphaOpt, col="black", lty=2)
    abline(h=object@valOpt, col="black", lty=2)
  }
)

#' @rdname beam.select-class
#' @aliases plotAdj
#' @param type character. Type of dependence to be displayed (marginal, conditional or both)
#' @param order character. "original" or "clust"
setMethod(
  f = "plotAdj",
  signature = "beam.select",
  definition = function(object, type=object@type, order = "original"){
    if(!type%in%c("marginal", "conditional", "both")){
      stop("type must be equal to 'marginal', 'conditional' or 'both' ")
    }else{
      # Set default
      if(type=="both"){
        type <- "conditional"
      }
      if(!order%in%c('original', 'clust')){
        stop("order must be equal to 'original' or 'clust'")
      }else{
        p <- object@dimX[2]
  
        # Extract correlations
        if(type=="marginal"){
          if(object@type %in% c('both', 'marginal')){
            marg <- marg(object)
            adjMat <- as.matrix(Matrix::sparseMatrix(i=c(marg$row,marg$col), j=c(marg$col,marg$row), dims=c(p,p)))
            
            if(order == 'clust'){
              bgr <- bgraph(object)
              clust <- igraph::cluster_louvain(bgr)
              memb <- igraph::membership(clust)
              ordered.idxs <- order(as.vector(memb))
              adjMat <- adjMat[ordered.idxs,ordered.idxs]
  
            }
          }else{
            stop('Method not available: no information about marginal dependence structure')
          }
        }
  
        if(type=="conditional"){
          if(object@type %in% c('both', 'conditional')){
            cond <- cond(object)
            adjMat <- as.matrix(Matrix::sparseMatrix(i=c(cond$row,cond$col), j=c(cond$col,cond$row), dims=c(p,p)))
            if(order == 'clust'){
              ugr <- ugraph(object)
              clust <- igraph::cluster_louvain(ugr)
              memb <- igraph::membership(clust)
              ordered.idxs <- order(as.vector(memb))
              adjMat <- adjMat[ordered.idxs,ordered.idxs]
  
            }
          }else{
            stop('Method not available: no information about conditional dependence structure')
          }
        }
  
        if(order == 'clust'){
          cat('Number of clusters identified with algorithm "cluster_louvain":', max(memb), '\n')
          cat('Number of nodes in each cluster:', igraph::sizes(clust), '\n')
          cat('Modularity =', igraph::modularity(clust), '\n')
        }
  
        # Function needed
        mirror <- function(mymat){
          xx <- as.data.frame(mymat);
          xx <- rev(xx);
          xx <- as.matrix(xx);
          xx;
        }
  
        # Plot
        par(mar=c(2, 2, 2, 2) + 0.1)
        image(1:nrow(adjMat), 1:ncol(adjMat), mirror(adjMat), zlim=c(0,1), col=c("white", "black"), breaks=c(0, 0.5, 1), xlab="", ylab="", main="", xaxt="n", yaxt="n")
        box()
      }
    }
  }
)

#' @rdname beam.select-class
#' @aliases bgraph
setMethod(
  f = "bgraph",
  signature = "beam.select",
  definition = function(object){
    if(object@type %in% c('marginal', 'both')){
      marg <- object@marginal
      edges <- .Idx2RowCol(as.numeric(rownames(marg))) # dataframe with all edges
      colnames(edges) <- c('node1','node2')
      if(length(object@varlabs)>0){
        vertices <- data.frame(id = 1:object@dimX[2], label = object@varlabs)
      }else{
        vertices <- data.frame(id = 1:object@dimX[2], label = as.character(1:object@dimX[2]))
      }
      myigraph <- igraph::graph_from_data_frame(d=edges, vertices=vertices, directed=FALSE)
      if(length(object@varlabs)>0){
        myigraph <- igraph::set.vertex.attribute(myigraph, "name", value=object@varlabs)
      }
      return(myigraph)
    }else{
      stop('Method not available: no information about marginal dependence structure')
    }
  }
)

#' @rdname beam.select-class
#' @aliases ugraph
setMethod(
  f = "ugraph",
  signature = "beam.select",
  definition = function(object){
    if(object@type %in% c('conditional', 'both')){
      cond <- object@conditional
      edges <- .Idx2RowCol(as.numeric(rownames(cond)))
      colnames(edges) <- c('node1','node2')
      if(length(object@varlabs)>0){
        vertices <- data.frame(id = 1:object@dimX[2], label = object@varlabs)
      }else{
        vertices <- data.frame(id = 1:object@dimX[2], label = as.character(1:object@dimX[2]))
      }
      myigraph <- igraph::graph_from_data_frame(d=edges, vertices=vertices, directed=FALSE)
      if(length(object@varlabs)>0){
        myigraph <- igraph::set.vertex.attribute(myigraph, "name", value=object@varlabs)
      }
      return(myigraph)
    }else{
      stop('Method not available: no information about conditional dependence structure')
    }
  }
)
