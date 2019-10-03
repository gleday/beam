#' Edge selection with multiple testing and error control
#'
#' @param object An object of class \code{\link{beam-class}}.
#' @param thres numeric. Significance threshold to be applied on adjusted tail probabilities.
#' @param method character. Method to use for multiple comparison adjustment of tail probabilities.
#' @param return.only character. Quantities to be returned.
#'
#' @description 
#' Infer graphical structures by multiple testing
#'
#' @details
#' The argument \code{method} allows to adjust the tail probabilities obtained from the null distributions of
#' the Bayes factors for multiple comparisons. Possible choices are: "holm", "bonferroni", "BH", "BY" and "HC".
#' Apart from "HC", these are passed onto the R function \code{p.adjust}
#' from package \pkg{stats} and we refer the user to its documentation for details. The method "HC" provides an
#' optimal decision threshold based on the Higher Criticism score which is computed using the R function \code{hc.thresh}
#' from package \pkg{fdrtool}. Again, we refer to the associated documentation for details.
#'
#' The argument \code{return.only} allows to decide which quantities have to be in the output:
#' it could be any subvector of c('cor', 'BF', 'prob', 'adj')
#' (provided that the requested quantities have been computed in the beam object, except for adjusted probabilities).
#' It can also be set to NULL: in this case, only the selected edges will be returned without any additional information.
#' The default value for this argument are the columns present in the beam object plus the adjusted probabilities.
#'
#' @return An object of class \code{\link{beam.select-class}}
#'
#' @author Gwenael G.R. Leday and Ilaria Speranza
#'
#' @references
#' Drton, M., & Perlman, M. D. (2007). Multiple testing and error control in Gaussian graphical model selection. Statistical Science, 430-449.\cr
#' Goeman, J. J., & Solari, A. (2014). Multiple hypothesis testing in genomics. Statistics in medicine, 33(11), 1946-1978.\cr
#' Donoho, D., & Jin, J. (2015). Higher criticism for large-scale inference, especially for rare and weak effects. Statistical Science, 30(1), 1-25.\cr
#' Klaus, B., & Strimmer, K. (2012). Signal identification for rare and weak features: higher criticism or false discovery rates?. Biostatistics, 14(1), 129-143.
#'
#' @export
beam.select <- function(object, thres = 0.1, method = "BH", return.only = c(object@return.only,'adj')){

  ###########################################
  #              PREPROCESSING              #
  ###########################################

  # Check arguments
  if(!is.numeric(thres)){
    stop("thres must be a numeric")
  }else{
    if((thres<=0) | (thres>=1)){
      stop("thres must be greater than 0 and lower than 1")
    }
  }
  if(!is.null(method)){
    if(is.character(method)){
      if(!method%in%c("holm", "bonferroni", "BH", "BY", "HC")){
        stop("method is not recognized")
      }
    }else{
      stop("method must be a character")
    }
  }

  if(is.character(return.only)){
    if(length(return.only)%in%c(1:4)){
      ind.return.only <- return.only%in%c("cor", "BF", "prob","adj")
      if(any(!ind.return.only)){
        stop("return.only contains characters that are not recognized")
      }
    }else{
      stop("return.only must contain at least 1 element and 4 elements at most")
    }
  }else{
    if(!is.null(return.only)){
      stop("return.only must be a character")
    }
  }

  if(any(c("cor", "BF", "prob") %in% setdiff(return.only, object@return.only))){
    stop("Some of the requested information is not available in the input beam object")
  }

  # get the data.frame with marginal and/or conditional estimations
  df <- object@table
  p <- object@dimX[2]  # get number of covariates
  df.cols <- colnames(df)

  ###########################################
  #        MARGINAL INDEPENDENCE GRAPH      #
  ###########################################

  if(object@type == 'conditional'){

    tableM <- as.data.frame(matrix(nrow=0, ncol=0))

  }else{
    marg.cols <- c()
    if('cor' %in% return.only){
      marg.cols <- c(marg.cols, 'm_cor')
    }
    if('BF' %in% return.only){
      marg.cols <- c(marg.cols, 'm_logBF')
    }

    if('prob' %in% return.only){
      marg.cols <- c(marg.cols, 'm_tail_prob')
    }

    if('adj' %in% return.only){ # add a column for adjusted probabilities
      tableM <- as.data.frame(df[, marg.cols])
      tableM <- cbind(tableM, vector(mode = 'numeric', length = p*(p-1)/2))
      names(tableM)[length(names(tableM))] <- paste0("m_tail_prob_", ifelse(method == 'blfdr',
                                                                            substr(method, 1, 5),
                                                                            substr(method, 1, 4))) # modify last column name
    }else{
      tableM <- as.data.frame(df[, marg.cols])
    }
print(class(tableM))
    if(!('m_tail_prob' %in% df.cols)){
      stop('Method ', method, ' requires tail probabilities, which are not currently included in the beam object')
    }

    if(method=="HC"){

      HCthres <- hc.thresh(df[,"m_tail_prob"], alpha0=1, plot=FALSE)
      if('adj' %in% return.only){
        tableM <- tableM[df[,"m_tail_prob"] < HCthres, -ncol(tableM)]
      }else{
        tableM <- tableM[df[,"m_tail_prob"] < HCthres, , drop = FALSE ]
      }

    }else{

      tailAdj <- p.adjust(df[,"m_tail_prob"], method=method)
      if('adj' %in% return.only){
        tableM[, ncol(tableM)] <- tailAdj
      }
      tableM <- tableM[tailAdj<thres,  , drop = FALSE]

    }

  }

  ###########################################
  #      CONDITIONAL INDEPENDENCE GRAPH     #
  ###########################################

  if(object@type == 'marginal'){

    tableC <- as.data.frame(matrix(nrow=0, ncol=0))

  }else{

    cond.cols <- c()
    if('cor' %in% return.only){
      cond.cols <- c(cond.cols, 'p_cor')
    }
    if('BF' %in% return.only){
      cond.cols <- c(cond.cols, 'p_logBF')
    }

    if('prob' %in% return.only){
      cond.cols <- c(cond.cols, 'p_tail_prob')
    }

    if('adj' %in% return.only){ # add a column for adjusted probabilities
      tableC <- as.data.frame(df[, cond.cols])
      tableC <- cbind(tableC, vector(mode = 'numeric', length = p*(p-1)/2))
      names(tableC)[length(names(tableC))] <- paste0("p_tail_prob_", ifelse(method == 'blfdr',
                                                                            substr(method, 1, 5),
                                                                            substr(method, 1, 4))) # modify last column name
    }else{
      tableC <- as.data.frame(df[,cond.cols])
    }

    if(!('p_tail_prob' %in% df.cols)){
      stop('Method ', method, 'requires tail probabilities, which are not currently included in the beam object')
    }

    if(method=="HC"){

      HCthres <- hc.thresh(df[,"p_tail_prob"], alpha0=1, plot=FALSE)
      if('adj' %in% return.only){  # put values in column if requested in output
        tableC <- tableC[df[,"p_tail_prob"] < HCthres, -ncol(tableC)]
      }else{
        tableC <- tableC[df[,"p_tail_prob"] < HCthres,  , drop = FALSE]
      }

    }else{

      tailAdj <- p.adjust(df[,"p_tail_prob"], method=method)
      if('adj' %in% return.only){  # put values in column if requested in output
        tableC[, ncol(tableC)] <- tailAdj
      }
      tableC <- tableC[tailAdj<thres,  , drop = FALSE]

    }

  }
  print(class(tableM))
  #########################
  #        OUTPUT         #
  #########################

  # List
  out <- new("beam.select",
             "marginal" = tableM,
             "conditional"= tableC,
             "dimX"= object@dimX,
             "type"= object@type,
             "varlabs" = object@varlabs,
             "alphaOpt"= object@alphaOpt,
             "gridAlpha" = object@gridAlpha,
             "valOpt" = object@valOpt,
             "method" = method,
             "thres" = thres
             )

  return(out)

}
