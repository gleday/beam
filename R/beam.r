#' Bayesian inference in large Gaussian graphical models
#'
#' @param X n by p data matrix
#' @param type character. Either "marginal", "conditional" or "both". See Details.
#' @param return.only character. Either "cor", "BF", "prob". See details.
#' @param verbose logical. Whether information on progress should be be printed.
#' @param D matrix. Prior marginal correlation matrix. Must be positive definite, well-conditioned and have unit variance.
#'
#' @description 
#' This function carries out covariance and inverse-covariance estimation within the Gaussian conjugate model.
#' The scale matrix parameter of the inverse-Wishart is, by default, set to the identity, whereas the
#' degree of freedom parameter is estimated by maximization of the marginal likelihood.
#' The function also computes Bayes factors and tail probabilities (p-values) to recover the marginal and/or
#' conditional independence structure between variables.
#'
#' @details
#' The arguments \code{type} and \code{return.only} have essentially been introduced for computational and memory savings.
#' Using argument \code{type} the user may indicate whether the marginal dependencies ("marginal"), conditional dependencies
#' ("conditional") or both ("both") are to be inferred. On the other hand, the argument \code{return.only} is used to indicate
#' whether the correlations ("cor"), Bayes factors ("BF") or tail probabilities ("prob") should be returned.
#' Default is to return all three quantities both for marginal and conditional dependencies.
#'
#' @return An object of class \code{\link{beam-class}}
#'
#' @author Gwenael G.R. Leday and Ilaria Speranza
#'
#' @references
#' Leday, G.G.R. and Richardson, S. (2019). Fast Bayesian inference in large Gaussian graphical models. \emph{Biometrics}.
#'
#' @examples
#' # Load data
#' data(TCPAprad)
#' 
#' # beam
#' fit <- beam(X = TCPAprad, type="both") 
#' 
#' # Print summary
#' summary(fit)
#' 
#' # Extract matrix of marginal correlations
#' mcor(fit)[1:5, 1:5]
#' 
#' # Extract matrix of partial correlations
#' pcor(fit)[1:5, 1:5]
#' 
#' # Plot log-marginal likelihood of the Gaussian conjugate model
#' plotML(fit)
#' 
#' # Plot heatmap of marginal (upper triangle) and/or
#' # partial (lower triangle) correlation estimates
#' plotCor(fit)
#' 
#' @export

beam <- function(X, type = "conditional", return.only = c("cor", "BF", "prob"), verbose=TRUE, D=NULL){

	time0 <- proc.time()

	###########################################
	#              PREPROCESSING              #
	###########################################

	if(verbose){
	   cat("Preprocessing... ")
	}

	# Check input argument X
	#stopifnot(is.matrix(X),  !any(is.na(X)), !any(is.infinite(X)))
	if(!is.matrix(X)){
		stop("X is not a matrix")
	}else{
		if(any(is.na(X))){
			stop("X contains missing values")
		}else{
		   if(!is.null(colnames(X))){
		     varlabs <- colnames(X)
		   }else{
				 varlabs <- character()
		   }
		}
	}
	# Check input argument type
	if(is.character(type)){
	  if(length(type)==1){
	    if(!type%in%c("both", "marginal", "conditional")){
	      stop("type is not recognized")
	    }
	  }else{
	    stop("type must be a character of length equal to 1")
	  }
	}else{
	  stop("type must be a character")
	}
	# Check input argument return.only
	if(is.character(return.only)){
	  if(length(return.only)%in%c(1:3)){
      ind.return.only <- return.only%in%c("cor", "BF", "prob")
	    if(any(!ind.return.only)){
	      stop("return.only contains characters that are not recognized")
	    }
	  }else{
	    stop("return.only must contain at least 1 element and 3 elements at most")
	  }
	}else{
	  stop("return.only must be a character")
	}
	# Check input argument D
	if(!is.null(D)){
		if(!is.matrix(D)){
			stop("D must be a matrix")
		}else{
			if(nrow(D)!=ncol(D)){
				stop("D must be a square matrix")
			}else{
				if(!isSymmetric(D)){
					stop("D must be a symmetric matrix")
				}else{
					if(any(is.na(D))){
						stop("D must not contain missing values")
					}
				}
			}
		}
	}else{
		D <- matrix(0, ncol(X), ncol(X))
	}

	#########################
	#         BEAM          #
	#########################
	
	ind <- c("cor", "BF", "prob") %in% return.only
	res <- .beam(X=X, type=type, ronly=ind*1, D = D, verbose = verbose)
	
	labs <- NULL
	if(type=="marginal" || type=="both"){
	  labs <- c(labs, paste0("m_", c("cor", "logBF", "tail_prob")[ind]))
	}
	if(type=="conditional" || type=="both"){
	  labs <- c(labs, paste0("p_", c("cor", "logBF", "tail_prob")[ind]))
	}
	colnames(res$table) <- labs
	
	time1 <- proc.time() - time0

	#########################
	#        OUTPUT         #
	#########################

	# List
	out <- new("beam",
	           "table" = res$table,
	           "deltaOpt" = res$deltaOpt,
	           "alphaOpt"= res$alphaOpt,
	           "dimX"= dim(X),
	           "type"= type,
	           "varlabs" = varlabs,
	           "gridAlpha" = res$gridAlpha,
	           "valOpt" = res$valOpt,
	           "return.only" = return.only,
	           "time" = time1[3])

	return(out)

}

