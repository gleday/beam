#' Bayesian inference in large Gaussian graphical models
#'
#' @param X n by p data matrix
#' @param type character. Either "marginal", "conditional" or "both". See Details.
#' @param return.only character. Either "cor", "BF", "prob". See details.
#' @param verbose logical. Whether information on progress should be be printed.
#' @param shrinkageMethod character. Either "ml" (marginal likelihood) or "cpo" (conditional predictive ordinate).
#' @param D matrix. Prior marginal correlation matrix. Must be positive definite, well-conditioned and have unit variance.
#'
#' @description 
#' This function carries out covariance and inverse-covariance estimation within the Gaussian conjugate model.
#' The scale matrix parameter of the inverse-Wishart is, by default, set to the identity, whereas the
#' degree of freedom parameter is estimated by maximization of the marginal likelihood (\code{shrinkageMethod = "ml"})
#' or conditional predictive ordinate (\code{shrinkageMethod = "cpo"}).
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
#' Leday, G.G.R. and Richardson, S. (2018). Fast Bayesian inference in large Gaussian graphical models. <arXiv:1803.08155>.
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

beam <- function(X, type = "conditional", return.only = c("cor", "BF", "prob"), verbose=TRUE, shrinkageMethod="ml", D=NULL){

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
					}else{
						if(any(diag(D)!=1)){
							D <- cov2cor(D)
							warning("the matrix D has been scaled to have unit diagonal")
						}
						cholD <- chol(D)
						if(all(diag(cholD)>0)){
							Dinv <- chol2inv(cholD)
							logdetD <- 2*sum(log(diag(cholD)))
						}else{
							stop("D is not positive definite")
						}
					}
				}
			}
		}
	}else{
		logdetD <- NULL
	}
	if(is.character(shrinkageMethod)){
		if(length(shrinkageMethod)!=1){
			stop("shrinkageMethod must be of length 1")
		}else{
			if(!shrinkageMethod%in%c("ml", "cpo")){
				stop("shrinkageMethod must equal 'ml' or 'cpo' ")
			}
		}
	}else{
		stop("shrinkageMethod must be a character")
	}

	# Dimension data
	n <- nrow(X)
	p <- ncol(X)

	# Standardize data
	X <- .standardize(X)

	# Cross-product
	S <- crossprod(X)

	# Compute eigenvalues
	XXT <- NULL
	if(is.null(D)){
		if(n>=p){
			resEigen <- eigen(S, only.values=TRUE)
			eigs <- resEigen$val
		}else{
			XXT <- tcrossprod(X)
			resEigen <- eigen(XXT, only.values=TRUE)
			eigs <- c(resEigen$val, rep(0,p-length(resEigen$val)))
		}
	}else{
		if(n>=p){
			resEigen <- eigen(Dinv%*%S, only.values=TRUE)
			eigs <- resEigen$val
		}else{
			resEigen <- eigen(X%*%Dinv%*%t(X), only.values=TRUE)
			eigs <- c(resEigen$val, rep(0,p-length(resEigen$val)))
		}
	}

	# Store in a vector the upper triangular part of S
	matidxs <- .upperTriIdxs(p)
	Svec <- S[matidxs]

	if(verbose){
	  cat("DONE \n")
	}
	
	###########################################
	#            OPTIMAL SHRINKAGE            #
	###########################################

	if(verbose){
	  cat("Optimal shrinkage... ")
	}
	
	# Range of values
	lowerVal <- .alphaToDelta(alpha=0.001, myn=n, myp=p)
	upperVal <- .alphaToDelta(alpha=0.999, myn=n, myp=p)
	initVal <- .alphaToDelta(alpha=0.5, myn=n, myp=p)

	if(shrinkageMethod=="ml"){

		# Maximize log-marginal likelihood
		resOpt <- optim(par=initVal, fn=.logML, myp=p, myn=n, myeigs=eigs, mylogdetD=logdetD, method='Brent', lower=lowerVal, upper=upperVal, control = list(fnscale=-1))

		# Optimum
		deltaOpt <- resOpt$par
		alphaOpt <- .deltaToAlpha(delta=deltaOpt, myn=n, myp=p)

		# Values of log-ML for a grid of alphas
		myalphas <- seq(from=0.01, to=0.99, length=100)
		mydeltas <- .alphaToDelta(alpha=myalphas, myn=n, myp=p)
		allmls <- sapply(mydeltas, function(x){.logML(x, myp=p, myn=n, myeigs=eigs, mylogdetD=logdetD)})
		gridAlpha <- cbind(mydeltas, myalphas, allmls)
		colnames(gridAlpha) <- c("delta", "alpha", "log-ML")

	}else{

		# Maximize log-cpo score
		resOpt <- optim(par=initVal, fn=.logCPO, myp=p, myn=n, myX=X, myeigs=eigs, myXXT=XXT, myS=S, myD=D, mylogdetD=NULL, method='Brent', lower=lowerVal, upper=upperVal, control = list(fnscale=-1))

		# Optimum
		deltaOpt <- resOpt$par
		alphaOpt <- .deltaToAlpha(delta=deltaOpt, myn=n, myp=p)

		# Values of log-cpo score for a grid of alphas
		myalphas <- seq(from=0.01, to=0.99, length=100)
		mydeltas <- .alphaToDelta(alpha=myalphas, myn=n, myp=p)
		allcpos <- sapply(mydeltas, .logCPO, myp=p, myn=n, myX=X, myeigs=eigs, myXXT=XXT, myS=S, myD=D, mylogdetD=NULL)
		gridAlpha <- cbind(mydeltas, myalphas, allcpos)
		colnames(gridAlpha) <- c("delta", "alpha", "sum log-CPO")

	}

	# Initialize data.frame with results (so that R already knows how much space we are going to use)
	nrows <- p*(p-1)/2
	ncols <- length(return.only)
	if(type == 'both'){
	  results <- as.data.frame(matrix(NA_real_, nrow = nrows, ncol = 2*ncols))
	} else{
	  results <- as.data.frame(matrix(NA_real_, nrow = nrows, ncol = ncols))
	}

	current_free_column <- 1

	if(verbose){
	  cat("DONE \n")
	}
	
	###########################################
	#          MARGINAL DEPENDENCIES          #
	###########################################

	# -> Compute marginal correlations, (scaled) log-Bayes factors and tail probabilities

	# Sample correlations
	rsij <- Svec/n

	# Relieve memory
	rm(Svec)

	if( type=="both" | type=="marginal" ){

	  if(verbose){
	     cat("Marginal dependencies:\n")
	  }

	  if(any( c("cor", "BF") %in% return.only)){

	    if(verbose){
	      cat("--> compute marginal correlation estimates...")
	    }
	    
	    # Prior and posterior marginal correlations
		if(is.null(D)){
			rtij <- (1-alphaOpt)*rsij
		}else{
			rfij <- D[matidxs]
			rtij <- alphaOpt*rfij + (1-alphaOpt)*rsij
		}

	    if("cor" %in% return.only){

	      # Update results table with marginal correlation
	      colnames(results)[current_free_column] <- 'm_cor'  # change col_name
	      results$m_cor <- rtij   # fill column
	      current_free_column = current_free_column + 1  # increment index of free column
	    }

	    if(verbose){
	      cat("DONE \n")
	    }

	    if("BF" %in% return.only){
	      
	      if(verbose){
	        cat("--> Compute scaled Bayes factors...")
	      }
	      
	      # part 1: log Gamma functions
	      part1 <- .lpvarGamma((deltaOpt+n-p+2)/2, p=2) - .lpvarGamma((deltaOpt-p+2)/2, p=2)
	      part1 <- part1 + 2*lgamma((deltaOpt-p+3)/2) - 2*lgamma((deltaOpt+n-p+3)/2)

	      # part 2: log ratio prior/posterior marginal correlations
			if(!is.null(D)){
	      	part2 <- ((deltaOpt-p+2)/2)*log(1-rfij^2) - ((deltaOpt+n-p+2)/2)*log(1-rtij^2)
			}else{
	      	part2 <- - ((deltaOpt+n-p+2)/2)*log(1-rtij^2)
			}

	      # log-BFs
	      logBFs <- part1 + part2

	      # Update results table with marginal log bayes factors
	      colnames(results)[current_free_column] <- 'm_logBF'  # change col_name
	      results$m_logBF <- logBFs   # fill column
	      current_free_column = current_free_column + 1  # increment index of free column

	      # Relieve memory
	      rm(logBFs, part1, part2)

	      if(verbose){
	        cat("DONE \n")
	      }
	      
	    }

	    # Relieve memory
	    rm(rtij)

	  }

	  if("prob" %in% return.only){

	    if(verbose){
	      cat("--> compute tail probabilities...")
	    }
	    
	    # Update table
	    colnames(results)[current_free_column] <- 'm_tail_prob'
	    results$m_tail_prob <- pbeta(rsij^2, 1/2, (n-1)/2, lower.tail=FALSE)
	    current_free_column <- current_free_column + 1

	    # Relieve memory
	    rm(rsij)

	  }

	  if(verbose){
	    cat("DONE \n")
	  }

	}


	###########################################
	#        CONDITIONAL DEPENDENCIES         #
	###########################################

	# -> Compute partial correlations, (scaled) log-Bayes factors and tail probabilities

	if( type=="both" | type=="conditional" ){

	if(verbose){
	    cat("Conditional dependencies:\n")
	    cat("--> compute partial correlation estimates...")
	  }

	  # Compute the inverse of T

	if(is.null(D)){
		if(n>=p){

			Tmat <- (deltaOpt-p-1)*diag(p) + S
			Tinv <- chol2inv(chol(Tmat))

			# Relieve memory
			rm(Tmat)

		}else{

			Tinv <- chol2inv(chol(diag(n)+(1/(deltaOpt-p-1))*XXT))
			Tinv <- crossprod(X, Tinv)/((deltaOpt-p-1)^2)
			Tinv <- Tinv%*%X
			Tinv <- (1/(deltaOpt-p-1))*diag(p) - Tinv

			# Relieve memory
			rm(XXT)

		}
	}else{
		Tmat <- (deltaOpt-p-1)*D + S
		Tinv <- chol2inv(chol(Tmat))

		# Relieve memory
		rm(Tmat)

	}

	# Relieve memory
	rm(S)

	# Prior and posterior partial correlations
	if(!is.null(D)){
		rgij <- - cov2cor(Dinv)
		rgij <- rgij[matidxs]
	}
	rqij <- - cov2cor(Tinv)
	rqij <- rqij[matidxs]

	if("cor" %in% return.only){
		colnames(results)[current_free_column] <- 'p_cor'
		results$p_cor <- rqij
		current_free_column = current_free_column + 1
	}

	if(verbose){
		cat("DONE \n")
	}

	# Compute scaled log-Bayes factors
	if("BF" %in% return.only){

		if(verbose){
			cat("--> compute Bayes factors...")
		}

		# part 1: log Gamma functions
		part1 <- lgamma((deltaOpt+n)/2) - lgamma(deltaOpt/2)
		part1 <- part1 + lgamma((deltaOpt+n-1)/2) - lgamma((deltaOpt-1)/2)
		part1 <- part1 + 2*lgamma((deltaOpt+1)/2) - 2*lgamma((deltaOpt+n+1)/2)

		# part 2: log ratio prior/posterior marginal correlations
		if(!is.null(D)){
			part2 <- 0.5*deltaOpt*log(1-rgij^2) - 0.5*(deltaOpt+n)*log(1-rqij^2)
		}else{
			part2 <- - 0.5*(deltaOpt+n)*log(1-rqij^2)
		}

		# Update table with log-BFs
		colnames(results)[current_free_column] <- 'p_logBF'  # change col_name
		results$p_logBF <- part1 + part2
		current_free_column = current_free_column + 1  # increment index of free column

		# Relieve memory
		rm(part1, part2)

		if(verbose){
			cat("DONE \n")
		}

	}

	  if("prob" %in% return.only){

	    if(verbose){
	      cat("--> compute tail probabilities...")
	    }

	    gii2gjj2 <- (deltaOpt-p-1)^2
	    diagTinv <- diag(Tinv)

	    kii2kjj2 <- tcrossprod(1/diagTinv, 1/diagTinv)
	    kii2kjj2 <- kii2kjj2[matidxs]
	    kii2kjj2 <- kii2kjj2/((1-rqij^2)^2)
	    Gii2Gjj2 <- tcrossprod(diagTinv,diagTinv)
	    Gii2Gjj2 <- Gii2Gjj2[matidxs]
	    tpdiag1 <- diagTinv[matidxs[,1]]
	    tpdiag2 <- diagTinv[matidxs[,2]]
	    Gii2plusGjj2 <- tpdiag1 + tpdiag2

	    # Relieve memory
	    rm(tpdiag1, tpdiag2, diagTinv, matidxs)

	    # Compute rfij
	    # rfij <- (sqrt(kii2kjj2)*tableC[,3])/sqrt(kii2kjj2 + gii2gjj2 - (deltaOpt-p-1)* ( Gii2plusGjj2/(Gii2Gjj2*(1-tableC[,3]^2)) ))
	    num <- sqrt(kii2kjj2)
	    num <- num*rqij
	    denom <- 1-rqij^2
	    denom <- denom*Gii2Gjj2
	    denom <- - (deltaOpt-p-1)*(Gii2plusGjj2/denom)
	    denom <- denom + kii2kjj2
	    denom <- denom + gii2gjj2
	    rfij <- num/sqrt(denom)

	    # Relieve memory
	    rm(Gii2plusGjj2, Gii2Gjj2, kii2kjj2, gii2gjj2)

	    # Update table
	    rfij2 <- rfij^2

	    colnames(results)[current_free_column] <- 'p_tail_prob'  # change col_name
	    #print(mean((rfij2)>qbeta(0.1, 1/2, (n-1)/2, lower.tail=FALSE)))
	    results$p_tail_prob <- pbeta(rfij2, 1/2, (n-1)/2, lower.tail=FALSE)

	    # Relieve memory
	    rm(rfij, rfij2)

	    if(verbose){
	      cat("DONE \n")
	    }
	    
	  }

	  # Relieve memory
	  rm(Tinv, rqij)

	}

	time1 <- proc.time() - time0
	#########################
	#        OUTPUT         #
	#########################

	# List
	out <- new("beam",
	           "table" = results,
	           "deltaOpt" = deltaOpt,
	           "alphaOpt"= alphaOpt,
	           "dimX"= c(n,p),
	           "type"= type,
	           "varlabs" = varlabs,
	           "gridAlpha" = gridAlpha,
	           "valOpt" = resOpt$value,
	           "return.only" = return.only,
	           "time" = time1[3])

	return(out)

}

