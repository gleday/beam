#' Bayesian inference for large-scale Gaussian graphical models
#'
#' @param X n by p data matrix
#' @param type character. Either "marginal", "conditional" or "both". See Details.
#' @param return.only character. Either "cor", "BF", "prob". See details.
#' @param verbose logical. Whether information on progress should be be printed.
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
#' @references Leday, G.G.R. and Richardson, S. (2018). Fast Bayesian inference in large Gaussian graphical models. Submitted.
#'
#' @export

beam <- function(X, type = "conditional", return.only = c("cor", "BF", "prob"), verbose=TRUE){

	time0 <- proc.time()

	###########################################
	#              PREPROCESSING              #
	###########################################

	if(verbose){
	   cat("Preprocessing... ")
	}

	# Check arguments
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

	# Dimension data
	n <- nrow(X)
	p <- ncol(X)

	# Standardize data
	X <- scale(X, center = TRUE, scale = FALSE)
	X <- scale(X, center = TRUE, scale = sqrt(colSums(X^2)/n))

	# Cross-product
	S <- crossprod(X)

	# Compute vector of eigenvalues
	if(n>=p){
		resEigen <- eigen(S, only.values=FALSE)
		eigs <- resEigen$val
	}else{
		resEigen <- eigen(tcrossprod(X), only.values=FALSE)
		eigs <- c(resEigen$val, rep(0,p-length(resEigen$val)))
	}

	# Store in a vector only the upper triangular part of S
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

	# Maximize log-marginal likelihood
	resOpt <- optim(par=initVal, fn=.logML, myp=p, myn=n, myeigs=eigs, method='Brent', lower=lowerVal, upper=upperVal, control = list(fnscale=-1))

	# Optimum
	deltaOpt <- resOpt$par
	alphaOpt <- .deltaToAlpha(delta=deltaOpt, myn=n, myp=p)

	# Keep values of log-ML for a grid of alphas
	myalphas <- seq(from=0.01, to=0.99, length=100)
	mydeltas <- .alphaToDelta(alpha=myalphas, myn=n, myp=p)
	allmls <- sapply(mydeltas, function(x){.logML(x, myp=p, myn=n, myeigs=eigs)})
	gridAlpha <- cbind(mydeltas, myalphas, allmls)
	colnames(gridAlpha) <- c("delta", "alpha", "logML")

	# Initialize data.frame with results (so that R already knows how much space we are going to use)
	nrows = p*(p-1)/2
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

	# All the off diag elements of M are null (M <- (deltaOpt-p-1)*diag(p))
	Hvec <- Svec

	if( type=="both" | type=="marginal" ){

	  if(verbose){
	     cat("Marginal dependencies:\n")
	  }

	  if(any( c("cor", "BF") %in% return.only)){

	    if(verbose){
	      cat("--> compute marginal correlation estimates...")
	    }
	    
	    # Marginal correlation estimates
	    rhij <- Hvec/(deltaOpt+n-p-1)

	    if("cor" %in% return.only){

	      # Update results table with marginal correlation
	      colnames(results)[current_free_column] <- 'm_cor'  # change col_name
	      results$m_cor <- rhij   # fill column
	      current_free_column = current_free_column + 1  # increment index of free column
	    }

	    if(verbose){
	      cat("DONE \n")
	    }

	    if("BF" %in% return.only){
	      
	      if(verbose){
	        cat("--> Compute Bayes factors...")
	      }
	      
	      # part 1: log Gamma functions
	      part1 <- .lpvarGamma((deltaOpt+n-p+2)/2, p=2) - .lpvarGamma((deltaOpt-p+2)/2, p=2)
	      part1 <- part1 + 2*lgamma((deltaOpt-p+3)/2) - 2*lgamma((deltaOpt+n-p+3)/2)

	      # part 2: log ratio prior/posterior marginal correlations
	      part2 <- - ((deltaOpt+n-p+2)/2)*log(1-rhij^2)

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
	    rm(rhij)

	  }

	  if("prob" %in% return.only){

	    if(verbose){
	      cat("--> compute tail probabilities...")
	    }
	    
	    # Sample correlations
	    rsij <- Svec/n

	    # Update table
	    colnames(results)[current_free_column] <- 'm_tail_prob'  # change col_name
      #print(mean((rsij^2)>qbeta(0.1, 1/2, (n-1)/2, lower.tail=FALSE)))
	    results$m_tail_prob <- pbeta(rsij^2, 1/2, (n-1)/2, lower.tail=FALSE)   # fill the free column
	    current_free_column = current_free_column + 1  # increment index of free column

	    # Relieve memory
	    rm(rsij)

	  }

	  if(verbose){
	    cat("DONE \n")
	  }

	}

	# Relieve memory
	rm(Svec, Hvec)


	###########################################
	#        CONDITIONAL DEPENDENCIES         #
	###########################################

	# -> Compute partial correlations, (scaled) log-Bayes factors and tail probabilities

	if( type=="both" | type=="conditional" ){

	  if(verbose){
	    cat("Conditional dependencies:\n")
	    cat("--> compute partial correlation estimates...")
	  }

	  # Compute scaled posterior expectation of the inverse covariance matrix (partial correlation estimates)
	  if(n>=p){
	    H <- (deltaOpt-p-1)*diag(p) + S
	    Hinv <- solve(H)
	    rm(H) # relieve memory
	  }else{
	    prodtp <- crossprod(X, resEigen$vectors)
	    Hinv <- tcrossprod(prodtp,diag(1/(1+(1/(deltaOpt-p-1))*resEigen$val)))
	    Hinv <- tcrossprod( Hinv, prodtp)/((deltaOpt-p-1)^2)
	    Hinv <- (1/(deltaOpt-p-1))*diag(p) - Hinv
	    #Hinv <- (1/(deltaOpt-p-1))*diag(p) - (1/(deltaOpt-p-1)^2)*t(X)%*%solve(diag(n)+(1/(deltaOpt-p-1))*tcrossprod(X))%*%X

	    # Relieve memory
	    rm(prodtp)

	  }

	  # Relieve memory
	  rm(S)

	  # Partial correlation estimates
	  rkij <- - cov2cor(Hinv)
	  rkij <- rkij[matidxs]  # keep only the upper triangular part, discard the rest

	  if("cor" %in% return.only){
	    colnames(results)[current_free_column] <- 'p_cor'  # change col_name
	    results$p_cor <- rkij  # fill the free column
	    current_free_column = current_free_column + 1  # increment index of free column
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
	    part2 <- - ((deltaOpt+n)/2)*log(1-rkij^2)

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
	    diagHinv <- diag(Hinv)

	    kii2kjj2 <- tcrossprod(1/diagHinv, 1/diagHinv)
	    kii2kjj2 <- kii2kjj2[matidxs]
	    kii2kjj2 <- kii2kjj2/((1-rkij^2)^2)
	    Gii2Gjj2 <- tcrossprod(diagHinv,diagHinv)
	    Gii2Gjj2 <- Gii2Gjj2[matidxs]
	    tpdiag1 <- diagHinv[matidxs[,1]]
	    tpdiag2 <- diagHinv[matidxs[,2]]
	    Gii2plusGjj2 <- tpdiag1 + tpdiag2

	    # Relieve memory
	    rm(tpdiag1, tpdiag2, diagHinv, matidxs)

	    # Compute rfij
	    # rfij <- (sqrt(kii2kjj2)*tableC[,3])/sqrt(kii2kjj2 + gii2gjj2 - (deltaOpt-p-1)* ( Gii2plusGjj2/(Gii2Gjj2*(1-tableC[,3]^2)) ))
	    num <- sqrt(kii2kjj2)
	    num <- num*rkij
	    denom <- 1-rkij^2
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
	  rm(Hinv, rkij)

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

