#######################
## Internal R functions
#######################

.upperTriIdxs <- function(p){
  z <- sequence(p)
  cbind(
    row = unlist(lapply(1:(p-1), function(x) 1:x), use.names = FALSE),
    col = rep(z[-1], times = tail(z, -1)-1))
}

.Idx2RowCol <- function(idx){
  col = ceiling((1 + sqrt(1+8*idx))/2)
  row = idx - (col^2 - 3*col +4)/2 + 1
  return(data.frame(row=row, col=col))
}

# .logCPO <- function(mydelta, myp, myn, myX, myeigs, myXXT, myS, myD, mylogdetD){
#   part1 <- -0.5*myp*log(pi) + .lpvarGamma((mydelta+myn)*0.5, p=myp) - .lpvarGamma((mydelta+myn-1)*0.5, p=myp) 
#   part2 <- -0.5*sum(log((mydelta-myp-1)+myeigs))
#   mycpo <- part1 + part2
#   if(!is.null(mylogdetD)){
#     mycpo <- mycpo - 0.5*mylogdetD
#     TinvDelta <- chol2inv(chol((mydelta-myp-1)*myD + myS))
#   }else{
#     if(is.null(myXXT)){
#       TinvDelta <- chol2inv(chol((mydelta-myp-1)*diag(myp) + myS))
#     }else{
#       TinvDelta <- chol2inv(chol(diag(myn)+(1/(mydelta-myp-1))*myXXT))
#       TinvDelta <- crossprod(myX, TinvDelta)/((mydelta-myp-1)^2)
#       TinvDelta <- TinvDelta%*%myX
#       TinvDelta <- (1/(mydelta-myp-1))*diag(myp) - TinvDelta
#     }
#   }
#   mycpo <- mycpo + 0.5*(mydelta+myn-1)*log(1-diag((myX%*%TinvDelta)%*%t(myX)))
#   return(sum(mycpo))
# }

