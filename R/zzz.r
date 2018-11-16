#######################
## Internal functions
#######################

.lpvarGamma <- function(z, p){# functions written by Harry Gray
  position = 1:p
  if (length(z) == 1){
    ans <- ((p * (p-1))/4) * log(pi) + (sum(lgamma(z + (1-position)/2)))
  }else{
    ans <- rep(0, length(z))
    for (i in 1:length(z)){
      ans[i] <- ((p * (p-1))/4) * log(pi) + (sum(lgamma(z[i] + (1-position)/2)))
    }
  }
  return(ans)
}

.alphaToDelta <- function(alpha, myn, myp){
  (alpha*myn+(1-alpha)*myp+(1-alpha))/(1-alpha)
}

.deltaToAlpha <- function(delta, myn, myp){
  (delta-myp-1)/(myn+delta-myp-1)
}

.logML <- function(mydelta, myp, myn, myeigs){
  part1 <- -0.5*myn*myp*log(pi) + .lpvarGamma((mydelta+myn)*0.5, p=myp) - .lpvarGamma(mydelta*0.5, p=myp)
  part2 <- 0.5*mydelta*myp*log(mydelta-myp-1) - 0.5*(mydelta+myn)*sum(log((mydelta-myp-1)+myeigs))
  part1+part2
}

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
