#' Probability that n initial cases lead to an outbreak that lasts at least g generations
#' of transmission AND has exactly j total cases after generation g
#'
#' @param g Number of generations of transmission
#' @param n Number of initial cases
#' @param j Total size of outbreak after generation g
#' @param R Reproduction number: mean of negative binomial offspring distribution
#' @param k Dispersion parameter of negative binomial offspring distribution
#' @author Damon Toth
#' @returns The probability of the given outbreak size at the given generation
#' @export
pSizeAtGen <- function(g,n,j,R,k){
  
  if(g==1){
    out <- pNextGenSize(n,j-n,R,k)
  }else if(g==2){
    out <- sum(pNextGenSize(n,1:(j-n-1),R,k) * pNextGenSize(1:(j-n-1),(j-n-1):1,R,k))
  }else{
    
    rs1 <- (j-n-g+1):1
    x1 <- rep(1:(j-n-g+1),choose(rs1+g-3,g-2))
    
    x <- matrix(0,length(x1),g-1)
    x[,1] <- x1
    
    pProd <- pNextGenSize(n,x1,R,k)
    
    rsA <- rs1
    for(i in 2:(g-1)){
      rsB <- sequence(rsA,rsA,-1)
      x[,i] <- rep(sequence(rsA),choose(rsB+g-2-i,g-1-i))
      pProd <- pProd * pNextGenSize(x[,i-1],x[,i],R,k)
      rsA <- rsB
    }
    xLast <- j-n-rowSums(x)
    pProd <- pProd * pNextGenSize(x[,g-1],xLast,R,k)
    out <- sum(pProd)
  }
  out
}