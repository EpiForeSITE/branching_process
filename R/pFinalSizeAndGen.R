#' Joint probability of outbreak final size and number of transmission generations
#'
#' @param g Number of generations.
#' @param n Number of initial cases
#' @param j Final size
#' @param R Reproduction number
#' @param k Dispersion parameter
#' @returns The joint probability of the final outbreak size and number of transmission generations
#' @author Damon Toth
#' @export
pFinalSizeAndGen <- function(g,n,j,R,k){
  
  if(g==0){
    out <- pNextGenSize(n,0,R,k)
  }else if(g==1){
    out <- pNextGenSize(n,j-n,R,k)*pNextGenSize(j-n,0,R,k)
  }else if(g==2){
    out <- sum(pNextGenSize(n,1:(j-n-1),R,k) * pNextGenSize(1:(j-n-1),(j-n-1):1,R,k) * pNextGenSize((j-n-1):1,0,R,k))
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
    pProd <- pProd * pNextGenSize(x[,g-1],xLast,R,k) * pNextGenSize(xLast,0,R,k)
    out <- sum(pProd)
  }
  out
}