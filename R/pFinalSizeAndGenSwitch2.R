#' Joint probability of outbreak final size and number of transmission generations with
#' offspring distribution parameters switched after generation two
#'
#' @param g Number of generations
#' @param n Number of initial cases
#' @param j Final size
#' @param R0 Mean of negative binomial offspring distribution from generation one and two
#' @param k0 Dispersion of negative binomial offspring distribution from generation one and two
#' @param Rc Mean of negative binomial offspring distribution from generation three on
#' @param kc Dispersion of negative binomial offspring distribution from generation three on
#' @author Damon Toth
#' @returns The joint probability of outbreak final size and number of transmission generations
#' @examples
#' # Probability that 1 initial infection leads to an outbreak of final size 20 over exactly
#' # 3 generations of transmission:
#' pFinalSizeAndGenSwitch2(g=3,n=1,j=20,R0=2,k0=0.1,Rc=0.5,kc=1)
#' @export
pFinalSizeAndGenSwitch2 <- function(g,n,j,R0,k0,Rc,kc){
  if(g==0){
    out <- pNextGenSize(n,0,R0,k0)
  }else if(g==1){
    out <- pNextGenSize(n,j-n,R0,k0)*pNextGenSize(j-n,0,R0,k0)
  }else if(g==2){
    out <- sum(pNextGenSize(n,1:(j-n-1),R0,k0) * pNextGenSize(1:(j-n-1),(j-n-1):1,R0,k0) * pNextGenSize((j-n-1):1,0,Rc,kc))
  }else{
    rs1 <- (j-n-g+1):1
    x1 <- rep(1:(j-n-g+1),choose(rs1+g-3,g-2))
    
    x <- matrix(0,length(x1),g-1)
    x[,1] <- x1
    
    pProd <- pNextGenSize(n,x1,R0,k0)
    
    rsA <- rs1
    rsB <- sequence(rsA,rsA,-1)
    x[,2] <- rep(sequence(rsA),choose(rsB+g-4,g-3))
    pProd <- pProd * pNextGenSize(x[,1],x[,2],R0,k0)

    if(g > 3){
      rsA <- rsB
      for(i in 3:(g-1)){
        rsB <- sequence(rsA,rsA,-1)
        x[,i] <- rep(sequence(rsA),choose(rsB+g-2-i,g-1-i))
        pProd <- pProd * pNextGenSize(x[,i-1],x[,i],Rc,kc)
        rsA <- rsB
      }
    }
    xLast <- j-n-rowSums(x)
    pProd <- pProd * pNextGenSize(x[,g-1],xLast,Rc,kc) * pNextGenSize(xLast,0,Rc,kc)
    out <- sum(pProd)
  }
  out
}