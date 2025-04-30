#' Probability of final outbreak size with offspring distribution parameters switched after
#' two generations of transmission
#'
#' pFinalSizeSwitch2 is the probability that n initial cases lead to an extinguished
#' outbreak of total size j after any number of transmission generations (j includes 
#' the n initial cases), for a branching process with offspring distribution parameters
#' (R0,k0) during the first 2 transmission generations and (Rc,kc) during subsequent 
#' generations
#'
#' @param n Number of initial cases in generation 0
#' @param j Total outbreak size (>= n).
#' @param R0 Mean of negative binomial offspring distribution before the switch
#' @param k0 Dispersion of negative binomial offspring distribution before the switch
#' @param Rc Mean of negative binomial offspring distribution after the switch
#' @param kc Dispersion of negative binomial offspring distribution before the switch
#' @returns The probability of the final outbreak size
#' @examples
#' #With 5 initial cases, the probability that the final outbreak size is 20
#' #(including the initial 5):
#' pFinalSizeSwitch2(n=5, j=20, R0=2, k0=0.1, Rc=0.2, kc=0.1)
#' @export
pFinalSizeSwitch2 <- function(n,j,R0,k0,Rc,kc){
  if(j==n){
    out <- pNextGenSize(n,0,R0,k0)
  }else{
    x <- 1:(j-n-1)
    xr <- rep(x,rev(x))
    yr <- sequence(rev(x))
    out <- pFinalSizeAndGen(1,n,j,R0,k0) + sum(pNextGenSize(n,x,R0,k0)[xr] * pNextGenSize(xr,yr,R0,k0) * pFinalSize(yr,j-n-xr,Rc,kc))
  }
  out
}
