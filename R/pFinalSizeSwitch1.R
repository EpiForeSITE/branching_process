#' Probability of final outbreak size with offspring distribution parameters switched after generation one
#'
#' pFinalSizeSwitch1 is the probability that n initial cases lead to an extinguished
#' outbreak of total size j after any number of transmission generations (j includes 
#' the n initial cases)
#'
#' @param n number of initial cases in generation 0
#' @param j total outbreak size (>= n).
#' @param R0 mean of negative binomial offspring distribution from generation one
#' @param k0 dispersion of negative binomial offspring distribution from generation one
#' @param Rc mean of negative binomial offspring distribution from generation two on
#' @param kc dispersion of negative binomial offspring distribution from generation two on
#' @examples
#' pFinalSizeSwitch1(5, 5:20, R0=2, k0=0.1, Rc=0.2, kc=0.1)
#'
#' @export
pFinalSizeSwitch1 <- function(n,j,R0,k0,Rc,kc){
  ifelse(j==n, pNextGen(n,0,R0,k0), sum(pNextGen(n,1:(j-n),R0,k0) * pFinalSize(1:(j-n),j-n,Rc,kc)))
}
