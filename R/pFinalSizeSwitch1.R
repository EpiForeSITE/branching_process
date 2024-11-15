#' Probability of final outbreak size with offspring distribution parameters switched after generation one
#'
#' pFinalSizeSwitch1 is the probability that n initial cases lead to an extinguished
#' outbreak of total size j after any number of transmission generations (j includes 
#' the n initial cases)
#'
#' @param n Number of initial cases in generation 0
#' @param j Total outbreak size (>= n).
#' @param R0 Mean of negative binomial offspring distribution from generation one
#' @param k0 Dispersion of negative binomial offspring distribution from generation one
#' @param Rc Mean of negative binomial offspring distribution from generation two on
#' @param kc Dispersion of negative binomial offspring distribution from generation two on
#' @returns The probability of the final outbreak size
#' @examples
#' With 5 initial cases, the probability that the final outbreak size is 5 to 20
#' (including the initial 5):
#' pFinalSizeSwitch1(n=5, j=5:20, R0=2, k0=0.1, Rc=0.2, kc=0.1)
#' @export
pFinalSizeSwitch1 <- function(n,j,R0,k0,Rc,kc){
  ifelse(j==n, pNextGenSize(n,0,R0,k0), sum(pNextGenSize(n,1:(j-n),R0,k0) * pFinalSize(1:(j-n),j-n,Rc,kc)))
}
