#' Probability that one initial case leads to an outbreak lasting
#' less than g generations of transmission, with offspring distribution 
#' parameters switched after generation two.
#'
#' @param gMax Maximum number of generations
#' @param R0 Basic reproduction number: mean of negative binomial offspring distribution from generation one and two
#' @param k0 Dispersion of negative binomial offspring distribution from generation one and two
#' @param Rc Control reproduction number: mean of negative binomial offspring distribution from generation three plus
#' @param kc Dispersion of negative binomial offspring distribution from generation three plus
#' @returns A vector of probabilities for each number of generations from 1 to gmax
#' @author Damon Toth
#' @examples
#' # Probability of outbreak lasting less than 1,2,3,...,10 generations:
#' pGenSwitch2(gMax=10, R0=3, k0=0.1, Rc=0.5, kc=1)
#' @export
pGenSwitch2 <- function(gMax,R0,k0,Rc,kc){
  pgl <- rep(0,gMax)
  pgl[1] <- ifelse(k0 < Inf, (1+R0/k0)^(-k0), exp(-R0))
  if(gMax > 1){
    pgl[2] <- ifelse(k0 < Inf, (1+R0/k0*(1-pgl[1]))^(-k0), exp(-R0*(1-pgl[1])))
  }
  if(gMax > 2){
    for(g in 3:gMax){
      pgl[g] <- ifelse(kc < Inf, (1+Rc/kc*(1-pgl[g-1]))^(-kc), exp(-Rc*(1-pgl[g-1])))
    }
  }
  pgl
}