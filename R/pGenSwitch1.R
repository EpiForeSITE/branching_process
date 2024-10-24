#' Probability that one initial case leads to an outbreak lasting
#' less than g generations of transmission.
#'
#' @param gMax maximum number of generations
#' @param R0 basic reproduction number: mean of negative binomial offspring distribution from generation one
#' @param k0 dispersion of negative binomial offspring distribution from generation one
#' @param Rc control reproduction number: mean of negative binomial offspring distribution from generation two plus
#' @param kc dispersion of negative binomial offspring distribution from generation two plus
#' @author Damon Toth
#' @export
pGenSwitch1 <- function(gMax,R0,k0,Rc,kc){
  pgl <- rep(0,gMax)
  pgl[1] <- (1+R0/k0)^(-k0)
  if(gMax > 1) for(g in 2:gMax) pgl[g] <- (1+Rc/kc*(1-pgl[g-1]))^(-kc)
  pgl
}