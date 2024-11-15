#' Probability that one initial case leads to an outbreak lasting
#' less than g generations of transmission.
#'
#' @param gMax maximum number of generations.
#' @param R reproduction number.
#' @param k dispersion parameter.
#' @author Damon Toth
#' @returns A vector of probabilities for each number of generations from 1 to gmax
#' @export
pGen <- function(gMax,R,k){
  pgl <- rep(0,gMax)
  pgl[1] <- (1+R/k)^(-k)
  if(gMax > 1) for(g in 2:gMax) pgl[g] <- (1+R/k*(1-pgl[g-1]))^(-k)
  pgl
}
