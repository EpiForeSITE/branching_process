#' Probability that one initial case leads to an outbreak lasting
#' less than g generations of transmission.
#'
#' @param gMax Maximum number of generations.
#' @param R Reproduction number.
#' @param k Dispersion parameter.
#' @author Damon Toth
#' @returns A vector of probabilities for each number of generations from 1 to gmax
#' @examples
#' # Probability of outbreak lasting less than 1,2,3,...,10 generations:
#' pGen(gMax=10, R=0.9, k=0.1)
#' @export
pGen <- function(gMax,R,k){
  pgl <- rep(0,gMax)
  pgl[1] <- (1+R/k)^(-k)
  if(gMax > 1) for(g in 2:gMax) pgl[g] <- (1+R/k*(1-pgl[g-1]))^(-k)
  pgl
}
