#' Probability that one initial case leads to an outbreak that eventually dies out (stochastic extinction)
#'
#' @param R Reproduction number: mean of negative binomial offspring distribution
#' @param k Dispersion parameter of negative binomial offspring distribution
#' @author Damon Toth
#' @returns The probability of outbreak extinction
#' @export
pExtinct <- function(R,k) ifelse(R > 1, uniroot(function(q) q - (1+R/k*(1-q))^(-k), c(0,0.999999), tol=1e-7)$root, 1)