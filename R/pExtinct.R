#' Probability that one initial case leads to an outbreak that eventually dies out (stochastic extinction)
#'
#' @param R Reproduction number: mean of negative binomial offspring distribution
#' @param k Dispersion parameter of negative binomial offspring distribution
#' @author Damon Toth
#' @returns The probability of outbreak extinction
#' @examples
#' # Probability that a single case leads to an extinct outbreak when the offspring distribution
#' # is negative binomial with mean 2 and dispersion parameter 0.5:
#' pExtinct(R=2,k=0.5)
#' @export
pExtinct <- function(R,k) ifelse(R > 1, uniroot(function(q) q - (1+R/k*(1-q))^(-k), c(0,0.999999), tol=1e-7)$root, 1)