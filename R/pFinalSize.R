#' Final outbreak size probability
#'
#' @param n Number of initial cases in generation 0
#' @param j Total outbreak size (>= n).
#' @param R Mean of negative binomial offspring distribution
#' @param k Dispersion of negative binomial offspring distribution
#' @returns The final size probability
#' @examples
#' # With 5 initial individuals and negative binomial offspring distribution with mean R=0.2
#' # and dispersion k=0.1, gives the probability of outbreak extinction with a total number
#' # final outbreak size of exactly 5 to 20 individuals (including the initial 5).
#' pFinalSize(5, 5:20, R=0.2, k=0.1)
#'
#' @export
pFinalSize <- function(n,j,R,k) pNextGenSize(j,j-n,R,k)*n/j
