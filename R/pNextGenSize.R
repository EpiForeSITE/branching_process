#' Probability of y total transmission directly from x independent infected individuals
#'
#' @param x number of infected individuals in generation n
#' @param y number of total transmissions in generation n+1
#' @param R mean of negative binomial offspring distribution
#' @param k dispersion of negative binomial offspring distribution
#' @returns The probability of the given number of transmissions
#' @examples
#' # With 5 individuals in this generation, what is the probability of
#' # 0 to 15 transmissions in the next generation?
#' pNextGenSize(x=5, y=0:15, R=0.2, k=0.1)
#' @export
pNextGenSize <- function(x,y,R,k) exp(lgamma(k*x+y)-lgamma(k*x)-lgamma(y+1)+y*log(R/k)-(k*x+y)*log(1+R/k))
