#' accelerationProbability
#'
#' @description 'accelerationProbability' Calculates the probability within a specified acceleration range for a Gaussian acceleration distribution
#'
#' @param lower The lower bound of the range.
#' @param upper The upper bound of the range.
#' @param params A list containing the mean ('mu') and variance ('sigma') of the Gaussian distribution.
#'
#' @return The probability within the specified range for the Gaussian distribution.
#'
#' @examples
#' # Calculate the probability of observing a value within the acceleration range [1, 3] for a Gaussian distribution with mean = 2 and variance = 1
#' accelerationProbability(upper = 1, lower = 3, params = list(mu = 2, sigma = 1))
#' 
#' @export

accelerationProbability <- function(lower, upper, params){
  sd = sqrt(params$sigma)   # Calculate the standard deviation (because params$sigma is variance)
  
  if (lower > 0 & is.finite(upper)){
    return(pnorm(upper, mean=params$mu, sd=sd) - 
             pnorm(lower, mean=params$mu, sd=sd))
  } 
  if (lower == 0 & is.finite(upper)){
    return(pnorm(upper, mean=params$mu, sd=sd) -
             pnorm(-Inf,  mean=params$mu, sd=sd))
  }	
  if (lower > 0 & !is.finite(upper)){
    return(pnorm(Inf,   mean=params$mu, sd=sd) - 
             pnorm(lower, mean=params$mu, sd=sd))
  }		
  if (lower == 0 & !is.finite(upper)){
    return(pnorm( Inf, mean=params$mu, sd=sd) - 
             pnorm(-Inf, mean=params$mu, sd=sd))
  }			
}