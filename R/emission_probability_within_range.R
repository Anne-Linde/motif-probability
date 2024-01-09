#' emission_probability_within_range
#'
#' @description 'emission_probability_within_range' Calculates the probability within a specified range for a Gaussian emission distribution
#'
#' @param lower_bound The lower bound of the range.
#' @param upper_bound The upper bound of the range.
#' @param params A list containing the mean ('mu') and variance ('sigma') of the Gaussian distribution.
#'
#' @return The probability within the specified range for the Gaussian distribution.
#'
#' @examples
#' # Calculate the probability within range [1, 3] for a Gaussian distribution with mean = 2 and variance = 1
#' params <- list(mu = 2, sigma = 1)
#' emission_probability_within_range(1, 3, params)
#' 
#' @export

emission_probability_within_range <- function(lower_bound, upper_bound, params) {
  emission_prob_lower <- pnorm(lower_bound, mean = params$mu, sd = sqrt(params$sigma))
  emisiion_prob_upper <- pnorm(upper_bound, mean = params$mu, sd = sqrt(params$sigma))
  emission_prob_within_range <- emisiion_prob_upper - emission_prob_lower
  return(emission_prob_within_range)
}