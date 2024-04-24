#' duration_probability
#'
#' @description 'duration_probability' Calculates the probability of a certain duration in a Poisson distribution, with a given duration and lambda parameter
#'
#' @param duration The duration for which the probability needs to be calculated.
#' @param lambda The parameter of the Poisson distribution.
#'
#' @return The probability of the given duration in the Poisson distribution.
#'
#' @examples
#' Calculate the probability of duration of 30 sec (model was trained on 15-sec epoch data) in a Poisson distribution with lambda = 0.5
#' #duration_probability(duration = 2, lambda = 0.5)
#' 
#' @importFrom stats dpois
#' @export

duration_probability <- function(duration, lambda) {
  duration_probability <- dpois(duration, lambda)
  return(duration_probability)
}