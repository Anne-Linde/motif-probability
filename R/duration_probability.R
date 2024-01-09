#' duration_probability
#'
#' @description 'duration_probability' Calculates the probability of a certain duration in an exponential distribution, with a given duration and lambda parameter
#'
#' @param duration The duration for which the probability needs to be calculated.
#' @param lambda The lambda parameter of the exponential distribution.
#'
#' @return The probability of the given duration in the exponential distribution.
#'
#' @examples
#' Calculate the probability of duration 2 in an exponential distribution with lambda = 0.5
#' duration_probability
#' 
#' @export

duration_probability <- function(duration, lambda) {
  duration_probability <- exp(-lambda * duration)
  return(duration_probability)
}