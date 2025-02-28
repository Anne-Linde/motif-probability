#' durationProbability
#'
#' @description `durationProbability` Calculates the probability of observing a specific duration for a (shifted) Poisson distribution.
#'
#' @param duration The duration for which the probability needs to be calculated.
#' @param shift The shift parameter of the Poisson distribution.
#' @param lambda The lambda parameter of the Poisson distribution.
#'
#' @return The probability of the given duration for the (shifted) Poisson distribution.
#'
#' @examples
#' Calculate the probability of duration of 30 sec (model was trained on 15-sec epoch data) in a Poisson distribution with lambda = 0.5
#' length = 30/15
#' #durationProbability(duration = length, shift = 1, lambda = 0.5)
#' 
#' @importFrom stats dpois
#' @export

durationProbability <- function(duration, shift, lambda){
  if((duration - min(shift)) < 0){
    return(0)
  } else{
    return(dpois(duration-shift, lambda))
  }
}