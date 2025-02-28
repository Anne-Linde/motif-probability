#' defineMotif
#'
#' @description `defineMotif` Generates a data frame for a user-specified sequence of bouts, with each rows representing a bout. Each bout is characterized by a range of accelerations (i.e., minimum and maximum acceleration value) and its duration.
#' 
#' @param Amin A vector containing the minimum acceleration values for the bouts in the sequence.
#' @param Amax A vector containing the maximum acceleration values for the bouts in the sequence.
#' @param duration A vector containing the durations in minutes for the bouts in the sequence.
#' @param nEpochsMin The number of epochs per minute (on which the fitted HSMM was trained) to consider when calculating the duration of the motif.
#'
#' @return A data frame with columns:
#' \itemize{
#'   \item \code{Amin}: The minimum acceleration.
#'   \item \code{Amax}: The maximum acceleration.
#'   \item \code{length}: The total number of epochs of the motif to represent the duration in minutes, calculated as \code{duration * nEpochsMin}.
#' }
#'
#' @examples
#' 
#' # Define a motif as the following sequence of bouts, for accelerometer data with epoch length = 15 sec:
#' # bout 1: Amin = 0, Amax = 2, duration = 5
#' # bout 2: Amin = 2, Amax = 10, duration = 1
#' 
#' Amin <- c(0, 2)
#' Amax <- c(2, 10)
#' duration <- c(5, 1)
#' nEpochsMin <- 60 / 15
#' # Create the motif
#' motif <- defineMotif(Amin, Amax, duration, nEpochsMin)
#' print(motif)
#' 
#' @export
defineMotif <- function(Amin, Amax, duration, nEpochsMin) {
  return(data.frame(Amin = Amin, Amax = Amax, length = duration * nEpochsMin))
}