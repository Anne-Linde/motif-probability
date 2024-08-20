#' createMotif
#'
#' @description 'createMotif' generates a data frame representing a motif based on input parameters for the minimum and maximum activity levels, 
#' the duration of the motif in minutes, and the number of epochs per minute.
#'
#' @param Amin The minimum acceleration for the motif.
#' @param Amax The maximum acceleration for the motif.
#' @param length The duration of the motif in minutes.
#' @param nEpochsMin The number of epochs per minute (on which the hsmm was trained) to consider when calculating the total length of the motif.
#'
#' @return A data frame with columns:
#' \itemize{
#'   \item \code{Amin}: The minimum acceleration.
#'   \item \code{Amax}: The maximum acceleration.
#'   \item \code{length}: The total number of epochs of the motif to represent the duration in minutes, calculated as \code{length * n_epochs_min}.
#' }
#'
#' @examples
#' Create a motif with a minimum acceleration of 0, acceleration of 8, motif duration of 5 minutes, and an epoch length of 10 sec.
#' #createMotif(Amin = 0, Amax = 8, length = 5, nEpochsMin = 60/10)
#' 
#' @export
createMotif <- function(Amin, Amax, length, nEpochsMin) {
  return(data.frame(Amin = Amin, Amax = Amax, length = length * nEpochsMin))
}