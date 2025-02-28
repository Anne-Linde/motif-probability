#' derive_initialparams
#' 
#' @description 'derive_initialparams' derives the initial parameters for the hsmm from the observations
#' 
#' @param data A 'hsmm.data' object that includes the accelerometer data.
#' @param n_states An integer indicating the number of states of the hsmm.
#' @param bout_durations A vector used to derive parameters for the sojourn duration (Poisson) (Default = c(1, 5, 10, 15)).
#' @param epoch_length An integer indicating the epoch length in seconds.
#' 
#' @return A list containing the following elements:
#'   \itemize{
#'     \item \code{initial_probability}: A vector of length \code{n_states} with the initial probabilities of the semi-Markov chain.
#'     \item \code{transition_matrix}: A transition matrix with the transition probabilities.
#'     \item \code{emission_distribution}: A list indicating the Gaussian observation distribution, consisting of \code{mu} and \code{sigma}.
#'     \item \code{sojourn_distribution}: A list with parameters for the Poisson runlength distribution.
#'   }
#' 
#' @details This function calculates initial parameters required for a hidden semi-Markov model (HSMM) based on the provided training data, number of states, bout durations, and epoch length.
#' 
#' @examples
#' # Define 'data', 'n_states', 'bout_durations', and 'epoch_length'
#' # result <- derive_initialparams(data, n_states = 3, bout_durations = c(60, 30, 15, 10), epoch_length = 5)
#' 
#' @importFrom stats quantile
#' @export

derive_initialparams <- function(data, n_states, bout_durations = c(60, 30, 15, 10), epoch_length){
  
  # Initial probabilities of the semi-Markov chain:
  init <- rep(1/n_states, n_states)
  
  # Transition probabilites: (Note: For two states, the matrix degenerates, taking 0 for the diagonal and 1 for the off-diagonal elements.)
  trans <- matrix(runif(n_states ^ 2), nrow = n_states) # transition matrix
  diag(trans) <- rep(0, n_states) # diagonal of the transition matrix must be 0
  for (j in 1:n_states) {
    trans[j, ] <- trans[j, ] / sum(trans[j, ])
  }
  
  # Observation/emission distribution (assume Gaussian distribution)
  # estimate means from data
  if (any(is.na(data))) {
    m0 <- seq(quantile(data, 0.05, na.rm = TRUE), 
              quantile(data, 0.95, na.rm = TRUE), length.out = n_states)
  } else {
    m0 <- seq(quantile(data, 0.05), quantile(data, 0.95), length.out = n_states)
  }
  #m0 <- c(0.75, 1 + 1:(n_states-1) * 0.5)^2
  #s0 <- sqrt(var0) 
  s0 <- 10 / (1 : n_states + 10)
  emis <- list(mu = m0, sigma = s0) 
  
  # Runlength/sojourn distibution
  soj_dur <- bout_durations * 60 / epoch_length # number of epochs for bout durations
  
  soj_poisson <- list(lambda = seq(soj_dur[1], soj_dur[length(n_states)], length.out = n_states),
                    shift = rep(1, n_states),
                   type = "poisson")
  # soj_gamma <- list(shape = rep((mean(data, na.rm = TRUE)^2)/var(data, na.rm = TRUE), n_states),
  #                   scale = rep((mean(data, na.rm = TRUE)^2)/var(data, na.rm = TRUE), n_states),
  #                   type = "gamma")
  
  return(list(initial_probability = init, transition_matrix = trans, emission_distribution = emis, 
              sojourn_distribution = soj_poisson))
}
