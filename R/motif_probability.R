#' motif_probability
#'
#' @description 'motif_probability' Calculates the probability of a motif within a hidden semi-Markov model (hsmm) using the forward algorithm. It requires a motif and the hsmm as inputs.
#'
#' @param motif A motif is defined as a sequence of activities characterized by their acceleration range and lengths. A data.frame object where the rows represent the activity, with columns Amin, Amax, and length.
#' @param hsmm An object containing the Hidden Semi-Markov Model (HSMM) parameters of class 'hsmm'.
#'
#' @return The logarithm of the probability of the entire motif within the hsmm
#'
#' @examples
#' # Define 'motif' and 'hsmm' objects
#' # motif <- data.frame(Amin = c(1, 2, 3), Amax = c(4, 5, 6), length = c(2, 3, 4))
#' # hsmm <- mhsmm::hsmmfit()
#' # Calculate motif probability
#' # motif_probability(motif, hsmm)
#'
#' @importFrom stats pnorm
#' @export

motif_probability <- function(motif, hsmm) {
  
  ### Define the hsmm parameters required to calculate the motif probability
  init_probs <- hsmm$model$init
  transition_matrix <- hsmm$model$transition # Transition matrix
  acceleration_params <- hsmm$model$parms.emission # Emission distribution parameters (state acceleration means and variance)
  lambda <- hsmm$model$sojourn$lambda # Lambda
  
  ### Initialize forward probabilities
  n_states <- nrow(transition_matrix)
  n_events <- nrow(motif)
  alpha <- matrix(-Inf, nrow = n_states, ncol = n_events)
  
  # Initialize the first column of alpha incorporating emission and duration probabilities
  emission_prob <- emission_probability_within_range(motif[1,]$Amin, motif[1,]$Amax, acceleration_params)
  duration_prob <- duration_probability(motif[1,]$length, lambda)
  alpha[, 1] <- log(init_probs) + log(emission_prob) + log(duration_prob)
  
  ### Forward algorithm recursion
  for (t in 2:n_events) { # For each time step
    emission_prob <- emission_probability_within_range(motif[t,]$Amin, motif[t,]$Amax, acceleration_params)
    duration_prob <- duration_probability(motif[t,]$length, lambda)
    
    for (j in 1:n_states) { # For each state
      # Calculate the sum for transitioning from all possible states to the current state (j) at the previous time step (t-1)
      sum_previous_alphas <- sum(exp(alpha[j, t-1] + log(transition_matrix[, j])))
      # Update the alpha value for state 'j' at time step 't'; calculate forward probability for that state (j) at time t considering the probabilities of transitioning from all previous states, multiplied by their respective sojourn distributions, and the emission probability of the current observation.
      alpha[j, t] <- log(sum_previous_alphas) + log(emission_prob[j]) + log(duration_prob[j])
    }
  }
  
  # Calculate the log probability of the entire motif (i.e. sum forward probabilities calculated at the last time point)
  log_motif_probability <- log(sum(exp(alpha[, n_events])))
  return(log_motif_probability)
}