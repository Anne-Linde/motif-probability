#' motifProbability
#'
#' @description 'motifProbability' Calculates the probability that a specific sequence of bouts (motif) occurs within the accelerometer data, based on the parameters of the fitted hidden semi-Markov model (hsmm) using the forward algorithm.
#'
#' @param motif A motif is defined as a sequence of bouts, each characterized by their acceleration range and duration. A data.frame object where each rows represent a bout, with columns Amin, Amax, and length.
#' @param hsmm An object containing the Hidden Semi-Markov Model (HSMM) parameters of class 'hsmm'.
#'
#' @return The logarithm of the motif probability 
#'
#' @examples
#' # Define 'motif' and 'hsmm' objects
#' # motif <- defineMotif(Amin = c(1, 2, 3), Amax = c(4, 5, 6), duration = c(1, 5, 10))
#' # hsmm <- mhsmm::hsmmfit()
#' # Calculate motif probability
#' # motifProbability(motif, hsmm)
#'
#' @importFrom stats pnorm
#' @export

motifProbability <- function(motif, hsmm){  
  ### Define the hsmm parameters required to calculate the motif probability
  iProbs  <- hsmm$model$init           # Initial state probabilities
  Pmatrix <- hsmm$model$transition     # Transition matrix
  Bparams <- hsmm$model$parms.emission # Emission distribution parameters (state acceleration means and variance)
  lambda  <- hsmm$model$sojourn$lambda # Lambda (rate parameter for Poisson distribution)
  shift   <- hsmm$model$sojourn$shift  # Shift parameter for shifted Poisson distribution
  
  ### Initialize forward probabilities
  nStates <- nrow(Pmatrix)
  nBouts <- nrow(motif)
  alpha   <- matrix(-Inf, nrow = nStates, ncol = nBouts)
  
  # Initialize the first column of alpha incorporating emission and duration probabilities
  eProb     <- accelerationProbability(motif[1,]$Amin, motif[1,]$Amax, Bparams)
  dProb     <- durationProbability(motif[1,]$length, shift, lambda)
  alpha[,1] <- log(iProbs) + motif[1,]$length * log(eProb) + log(dProb)
  
  ### Forward algorithm recursion
  if(nBouts > 1){ 
    # if motif has multiple bouts
    for (t in 2:nBouts) { # For each time step
      eProb <- accelerationProbability(motif[t,]$Amin, motif[t,]$Amax, Bparams)
      dProb <- durationProbability(motif[t,]$length, shift, lambda)
      
      for (j in 1:nStates) { # Iterate through each state of the fitted hsmm
        # Calculate the sum for transitioning from all possible states to the current state (j) at the previous time step (t-1)
        sumOldAlphas <- sum(exp(alpha[j, t-1] + log(Pmatrix[, j])))
        
        # Update the alpha value for state 'j' at time step 't'; calculate forward probability for that 
        # state (j) at time t considering the probabilities of transitioning from all previous states, 
        # multiplied by their respective sojourn distributions, and the emission probability of the current observation.
        alpha[j, t] <- log(sumOldAlphas) + motif[t,]$length * log(eProb[j]) + log(dProb[j])
      }
    }
  }
  
  # Calculate the log probability of the entire motif (i.e. sum forward probabilities calculated at the last time point)
  logMotifProb <- log(sum(exp(alpha[, nBouts])))
  
  #logMotifProb <- log(sum(exp(alpha[, nBouts]), na.rm = T))
  return(logMotifProb)
}
