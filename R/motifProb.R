#' motifProb
#'
#' @description 'motifProb' Calculates the probability of a motif within a hidden semi-Markov model (hsmm) using the forward algorithm. It requires a motif and the hsmm as inputs.
#'
#' @param motif A motif is defined as a sequence of activities characterized by their acceleration range and lengths. A data.frame object where the rows represent the activity, with columns Amin, Amax, and length.
#' @param hsmm An object containing the Hidden Semi-Markov Model (HSMM) parameters of class 'hsmm'.
#'
#' @return The logarithm of the probability of the entire motif within the hsmm
#'
#' @examples
#' # Define 'motif' and 'hsmm' objects
#' # motif <- data.frame(Amin = c(1, 2, 3), Amax = c(4, 5, 6), length = c(2, 3, 4), upper_bound = c(10, 25, 30))
#' # hsmm <- mhsmm::hsmmfit()
#' # Calculate motif probability
#' # motifProb(motif, hsmm)
#'
#' @importFrom stats pnorm
#' @export

motifProb <- function(motif, hsmm){  
  ### Define the hsmm parameters required to calculate the motif probability
  iProbs  <- hsmm$model$init           # Initial state probabilities
  Pmatrix <- hsmm$model$transition     # Transition matrix
  Bparams <- hsmm$model$parms.emission # Emission distribution parameters (state acceleration means and variance)
  lambda  <- hsmm$model$sojourn$lambda # Lambda (rate parameter for Poisson distribution)
  shift   <- hsmm$model$sojourn$shift  # Shift parameter for shifted Poisson distribution
  
  ### Initialize forward probabilities
  nStates <- nrow(Pmatrix)
  nEvents <- nrow(motif)
  alpha   <- matrix(-Inf, nrow = nStates, ncol = nEvents)
  
  # Initialize the first column of alpha incorporating emission and duration probabilities
  eProb     <- emissionProb(motif[1,]$Amin, motif[1,]$Amax, Bparams)
  dProb     <- durationProb(motif[1,]$length, shift, lambda)
  alpha[,1] <- log(iProbs) + motif[1,]$length * log(eProb) + log(dProb)
  
  ### Forward algorithm recursion
  if(nEvents > 1){ 
    # if motif has multiple events
    for (t in 2:nEvents) { # For each time step
      eProb <- emissionProb(motif[t,]$Amin, motif[t,]$Amax, Bparams)
      dProb <- durationProb(motif[t,]$length, shift, lambda)
      
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
  logMotifProb <- log(sum(exp(alpha[, nEvents])))
  
  #logMotifProb <- log(sum(exp(alpha[, nEvents]), na.rm = T))
  return(logMotifProb)
}
