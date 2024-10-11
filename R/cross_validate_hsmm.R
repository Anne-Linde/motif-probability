#' Configuration Parameters for HSMM
#'
#' This section defines the configuration parameters for the Hidden Semi-Markov Model (HSMM) 
#' fitting process. Adjust these parameters as necessary based on the specifics of your data 
#' and the requirements of your analysis.
#' 
#' @param data A list in `mhsmmdata` format, which should contain:
#'   \itemize{
#'     \item{Y:}{A numeric vector of observations.}
#'     \item{N:}{A numeric vector indicating the length of each observation sequence.}
#'   }
#'   This format is essential for fitting the HSMM and performing cross-validation. 
#'   Ensure that the data is correctly structured to avoid errors during processing.
#'                 
#' @param k_folds An integer specifying the number of folds for cross-validation. 
#'                 The default value is set to 5. If the dataset contains fewer valid 
#'                 days than the specified number of folds, the number of folds 
#'                 will be adjusted to match the number of data points available. 
#'                 
#' @param nStates An integer representing the number of states in the HSMM. 
#'                 Default is set to 3. 
#'                 
#' @param qMax An integer specifying the maximum number of iterations for the HSMM fitting 
#'              process. Default is set to 15.
#'              
#' @param dMax An integer that determines the maximum state duration, calculated as the 
#'              number of epochs per hour. In this case, it is set to 240, which corresponds 
#'              to a maximum duration of 60 minutes based on the epoch length (in seconds). 
#'              This is in line with the methodology outlined in van Kuppevelt et al. (2019).
#'              
#' @param boutDurations A numeric vector indicating the durations of each bout for the HSMM. 
#'                      The first bout must have a non-zero duration as the parameter lambda 
#'                      requires it to be greater than zero. This example sets the bout durations 
#'                      to 30 seconds, 10 seconds, and 5 seconds.
#'
#' @examples
#' # Setting up the parameters for HSMM fitting
#' nStates <- 3
#' qMax <- 15
#' epoch <- 15 # Example epoch length in seconds
#' dMax <- (60 * 60) / epoch
#' boutDurations <- c(30, 10, 5)
#'
#' @export

cross_validate_hsmm <- function(data, k_folds = 5, nStates = 3, boutDurations = c(30, 10, 5), epoch, qMax = 15, dMax = 240) {
  
  # Set folds
  set.seed(123)  # For reproducibility
  if(length(data$N) < k_folds){ # If there are less than 5 days of valid data, use the available number of days
    k_folds = length(data$N)
  }
  folds <- sample(rep(1:k_folds, length.out = k_folds))
  
  # Best model parameters
  log_likelihoods <- numeric(k_folds)  # Store log-likelihoods for each fold
  best_model <- NULL
  best_log_likelihood <- -Inf
  
  # Loop through each fold
  for (i in 1:k_folds) {
    cat("Processing Fold", i, "\n")

    # Prepare training and test sets
    test_indices <- which(folds == i, arr.ind = TRUE)
    train_indices <- which(data$N != data$N[test_indices])
  
    if (test_indices == 1) { # test set
      test_data <- data$Y[1:data$N[test_indices]]
    } else {
      # General case for test_indices > 1
      test_data <- data$Y[data$N[test_indices - 1] + 1 : data$N[test_indices]]
    }
    
    train_data <- c() # training set
    for (td in 1:length(train_indices)) {
      if(train_indices[td] == 1){
        train_data <- c(train_data, data$Y[1:data$N[train_indices[td]]])
      } else{
        train_data <- c(train_data, data$Y[data$N[train_indices[td-1]]+1:data$N[train_indices[td]]])
      }
    }
    
    # Fit HSMM to training data
    initParams <- derive_initialparams(train_data, nStates, boutDurations, epoch)
    startmodel <- hsmmspec(
      initParams$initial,
      initParams$transition,
      initParams$emission,
      initParams$sojourn_distribution,
      dens.emis = dnorm.hsmm
    )
    
    hsmms <- hsmmfit(train_data, startmodel, mstep = mstep.norm, maxit = qMax, M = dMax)
    
    # Calculate log-likelihood of the test data under the fitted model
    predictions <- predict(hsmms, test_data)
    log_likelihoods[i] <- predictions$loglik
    
    # Update the best model id the current log-likelihood is better
    if(log_likelihoods[i] > best_log_likelihood){
      best_log_likelihood <- log_likelihoods[i]
      best_model <- hsmms
    }
  }
  
  # Return log-likelihoods and the best model, and the sequence with original data
  mhsmmdata$s <- predict(best_model, mhsmmdata$Y)$s
  mean_loglik <- mean(log_likelihoods)
  
  return(list(log_likelihoods = log_likelihoods, best_model = best_model, sequence_data = mhsmmdata))
}