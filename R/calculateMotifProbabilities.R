#' calculateMotifProbabilities
#'
#' This function calculates the probabilities of specified motifs for each participant 
#' using fitted Hidden Semi-Markov Models (HSMMs). It loads HSMMs from specified files 
#' and computes motif probabilities for each participant based on their respective model.
#'
#' @param motifs A character vector containing the names of motifs for which probabilities 
#'               will be calculated. Each name should correspond to an object in the 
#'               current environment.
#' @param filepath A character string specifying the directory path where the fitted HSMM 
#'                 files are located. The function will load all files in this directory.
#'
#' @return A data frame containing the following columns:
#' \describe{
#'   \item{id}{Participant ID extracted from the filename.}
#'   \item{motifs}{Columns corresponding to each motif, containing the calculated probabilities.}
#' }
#' The first column of the data frame contains participant IDs, and subsequent columns contain 
#' the motif probabilities as numeric values.
#'
#' @examples
#' # Assuming 'motifs_list' contains names of motifs and the HSMM files are in 'hsmm_models/' directory
#' motif_probs <- calculateMotifProbabilities(motifs = motifs_list, filepath = "hsmm_models/")
#'
#' @export

calculateMotifProbabilities <- function(motifs, filepath) {
  #filelist......
  filelist <- list.files(filepath)
  probabilities <- matrix(NA, nrow = length(filelist), ncol = length(motifs) + 1)
  # Load fitted hsmms and extract participant IDs
  for (pp in 1:length(filelist)) {
    load(paste0(filepath, "/", filelist[pp]))# Load the fitted hsmm
    probabilities[pp, 1] <- strsplit(filelist[pp], "_")[[1]][1] # Extract participant ID
    
    # Calculate motif probabilities for each participant
    for (motif_idx in 1:length(motifs)) {
      motif <- get(motifs[motif_idx])
      
      # Calculate motif probability and sum probabilities for each length
      ps <- exp(motifProb(motif, hsmms))
      
      # Store the total probability for this motif
      probabilities[pp, motif_idx + 1] <- ps
    }
    
  }
  colnames(probabilities) <- c("id", motifs)
  probabilities <- as.data.frame(probabilities)
  probabilities[, 2:length(motifs)+1] <- lapply(probabilities[, 2:length(motifs)+1], function(x) as.numeric(as.character(x))) # Convert motif probabilities to numeric
  
  return(probabilities)
}
