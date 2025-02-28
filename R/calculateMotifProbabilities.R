#' calculateMotifProbabilities
#'
#' This function calculates the probabilities of specified motifs for each participant 
#' using fitted Hidden Semi-Markov Models (HSMMs). It loads HSMMs from specified files 
#' and computes motif probabilities for each participant based on their respective model.
#' @description calculateMotifProbabilities Calculates the probabilities of the specified motifs for each individual using the individually fitted Hidden Semi-Markov Models (HSMMs). For each file in the specified folder, it loads the HSMM and computes the specified motif probabilities.
#'
#' @param motifs A character vector containing the names of motifs for which the motif probability will be calculated. Each name should correspond to an object in the current R environment.
#' @param filepath A character string specifying the directory path where the fitted HSMM objects are located. The function will calculate the motif probabilities for all files in this directory.
#'
#' @return A data frame containing the following columns:
#' \describe{
#'   \item{id}{The first column corresponds to the filename}
#'   \item{motif probabilities}{Columns corresponding to each motif, containing the calculated probabilities.}
#' }
#' @examples
#' motifs_list <- c("motif1", "motif2")
#' motif_probabilities <- calculateMotifProbabilities(motifs = motifs_list, filepath = "path/to/hsmm/files")
#' 
#' @export

calculateMotifProbabilities <- function(motifs, filepath) {
  filelist <- list.files(filepath)
  probabilities <- matrix(NA, nrow = length(filelist), ncol = length(motifs) + 1)
  # Load fitted hsmms
  for (pp in 1:length(filelist)) {
    load(paste0(filepath, "/", filelist[pp]))# Load the fitted hsmm
    # Calculate motif probabilities for file
    for (motif_idx in 1:length(motifs)) {
      motif <- get(motifs[motif_idx])
      
      # Calculate motif probability and sum probabilities for each length
      ps <- exp(motifProbability(motif, hsmms))
      
      # Store the probability for this motif
      probabilities[pp, motif_idx + 1] <- ps
    }
    
  }
  probabilities <- cbind(filelist[pp], probabilities)
  colnames(probabilities) <- c("id", motifs)
  probabilities <- as.data.frame(probabilities)
  probabilities[, 2:length(motifs)+1] <- lapply(probabilities[, 2:length(motifs)+1], function(x) as.numeric(as.character(x))) # Convert motif probabilities to numeric
  
  return(probabilities)
}
