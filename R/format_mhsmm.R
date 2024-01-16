#' format_mhsmm
#'
#' @description 'format_mhsmm' formats the valid data list into a hsmm.data object
#' 
#' @param data A list containing participant data with elements that include a valid day data.frame object containing metrics.
#' @param variable_name A string containing the name of the metric to be selected.
#' @return A 'hsmm.data' object consisting of:
#'   \itemize{
#'     \item \code{df}: A data.frame containing the observed sequence, i.e., all vectors with observations pasted.
#'     \item \code{N}: A vector indicating the length of each observed sequence.
#'   }
#' 
#' @details This function takes a list of participant data and a vector specifying metric names and constructs an 'hsmm.data' object. It extracts the specified metrics from each valid day data frame and concatenates them into a single data frame, also recording the length of each observed sequence.
#' 
#' @examples
#' # Define 'data' and 'variable_name'
#' # result_hsmm <- format_mhsmm(data, variable_name = "ENMO")
#' 
#' @export

format_mhsmm <- function(data, variable_name) {
  N <- c()
  df <- data.frame()
  # Vectors including the data observations:
  for (day in 1:length(data)) {
    col = which(colnames(data[[day]]) == variable_name)
    df <- rbind(df, data[[day]][col])
    N <- c(N, nrow(data[[day]][col]))
  }
  colnames(df) <- c("Y")
  object <- list(Y = df[["Y"]], N = N)
  class(object) <- "hsmm.data"
  return(object)
}