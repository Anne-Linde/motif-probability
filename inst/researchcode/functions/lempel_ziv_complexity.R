# Function to calculate Lempel-Ziv complexity
lempel_ziv_complexity <- function(seq) {
  # Convert sequence to a string
  seq <- as.character(seq)
  n <- length(seq)
  if (n == 0) return(0)
  
  # Initialize variables
  lz_complexity <- 0
  substring_set <- character()
  
  # Loop through the sequence
  i <- 1
  while (i <= n) {
    for (j in i:n) {
      current_substring <- paste(seq[i:j], collapse = "")
      
      if (!(current_substring %in% substring_set)) {
        substring_set <- c(substring_set, current_substring)
        lz_complexity <- lz_complexity + 1
        break
      }
    }
    i <- j + 1
  }
  
  return(lz_complexity)
}