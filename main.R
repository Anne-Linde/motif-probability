# This script illustrates how the functions in this directory can be applied

devtools::install_github("Anne-Linde/motif-probability", force = TRUE)
library(motifprobability)

### 1 ) Define motif (i.e. a sequence of bouts, each characterized by their own acceleration range and duration)
n_epochs_min <- 60 / 15 # hsmms were trained based on 15-sec epochs
motifA <- defineMotif(Amin = c(3, 0), Amax = c(4, 1.5), duration = c(30, 30), n_epochs_min) # Longer uninterrupted periods of high en low intensity bouts
motifB <- defineMotif(Amin = rep(c(3, 0), 5), Amax = rep(c(4, 1.5), 5), duration = rep(5, 10), n_epochs_min) # Frequent alternation of shorter high en low intensity bouts

### 2) Apply forward algorithm for motif probability calculation
motifs <- c( "motifA",  "motifB")
hsmmDir <- "/Users/annelindelettink/GECKO/preprocessing/manuscript/subset/output_rawinput_5years_complete_anthro/epochdata/validdata/mhsmmdata/models" # Contains the fitted HSMMs (obtained by running script: inst/1_preprocess_segment_accelerometer_data.R)

probabilities <- calculateMotifProbabilities(motifs, hsmmDir)