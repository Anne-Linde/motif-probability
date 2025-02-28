# This script illustrates how the functions in this directory can be applied

# Load all scripts
# This can be replaced by the package install later on...
rDir <- "/Users/annelindelettink/Documents/Work MacBook Pro Annelinde/My Little Moves (MLM)/Sequence mapping/Physical behavior patterns/motif-probability/R/"
setwd(rDir)
source("durationProbability.R")
source("accelerationProbability.R")
source("motifProbability.R")
source("defineMotif.R")
source("calculateMotifProbabilities.R")


### 1 ) Define motif (i.e. a sequence of bouts, each characterized by their own acceleration range and duration)
n_epochs_min <- 60 / 15 # hsmms were trained based on 15-sec epochs
motif.A <- defineMotif(Amin = c(3, 0), Amax = c(4, 1.5), duration = c(30, 30), n_epochs_min) # Longer uninterrupted periods of high en low intensity bouts
motif.B <- defineMotif(Amin = rep(c(3, 0), 5), Amax = rep(c(4, 1.5), 5), duration = rep(5, 10), n_epochs_min) # Frequent alternation of shorter high en low intensity bouts

### 2) Apply forward algorithm for motif probability calculation
motifs <- c( "motif.A",  "motif.B")
hsmmDir <- "/Users/annelindelettink/GECKO/preprocessing/manuscript/subset/output_rawinput_5years_complete_anthro/epochdata/validdata/mhsmmdata/models" # Contains the fitted HSMMs (obtained by running script: inst/1_preprocess_segment_accelerometer_data.R)

probabilities <- calculateMotifProbabilities(motifs, hsmmDir)