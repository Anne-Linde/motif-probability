# Tabula rasa
rm(list=ls())
gc()

# user input required:
datadir <- "/Users/annelindelettink/Documents/Work MacBook Pro Annelinde/My Little Moves (MLM)/Sequence mapping/Physical behavior patterns/sequence-probability"
sequencedir <- paste0(datadir, "/data/sequences")

# Load all scripts
# This can be replaced by the package install later on...
my_functions_folder = "/Users/annelindelettink/Documents/Work MacBook Pro Annelinde/My Little Moves (MLM)/Sequence mapping/Physical behavior patterns/motif-probability/R"
for (function_file in dir(my_functions_folder, full.names = F)) source(function_file) #load functions

### 1 ) Define motif (i.e. a sequence of events, each characterised by their acceleration range and length)
# Article example input
Ab.odd <- c(3, 4)
Ab.even <- c(0, 1.5)
Lb.alternating <- rep(5, 10)
Lb.long <- rep(30, 2)
# Length from minutes to epoch level, lambda was based on 15-sec epoch time series
epoch_length = 15
n_epochs_min = 60/epoch_length
Lb.alternating <- Lb.alternating * n_epochs_min
Lb.long <- Lb.long * n_epochs_min

motif.alternating <- data.frame()
motif.long <- data.frame()

for(acc in 1:length(Lb.alternating)){
  if (acc %% 2 == 0) { # Even
    motif.alternating <- rbind(motif.alternating, c(3, 4, Lb.alternating[acc]))
  } else {
    motif.alternating <- rbind(motif.alternating, c(0, 1.5, Lb.alternating[acc]))
  }
}
for(acc in 1:length(Lb.long)){
  if (acc %% 2 == 0) { # Even
    motif.long <- rbind(motif.long, c(3, 4, Lb.long[acc]))
  } else {
    motif.long <- rbind(motif.long, c(0, 1.5, Lb.long[acc]))
  }
}
colnames(motif.alternating) <- c("Amin", "Amax", "length")
colnames(motif.long) <- c("Amin", "Amax", "length")

### 2) Apply forward algorithm for motif probability calculation

# For each participant load the fitted hsmm
filelist <- list.files(sequencedir, pattern = ".RData")

prob_long <- c()
prob_alternating <- c()

for(pp in 1:length(filelist)){
  load(paste0(sequencedir, "/", filelist[pp]))
  prob_long <- c(prob_long, motif_probability(motif.long, hsmms))
  prob_alternating <- c(prob_alternating, motif_probability(motif.alternating, hsmms))
}

