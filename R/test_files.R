# Tabula rasa
rm(list=ls())
gc()

library(mhsmm)

# user input required:
## TO DO: add test data file to repository...
datadir <- "/Users/annelindelettink/Documents/Work MacBook Pro Annelinde/My Little Moves (MLM)/Sequence mapping/Physical behavior patterns"
sequencedir <- paste0(datadir, "/sequences/all")

#datadir <- "/home/wessel/Research/AnnelindeMotif/Data/"
#Rdir    <- "/home/wessel/Research/AnnelindeMotif/Motif probability/"
#setwd(Rdir)
source("R/createMotif.R")
source("R/durationProb.R")
source("R/emissionProb.R")
source("R/motifProb.R")

# load a fitted HSMM
#filelist <- list.files(datadir, pattern = ".RData")
#load(paste0(datadir, "/", filelist[2]))
load(paste(sequencedir, list.files(sequencedir, pattern = ".RData")[1], sep = "/"))

############################################################
# Check 1 : Basic Probability Sums to 1
############################################################
# Situation: Ensure that the sum of motif probabilities with complete acceleration range over time is 1 
epochMin = 60/15
motif <- createMotif(Amin = 0, Amax = 1000, length = 1, nEpochsMin = epochMin)
ps <- 0 
for (t in 1:1000){ # to timestep 1000 to ensure that the motif is over the complete length
  motif$length <- t
  ps <- ps + exp(motifProb(motif, hsmms))
}
ps 

#############################################################
# Check 2: Splitting Probabilities over 1 Time Step 
#############################################################
#Situation: Ensure that splitting the acceleration range over one time step results in probabilities that sum correctly.

## Using one cutoff
cutoff     <- 2
Tlength    <- 1
motifL     <- motifU <- motifA   <- motif
motifL[1,] <- c(0,      cutoff, Tlength)
motifU[1,] <- c(cutoff,   1000, Tlength)
motifA[1,] <- c(0,        1000, Tlength)

# Summing the upper and lower cutoff acceleration motif probability should be very close to the motif probability of the complete acceleration range (e.g. within 1e-05)
split1     <- exp(motifProb(motifL, hsmms)) + exp(motifProb(motifU, hsmms))
ps1        <- exp(motifProb(motifA, hsmms))
split1 - ps1
abs(split1 - ps1) < 1e-5

# Using multiple cutoffs
cutoff     <- 2
Tlength    <- 1
motifL     <- motifM <- motifU <- motifA   <- motif
motifL[1,] <- c(0,        cutoff,   Tlength)
motifM[1,] <- c(cutoff,   2*cutoff, Tlength)
motifU[1,] <- c(2*cutoff, 1000,     Tlength)
motifA[1,] <- c(0,        1000,     Tlength)

# Summing the upper, middle, and lower cutoff acceleration motif probability should be very close to the motif probability of the complete acceleration range (e.g. within 1e-05)
split2     <- exp(motifProb(motifL, hsmms)) + exp(motifProb(motifM, hsmms)) + exp(motifProb(motifU, hsmms))
split2 - ps1
abs(split2 - ps1) < 1e-5
