# Tabula rasa
rm(list=ls())
gc()

# User input required:
rDir <- "/Users/annelindelettink/Documents/Work MacBook Pro Annelinde/My Little Moves (MLM)/Sequence mapping/Physical behavior patterns/motif-probability/R/"
dataDir <- "/Users/annelindelettink/Documents/Work MacBook Pro Annelinde/My Little Moves (MLM)/Sequence mapping/Physical behavior patterns/hsmms/all/"
#datadir <- "/home/wessel/Research/AnnelindeMotif/Data/"
#Rdir    <- "/home/wessel/Research/AnnelindeMotif/Motif probability/"
# Load package functions
## TO DO: install package?
setwd(rDir)
source("durationProb.R")
source("emissionProb.R")
source("motifProb.R")
source("createMotif.R")

# Load fitted HSMM
## TO DO: add example fitted? as test file
load(paste0(dataDir, list.files(dataDir, pattern = ".RData")[1]))

############################################################
# Check 1 : Basic Probability Sums to 1
############################################################
# Situation: Ensure that the sum of motif probabilities with complete acceleration range over time is 1 

motif <- data.frame(0, 1000, 1)
colnames(motif) <- c("Amin", "Amax", "length")
ps <- 0 
for (t in 1:1000){
  motif[3] <- t
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

#############################################################
# Check 3: Splitting Probabilities over Multiple Time Steps
#############################################################
# Situation: Ensure that splitting the same motif over multiple time steps maintains the correct probability.

#############################################################
## Over 2 time steps
cutoff      <- 2
Tlength     <- 2
motifL2     <- motifU2 <- motifA2 <- motifL
motifL2[1,] <- c(0,      cutoff, Tlength)
motifU2[1,] <- c(cutoff,   1000, Tlength)
motifA2[1,] <- c(0,        1000, Tlength)

split2      <- exp(motifProb(motifL2, hsmms)) + exp(motifProb(motifU2, hsmms))
ps2         <- exp(motifProb(motifA2, hsmms))
split2 - ps2
abs(split2 - ps2) < 1e-5

# Situation: Ensure that split probabilities over two time steps maintains quite close correct probability
motifM2a     <- rbind(motifL2, motifU2)
motifM2b     <- rbind(motifU2, motifL2)
motifM2a[,3] <- motifM2b[,3] <- 1 
exp(motifProb(motifM2a, hsmms)) - exp(motifProb(motifM2b, hsmms))

motifL11     <- rbind(motifL2, motifL2)
motifU11     <- rbind(motifU2, motifU2)
motifL11[,3] <- motifU11[,3] <- 1
exp(motifProb(motifL11, hsmms)) - exp(motifProb(motifU11, hsmms))
# These should be different

ps2 - (exp(motifProb(motifM2a, hsmms)) + exp(motifProb(motifM2b, hsmms)) +
         exp(motifProb(motifL11, hsmms)) + exp(motifProb(motifU11, hsmms)) +
         exp(motifProb(motifL2,  hsmms)) + exp(motifProb(motifU2, hsmms)))

abs(ps2 - (exp(motifProb(motifM2a, hsmms)) + exp(motifProb(motifM2b, hsmms)) +
             exp(motifProb(motifL11, hsmms)) + exp(motifProb(motifU11, hsmms)) +
             exp(motifProb(motifL2,  hsmms)) + exp(motifProb(motifU2, hsmms)))) < 1e-5

exp(motifProb(motifL11, hsmms)) - exp(motifProb(rbind(motifL), hsmms))^2

#############################################################
## > 2 time steps

Tlength    <- 40 #doel 40
cutoff     <- 2 # zodra cutoff 5 of hoger is, wordt het verschil tolerabel (< 1e-5) 

# Determine unique step combinations (ignoring order)
stepSizes <- list()
stepSizes[[1]] <- rep(1, Tlength)
if (Tlength > 2){
  for (s in 2:(Tlength-1)){
    for (u in 1:100000){
      steps <- sample(1:Tlength, s, replace=TRUE)
      if (sum(steps) == Tlength){ 
        stepSizes[[length(stepSizes)+1]] <- sort(steps)
      }
    }
  }
}
stepSizes[[length(stepSizes)+1]] <- Tlength
stepSizes                        <- unique(stepSizes)

# Calculate frequency of the step sizes
nCombStep <- numeric()
for (u in 1:length(stepSizes)){
  slh          <- as.numeric(table(stepSizes[[u]]))
  nCombStep[u] <- factorial(sum(slh)) / prod(factorial(slh))
}

# Create motif list
motifList <- list()
for (u in 1:length(stepSizes)){
  steps <- stepSizes[[u]] 
  if (length(steps) > 1){
    LU      <- matrix(c(0, cutoff, cutoff, 1000), ncol=2, byrow=TRUE)
    LUcombi <- matrix(as.numeric(as.matrix(DoE.base::fac.design(2, length(steps)))), ncol=length(steps))
    for (w in 1:nrow(LUcombi)){
      motif                            <- data.frame(cbind(LU[LUcombi[w,],], steps, nCombStep[u]))
      colnames(motif)                  <- c("Amin", "Amax", "length", "freq")
      motifList[[length(motifList)+1]] <- motif
    }
  }
}
# Add last to motifs
motif                            <- data.frame(matrix(c(0,    cutoff, Tlength, 1), nrow=1))
colnames(motif)                  <- c("Amin", "Amax", "length", "freq")
motifList[[length(motifList)+1]] <- motif
motif                            <- data.frame(matrix(c(cutoff,  1000, Tlength, 1), nrow=1))
colnames(motif)                  <- c("Amin", "Amax", "length", "freq")
motifList[[length(motifList)+1]] <- motif

# Now calculate probabilities of all motifs in list
probTotal <- 0
for (m in 1:length(motifList)){
  probTotal <- probTotal + motifList[[m]][1,4] * exp(motifProb(motifList[[m]][,1:3], hsmms))
}

# Calculate check probability
motifAll           <- data.frame(matrix(c(0,    1000, Tlength), nrow=1))
colnames(motifAll) <- c("Amin", "Amax", "length")
exp(motifProb(motifAll, hsmms)) - probTotal
abs(exp(motifProb(motifAll, hsmms)) - probTotal) < 1e-5

#############################################################
# Check 4: Order of bouts
#############################################################

# The order does not affect the probabilities (bureau / wandelen - wandelen / bureau)
exp(motifProb(motifM2a, hsmms)) 
exp(motifProb(motifM2b, hsmms)) 

# But if the length of the motif is longer, the order should matter (bureau bureau / wandelen wandelen - bureau wandelen bureau wandelen)
exp(motifProb(rbind(motifM2a, motifM2a), hsmms)) 
exp(motifProb(rbind(motifL2, motifU2), hsmms)) 
# it does!

rbind(motifM2a, motifM2a)