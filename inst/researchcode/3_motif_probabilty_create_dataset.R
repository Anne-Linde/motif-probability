# Tabula rasa
rm(list=ls())
gc()

# Load required libraries
library(dplyr)
library(TSEntropies)

# Load package functions
## TO DO: install package!
rDir <- "/Users/annelindelettink/Documents/Work MacBook Pro Annelinde/My Little Moves (MLM)/Sequence mapping/Physical behavior patterns/motif-probability/R/"
setwd(rDir)
source("durationProb.R")
source("emissionProb.R")
source("motifProb.R")
source("createMotif.R")
source("calculateMotifProbabilities.R")
source("lempel_ziv_complexity.R")

### User input required:
dataDir <- "/Users/annelindelettink/Documents/Work MacBook Pro Annelinde/My Little Moves (MLM)/Sequence mapping/Physical behavior patterns/sequences_old/all"
#modelDir <- paste0(dataDir, "/trained_models/")
#accelerationSequenceDir <- paste0(dataDir, "/segmented_data/")

#Load antrophometric data and GGIR day summary
spss <- read.csv("/Users/annelindelettink/GECKO/deel1/anthro_acc_measurement.csv")
ggir <- read.csv("/Users/annelindelettink/GECKO/preprocessing/part2_summary.csv")

# Define acceleration ranges for PA and SB
Amin.PA <- 2
Amax.PA <- 9
Amin.SB <- 0
Amax.SB <- 2
# Define epoch length
epoch_length <- 15 # hsmms were trained based on 15-sec epochs
n_epochs_min <- 60 / epoch_length

# Define motifs
motifs <- c( "motif.PA_30min",  "motif.PA_15min", "motif.PA_10min", "motif.PA_7min", "motif.PA_5min", "motif.PA_1min", "motif.SB_5min", "motif.SB_10min", "motif.SB_30min", "motif.SB_60min")
motif.PA_30min <- createMotif(Amin.PA, Amax.PA, 30, n_epochs_min)
motif.PA_15min <- createMotif(Amin.PA, Amax.PA, 15, n_epochs_min)
motif.PA_10min <- createMotif(Amin.PA, Amax.PA, 10, n_epochs_min)
motif.PA_7min <- createMotif(Amin.PA, Amax.PA, 7, n_epochs_min)
motif.PA_5min <- createMotif(Amin.PA, Amax.PA, 5, n_epochs_min)
motif.PA_1min <- createMotif(Amin.PA, Amax.PA, 1, n_epochs_min)
motif.SB_5min <- createMotif(Amin.SB, Amax.SB, 5, n_epochs_min)
motif.SB_10min <- createMotif(Amin.SB, Amax.SB, 10, n_epochs_min)
motif.SB_30min <- createMotif(Amin.SB, Amax.SB, 30, n_epochs_min)
motif.SB_60min <- createMotif(Amin.SB, Amax.SB, 60, n_epochs_min)

######################
# Calculate the probabilties for all individuals: 
probabilities <- calculateMotifProbabilities(motifs, dataDir)

# Merge anthropometrics with probabilities
colnames(probabilities)[colnames(probabilities) == "id"] <- "GKNO"
merged_data <- merge(probabilities, spss[,2:10], by = "GKNO", all = FALSE)  # inner join
rm(probabilities, spss)

# Add useful accelerometer info to dataframe

#match id's
GKNO <- c()
for(r in 1:nrow(ggir)){
  GKNO <- c(GKNO, strsplit(strsplit(ggir$ID, split = " ")[[r]][1], split = "_")[[1]][1])
}
ggir$GKNO <- GKNO
matching_rows <- which(ggir$GKNO %in% merged_data$GKNO)
ggir <- ggir[matching_rows,]

#AD_mean_HFENplus_mg_1-6am # average acceleration between 1-6 am
#AD_mean_HFENplus_mg_0-24hr # average daily acceleration 
#AD_L5_HFENplus_mg_0-24hr # average acceleration value during least active 5 hours
#AD_M5_HFENplus_mg_0-24hr # average acceleration value during most active 5 hours
#AD_MVPA_E15S_T890_HFENplus_0-24hr # Time spent in moderate-to-vigorous based on 15 second epoch size and an HFENplus metric
#AD_MVPA_E15S_T890_NeishabouriCount_y_0-24hr # Time spent in moderate-to-vigorous based on 15 second epoch size and Sirard cut-point

# Select specific columns from the data frame
selected_acc_data <- ggir %>%
  select(
    AD_mean_HFENplus_mg_1.6am,  # average acceleration between 1-6 am
    AD_mean_HFENplus_mg_0.24hr, # average daily acceleration 
    AD_L5_HFENplus_mg_0.24hr,   # average acceleration value during least active 5 hours
    AD_M5_HFENplus_mg_0.24hr,   # average acceleration value during most active 5 hours
    AD_MVPA_E15S_T890_HFENplus_0.24hr,  # Time spent in moderate-to-vigorous based on 15 second epoch size and an HFENplus metric
    AD_MVPA_E15S_T890_NeishabouriCount_y_0.24hr  # Time spent in moderate-to-vigorous based on 15 second epoch size and Sirard cut-point
  )

merged_data <- cbind(merged_data, selected_acc_data)

# Calculate complexity metrics
lz <- c()
sampleEntropy <- c()

filelist <- list.files(dataDir)
for(pp in 1:length(filelist)){
  load(paste0(dataDir, "/", filelist[pp]))# Load the sequences
  lz <- c(lz, lempel_ziv_complexity(hsmms$yhat))
  sampleEntropy <- c(sampleEntropy, TSEntropies::SampEn(hsmms$yhat))
}
merged_data$sample_entropy <- sampleEntropy
merged_data$lempel_ziv <- lz

# Ensure correct variable types
merged_data$motif.PA_30min <- as.numeric(as.character(merged_data$motif.PA_30min))
merged_data <- merged_data[order(merged_data$gender), ]

save(merged_data, file = paste0(dataDir, "/merged_data.RData"))
