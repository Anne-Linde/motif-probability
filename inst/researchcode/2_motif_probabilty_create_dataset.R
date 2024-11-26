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
dataDir <- "/Users/annelindelettink/GECKO/preprocessing/manuscript/subset/output_rawinput_5years_complete_anthro/epochdata/validdata/mhsmmdata/models"
#dataDir <- "/Users/annelindelettink/Documents/Work MacBook Pro Annelinde/My Little Moves (MLM)/Sequence mapping/Physical behavior patterns/sequences_old/all"
#modelDir <- paste0(dataDir, "/trained_models/")
#accelerationSequenceDir <- paste0(dataDir, "/segmented_data/")

#Load antrophometric data and GGIR day summary
spss <- read_sav("/Users/annelindelettink/GECKO/preprocessing/manuscript/subset/complete_spss_GECKO_5years.sav")
ggir <- read.csv("/Users/annelindelettink/GECKO/preprocessing/manuscript/subset/output_rawinput_5years_complete_anthro/results/part2_summary.csv")

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
colnames(probabilities)[colnames(probabilities) == "id"] <- "GKNO"

#### Select anthropometrics
who <- anthro::anthro_zscores(sex = as.numeric(spss$sex),
                              age = as.integer(spss$measure_age_217 / 30.4375),
                              weight = as.numeric(spss$GEWICHT_217/1000),
                              lenhei = as.numeric(spss$LENGTE_217))
spss_subset <- as.data.frame(cbind(spss$GKNO, gsub(".gt3x", ".csv.RData", spss$Filename), spss$sex, who$zbmi, (spss$measure_age_217/ 30.4375)/12, spss$SDQ_HYP_217, spss$BUIKOM_217_QC))
colnames(spss_subset) <- c("GKNO", "filename", "sex", "zbmi", "measure_age", "SDQ_hyp", "buikomvang")

# Merge anthropometrics with probabilities
merged_data <- merge(probabilities, spss_subset, by = "GKNO", all = FALSE)  # inner join

#### Add accelerometer volume data
#match id's
GKNO <- c()
for(r in 1:nrow(ggir)){
  GKNO <- c(GKNO, strsplit(strsplit(ggir$ID, split = " ")[[r]][1], split = "_")[[1]][1])
}
ggir$GKNO <- GKNO
matching_rows <- which(ggir$GKNO %in% merged_data$GKNO)
ggir <- ggir[matching_rows,]

# Select specific columns from the data frame
selected_acc_data <- ggir %>%
  select(
    AD_mean_HFENplus_mg_0.24hr, # average daily acceleration
    AD_L5_HFENplus_mg_0.24hr,   # average acceleration value during least active 5 hours
    AD_M5_HFENplus_mg_0.24hr,   # average acceleration value during most active 5 hours
    AD_mean_NeishabouriCount_vm_mg_0.24hr, # average daily VM count
    #AD_MVPA_E15S_T890_HFENplus_0.24hr,  # Time spent in moderate-to-vigorous based on 15 second epoch size and an HFENplus metric
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
merged_data <- merged_data[order(merged_data$sex), ]
merged_data$sex[merged_data$sex == "2"] <- "0"
merged_data$sex <- as.numeric(merged_data$sex)

names(merged_data) <- c(names(merged_data)[1:17], "mean_HFEN+", "L5_HFEN+", "M5_HFEN+", "mean_NeishabouriCount_VM", "MVPA_T890_NeishabouriCount_y", names(merged_data)[23:24])

save(merged_data, file = paste0("/Users/annelindelettink/GECKO/preprocessing/manuscript/subset/analyses", "/merged_data.RData"))
