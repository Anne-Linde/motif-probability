### This script was used to create the data set corresponding to the following article:
## Unveiling hidden temporal physical behavior patterns: 
# A forward algorithm of the hidden semi-Markov model to capture motif probabilty beyond total volume and complexity
## Annelinde Lettink, Mai JM Chin A Paw, Eva Corpeleijn, Teatske M Altenburg, & Wessel N van Wieringen

# Load required libraries
library(haven)
library(dplyr)
library(TSEntropies)
devtools::install_github("Anne-Linde/motif-probability", force = TRUE)
library(motifprobability)
devtools::source_url("https://raw.githubusercontent.com/Anne-Linde/motif-probability/refs/heads/main/inst/researchcode/functions/lempel_ziv_complexity.R")

### User input required:
hsmmDir <- "/Users/annelindelettink/GECKO/preprocessing/manuscript/subset/output_rawinput_5years_complete_anthro/epochdata/validdata/mhsmmdata/models" # Contains the fitted HSMMs (obtained by running script: 1_preprocess_segment_accelerometer_data.R)

#Load anthropometric data and GGIR day summary
spss <- read_sav("/Users/annelindelettink/GECKO/preprocessing/manuscript/subset/complete_spss_GECKO_5years.sav") 
ggir <- read.csv("/Users/annelindelettink/GECKO/preprocessing/manuscript/subset/output_rawinput_5years_complete_anthro/results/rerun/part2_summary.csv")
ggir2 <- read.csv("/Users/annelindelettink/GECKO/preprocessing/manuscript/subset/output_rawinput_5years_complete_anthro/results/rerun/part5_personsummary_WW_L398M890V1254_T5A5.csv")

# Define acceleration ranges for PA and SB
# read in descriptive file with average accelerations in data 
desc <- read.csv("/Users/annelindelettink/Documents/Work MacBook Pro Annelinde/My Little Moves (MLM)/Sequence mapping/Physical behavior patterns/data/sequences/mhsmmdata/descriptives.csv")
M = mean(desc$Median)
Q3 = mean(desc$Q3)
threshold <- round((M+Q3)/2)
Amin.PA <- threshold
Amax.PA <- 9 # maximum value, note this value is a bit bigger to ensure all data was considered
Amin.SB <- 0 # minimum value
Amax.SB <- threshold
# Define epoch length
epoch_length <- 15 # hsmms were trained based on 15-sec epochs
n_epochs_min <- 60 / epoch_length

# Define motifs
motif.PA_30min <- defineMotif(Amin.PA, Amax.PA, 30, n_epochs_min) 
motif.PA_15min <- defineMotif(Amin.PA, Amax.PA, 15, n_epochs_min)
motif.PA_10min <- defineMotif(Amin.PA, Amax.PA, 10, n_epochs_min)
motif.PA_7min <- defineMotif(Amin.PA, Amax.PA, 7, n_epochs_min)
motif.PA_5min <- defineMotif(Amin.PA, Amax.PA, 5, n_epochs_min)
motif.PA_1min <- defineMotif(Amin.PA, Amax.PA, 1, n_epochs_min)
motif.SB_5min <- defineMotif(Amin.SB, Amax.SB, 5, n_epochs_min)
motif.SB_10min <- defineMotif(Amin.SB, Amax.SB, 10, n_epochs_min)
motif.SB_30min <- defineMotif(Amin.SB, Amax.SB, 30, n_epochs_min)
motif.SB_60min <- defineMotif(Amin.SB, Amax.SB, 60, n_epochs_min)
motifs <- c( "motif.PA_30min",  "motif.PA_15min", "motif.PA_10min", "motif.PA_7min", "motif.PA_5min", "motif.PA_1min", "motif.SB_5min", "motif.SB_10min", "motif.SB_30min", "motif.SB_60min")

######################
# Calculate the probabilties for all individuals: 
probabilities <- calculateMotifProbabilities(motifs, hsmmDir)
for(pp in 1:nrow(probabilities)){
  probabilities$id[pp] <- strsplit(probabilities$id[pp], "_")[[1]][1] # Extract participant ID
}
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
    #AD_L5_HFENplus_mg_0.24hr,   # average acceleration value during least active 5 hours
    AD_mean_NeishabouriCount_vm_mg_0.24hr, # average daily VM count
    AD_MVPA_E15S_T890_NeishabouriCount_y_0.24hr  # Time spent in moderate-to-vigorous based on 15 second epoch size and Sirard cut-point
  )
merged_data <- cbind(merged_data, selected_acc_data)

#match id's
GKNO <- c()
for(r in 1:nrow(ggir2)){
  GKNO <- c(GKNO, strsplit(strsplit(ggir2$ID, split = " ")[[r]][1], split = "_")[[1]][1])
}
ggir2$GKNO <- GKNO
matching_rows <- which(ggir2$GKNO %in% merged_data$GKNO)
ggir2 <- ggir2[matching_rows,]

selected_acc_data2 <- ggir2 %>%
  select(
    dur_day_total_IN_min_pla,
    dur_day_total_LIG_min_pla, 
    dur_day_total_MOD_min_pla,
    dur_day_total_VIG_min_pla,
    quantile_mostactive30min_mg_pla,
    quantile_mostactive60min_mg_pla,
    L5VALUE_pla,
    M5VALUE_pla,
    GKNO
  )

merged_data <- left_join(merged_data, selected_acc_data2, by = "GKNO")
merged_data$MVPA <- merged_data$dur_day_total_MOD_min_pla + merged_data$dur_day_total_VIG_min_pla

# Calculate complexity metrics
lz <- c()
sampleEntropy <- c()

filelist <- list.files(hsmmDir)
for(pp in 1:length(filelist)){
  load(paste0(hsmmDir, "/", filelist[pp]))# Load the sequences
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

#names(merged_data) <- c(names(merged_data)[1:17], "mean_HFEN+", "L5_HFEN+", "M5_HFEN+", "mean_NeishabouriCount_VM", "MVPA_T890_NeishabouriCount_y", names(merged_data)[23:24])
names(merged_data) <- c(names(merged_data)[1:17], "avg_acc_mg", "L5_all", "avg_acc_VMcts", "MVPA_all_VMcts", "SB", "LPA", "MPA", "VPA", "M30", "M60", "L5", "M5", names(merged_data)[30:32])

#save(merged_data, file = paste0("/Users/annelindelettink/GECKO/preprocessing/manuscript/subset/analyses", "/merged_data.RData"))
save(merged_data, file = paste0("/Users/annelindelettink/GECKO/preprocessing/manuscript/subset/analyses", "/merged_data_rerun.RData"))
