### STEP 2: Trains hidden-semi Markov models to segment the preprocessed accelerometer data, by 
## a) Generating the mhsmm data format for the HFENplus variable
## b) Fitting hidden semi-Markov models for each file

# Libraries
library(mhsmm)

rDir <- "/Users/annelindelettink/Documents/Work MacBook Pro Annelinde/My Little Moves (MLM)/Sequence mapping/Physical behavior patterns/motif-probability/R/"
setwd(rDir)
source("format_mhsmm.R")
source("cross_validate_hsmm.R")


## User settings
epochdir <- paste0("C://Users//P070400//OneDrive - Amsterdam UMC//Documenten/", "epochdata/")
validdatadir <- paste0(epochdir, "validdata/")
sequencedir <- "C://Users//P070400//OneDrive - Amsterdam UMC//Documenten/hsmms/"
nStates <- 3 #For now 3 (if more states need to be trained look at: ~/Documents/Werk/PROGRAMMING/script/hsmm_functions.R)
qMax <- 15 # Maximum number of iterations
dMax <- (60*60) / epoch  # Maximum state duration of 240 15-sec epochs, corresponds to 60 minutes - as was set in van Kuppevelt et al. (2019)
boutDurations <- c(30, 10, 5) # Note that the duration of the first bout has to be non-zero (as lambda is required to be > 0)
variable_name <- "HFENplus" #metric for training hsmm
# Load functions required for hsmm
r_files <- list.files(path = "//vumc.nl/afd$/DIV10/POH/Sectie_2/2016_SB and PA pattern analysis/Team Annelinde/Physcial behavior patterns/Physical behavior patterns/R"
                      , pattern = "\\.R$", full.names = TRUE)
lapply(r_files, source)

## First, check if there are still double files
validfiles <- list.files(validdatadir, pattern = ".RData")
f_name <- c()
for (file in 1:length(validfiles)) {
  f_name <- c(f_name, strsplit(validfiles[file], "_")[[1]][1])
}

validfiles[which(duplicated(f_name))] # No double files

# Create a vector of the files that are associated with multiple measurements
f_name[which(grepl("-2", f_name))]

load(paste0(validdatadir, "1155-1_19mei11.RData"))
load(paste0(validdatadir, "1155-2_19mei11.RData")) # Both equal number of days data, so keep the first
validfiles <- validfiles[-which(validfiles == "1155-2_19mei11.RData")]

load(paste0(validdatadir, "1770-1_27feb12.RData"))
load(paste0(validdatadir, "1770-2_27feb12.RData")) # Both equal number of days data, so keep the first
validfiles <- validfiles[-which(validfiles == "1770-2_27feb12.RData")]

load(paste0(validdatadir, "1847-1_26okt11.RData"))
load(paste0(validdatadir, "1847-2_25okt11.RData")) # Both equal number of days data, so keep the first
validfiles <- validfiles[-which(validfiles == "1847-2_25okt11.RData")]

load(paste0(validdatadir, "2774-1_19APR12.RData"))
load(paste0(validdatadir, "2774-2_19APR12.RData")) # Both equal number of days data, so keep the first
validfiles <- validfiles[-which(validfiles == "2774-2_19APR12.RData")]

load(paste0(validdatadir, "2934-1_15juni11.RData"))
load(paste0(validdatadir, "2934-2_15juni11.RData")) # Both equal number of days data, so keep the first
validfiles <- validfiles[-which(validfiles == "2934-2_15juni11.RData")]


load(paste0(validdatadir, "3497-1_7dec11.RData"))
load(paste0(validdatadir, "3497-2_7dec11.RData")) # Both equal number of days data, so keep the first
validfiles <- validfiles[-which(validfiles == "3497-2_7dec11.RData")]

load(paste0(validdatadir, "3589-1_3APR12.RData"))
load(paste0(validdatadir, "3589-2_3APR12.RData")) # Both equal number of days data, so keep the first
validfiles <- validfiles[-which(validfiles == "3589-2_3APR12.RData")]



## a) Generate mhsmm data for validfiles
for(f in 1:length(validfiles)){
  cat(f)
  load(paste0(validdatadir, validfiles[f]))
  data <- validData
  filename = validfiles[f]
  
  #format data for hsmm
  mhsmmdata <- format_mhsmm(data)
  mhsmmdata$Y <- asinh(mhsmmdata$Y*1000) # Scale the observations 
  mhsmmdata$Y[which(mhsmmdata$Y == 0)] <- rnorm(n = sum(mhsmmdata$Y == 0), mean = 0, sd = 0.05) # Add a little noise to zeros
  save(file = paste0(paste0(validdatadir, "mhsmmdata2/"), filename), mhsmmdata) 
}

## b) Fit Hsmms, but using k-fold cross-validation per individual to avoid overfitting
files <- list.files(path = paste0(validdatadir, "mhsmmdata/"), pattern = ".RData")

for(file in 1:length(files)){
  # Load the data from the file
  load(paste0(paste0(validdatadir, "mhsmmdata/"), files[file]))
  #Perform cross-validation
  results <- cross_validate_hsmm(mhsmmdata, k_folds = 5, nStates, boutDurations, epoch, qMax, dMax)
  
  # Save the best model as an .RData file
  best_model_filename <- paste0(sequencedir, file.path("trained_models", files[file]))
  best_model <- results$best_model
  save(best_model, file = best_model_filename)
  sequence_data <- results$sequence_data
  sequence_filename <- paste0(sequencedir, file.path("segmented_data", files[file]))
  save(sequence_data, file = sequence_filename)
}
