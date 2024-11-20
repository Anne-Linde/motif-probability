### STEP 3: Train hidden-semi Markov models to segment the preprocessed accelerometer data into movement sequence maps, by 
## a) Generating the mhsmm data format for the HFENplus variable
## b) Fitting hidden semi-Markov models for each file

# Libraries
library(mhsmm)

rDir <- "/Users/annelindelettink/Documents/Work MacBook Pro Annelinde/My Little Moves (MLM)/Sequence mapping/Physical behavior patterns/motif-probability/R/"
setwd(rDir)
source("format_mhsmm.R")
source("derive_initialparams.R")

## User settings
epochdir <- "/Users/annelindelettink/GECKO/preprocessing/manuscript/output_60min0values/epochdata"
#epochdir <- "/Users/annelindelettink/Documents/Work MacBook Pro Annelinde/My Little Moves (MLM)/Accelerometer data/Measurement1/Motif probability/epochdata"
#epochdir <- paste0("C://Users//P070400//OneDrive - Amsterdam UMC//Documenten/", "epochdata/")
validdatadir <- paste0(epochdir, "/validdata/")
#sequencedir <- "C://Users//P070400//OneDrive - Amsterdam UMC//Documenten/hsmms/"
sequencedir <- paste0(validdatadir, "mhsmmdata/sequences/")

nStates <- 3 #For now 3 (if more states need to be trained look at: ~/Documents/Werk/PROGRAMMING/script/hsmm_functions.R)
qMax <- 15 # Maximum number of iterations
dMax <- (60*60) / epoch  # Maximum state duration of 240 15-sec epochs, corresponds to 60 minutes - as was set in van Kuppevelt et al. (2019)
boutDurations <- c(30, 10, 5) # Note that the duration of the first bout has to be non-zero (as lambda is required to be > 0)
variable_name <- "HFENplus" #metric for training hsmm


## First, check if there are still double files/multiple measurements
validfiles <- list.files(validdatadir, pattern = ".RData")
f_name <- c()
for (file in 1:length(validfiles)) {
  f_name <- c(f_name, strsplit(validfiles[file], "_")[[1]][1])
}

validfiles[which(duplicated(f_name))] # No double files

# Create a vector of the files that are associated with multiple measurements
f_name[which(grepl("-2", f_name))]

load(paste0(validdatadir, "1155-1_19mei11.RData"))
load(paste0(validdatadir, "1155-2_19mei11.RData")) # Both equal number of days data, but more data for last so keep the last
validfiles <- validfiles[-which(validfiles == "1155-1_19mei11.RData")]

load(paste0(validdatadir, "1770-1_27feb12.RData"))
load(paste0(validdatadir, "1770-2_27feb12.RData")) # Both equal number of days data, but more data for last so keep the last
validfiles <- validfiles[-which(validfiles == "1770-1_27feb12.RData")]

load(paste0(validdatadir, "1847-1_26okt11.RData"))
load(paste0(validdatadir, "1847-2_25okt11.RData")) # Both equal number of days data, but more data for last so keep the last
validfiles <- validfiles[-which(validfiles == "1847-1_26okt11.RData")]

load(paste0(validdatadir, "2774-1_19APR12.RData"))
load(paste0(validdatadir, "2774-2_19APR12.RData")) # Both equal number of days data, so keep the first
validfiles <- validfiles[-which(validfiles == "2774-2_19APR12.RData")]

#load(paste0(validdatadir, "2934-1_15juni11.RData"))
load(paste0(validdatadir, "2934-2_15juni11.RData")) # First measurement is not valid
#validfiles <- validfiles[-which(validfiles == "2934-2_15juni11.RData")]

load(paste0(validdatadir, "3497-1_7dec11.RData"))
load(paste0(validdatadir, "3497-2_7dec11.RData")) # Both equal number of days data, but more data for last so keep the last
validfiles <- validfiles[-which(validfiles == "3497-1_7dec11.RData")]

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
  mhsmmdata <- format_mhsmm(data, variable_name)
  mhsmmdata$Y <- asinh(mhsmmdata$Y*1000) # Scale and normalize the observations 
  mhsmmdata$Y[which(mhsmmdata$Y == 0)] <- rnorm(n = sum(mhsmmdata$Y == 0), mean = 0, sd = 0.05) # Add a little noise to zeros
  save(file = paste0(paste0(validdatadir, "mhsmmdata/"), filename), mhsmmdata) 
}

## b) Fit Hsmms per individual
epoch = 15
variable_name <- "HFENplus"
for(file in 1:length(validfiles)){
  # Load the data from the file

  if(!file.exists(paste0(sequencedir, validfiles[file]))){
    load(paste0(paste0(validdatadir, "mhsmmdata/"), validfiles[file]))
    
    # Fit the hsmm
    initParams <- derive_initialparams(mhsmmdata$Y, nStates, boutDurations, epoch)
    startmodel <- hsmmspec(initParams$initial, initParams$transition, initParams$emission, 
                           initParams$sojourn_distribution, dens.emis = dnorm.hsmm)
    hsmms <- hsmmfit(mhsmmdata$Y, startmodel, mstep = mstep.norm, maxit = qMax, M = dMax) 
    save(file = paste0(sequencedir, validfiles[file]), hsmms)
  }

}
