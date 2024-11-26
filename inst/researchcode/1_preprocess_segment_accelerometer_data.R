# Tabula rasa
rm(list=ls())
gc()

## This script:
### STEP 1: Subsets the GECKO data for accelerometer measurement at 5 years and complete data on BMI

### STEP 2: Preprocesses the accelerometer data, by
## a) Load raw data and calculate metrics for 15-sec epochs
## b) Process accelerometer files to epoch level and store valid data

### STEP 3: Trains hidden-semi Markov models to segment the accelerometer data, by
## a) Generating the mhsmm data format for the HFENplus variable
## b) Fitting hidden semi-Markov models for each file


### STEP 1) Subset the GECKO data:
library(haven)

## User settings
rawaccdir <- "//vumc.nl/afd$/DIV10/POH/Sectie_2/2016_SB and PA pattern analysis/Team Annelinde/Physcial behavior patterns/GECKO data/rawinput/gt3x"
rawspssdir <- "//vumc.nl/afd$/DIV10/POH/Sectie_2/2016_SB and PA pattern analysis/Team Annelinde/Data/GECKO Drenthe"
spssdata <- haven::read_sav(paste0(rawspssdir, "/Deel 1/Data_GECKO_MyLittleMoves_Deel1_2020-10-22.sav")) # Antropometrics
accspss <- haven::read_sav(paste0(rawspssdir, "/Data opschonen/GECKONRS_ALL_Butte2014VM_RW_opgeschoond.sav"))
storedir <- "//vumc.nl/afd$/DIV10/POH/Sectie_2/2016_SB and PA pattern analysis/Team Annelinde/Physcial behavior patterns/GECKO data/subset"

# Check n participants with accelerometer data and in spss file with antropometrics
accfiles <- list.files(rawaccdir)
GECKOid <- c() # Abbreviate participant id from accelerometer data to correspond with SPSS data
numberPP <- c()
for (file in 1:length(accfiles)) {
  GECKOid <- c(GECKOid, strsplit(accfiles[file], "_")[[1]][1])
  numberPP <- c(numberPP, strsplit(accfiles[file], "-")[[1]][1])
}
length(unique(numberPP)) #1518 participants with acceleration files

spssDataAcc <- spssdata[which(spssdata$GKNO %in% GECKOid),] # select the data
rm(GECKOid)
idFive <- accspss$GKNO[accspss$Age == 5] # only children under 5 years
spss_merged <- merge(spssDataAcc, accspss)

which(spssdata$GKNO == "4292-1")

spss_merged_5years <- spss_merged[which(spss_merged$Age == 5),] # n = 619

numberPP <- c()
for (file in 1:nrow(spss_merged_5years)) {
  numberPP <- c(numberPP, strsplit(spss_merged_5years$Filename[file], "-")[[1]][1])
}
length(unique(numberPP)) #611 participants with acceleration files
#1518-611 = 907, 907-142 = 765 pp age != 5 years

# Check which participants have antropometric data
Rsex <- which(is.na(spss_merged_5years$sex)) # all complete
Rage <- which(is.na(spss_merged_5years$measure_age_217)) # 43 missings
Rweight <- which(is.na(spss_merged_5years$GEWICHT_217)) # 51 missings
Rlenght <- which(is.na(spss_merged_5years$LENGTE_217)) # 51 missings
Rbmi <- unique(c(Rsex,Rage, Rweight, Rlenght))

idMissings <- spss_merged_5years$GKNO[Rbmi]
numberPP <- numberPP[-Rbmi]
length(unique(numberPP)) #561 pp with complete data

spss_complete <- spss_merged_5years[-which(spss_merged_5years$GKNO %in% idMissings),]
spss_complete$Filename <- gsub("15sec\\.agd", ".gt3x",spss_complete$Filename)

write_sav(spss_complete, paste0(storedir, "/complete_spss_GECKO_5years.sav")) #save spss file
rm(accspss, spssdata, spssDataAcc, spssDataAccAge, spss_merged, spss_merged_5years)

# Check which participants have multiple files (either double or multiple measurements)
accfiles_complete <- spss_complete$Filename
multiple_days <- c() # measured multiple days; these data can be included
double_file <- c() # measured twice or more for 1 measurement point

start_pp <- spss_complete$GKNO[1]
measurement <- strsplit(strsplit(accfiles_complete[1], "_")[[1]][1], "-")[[1]][2]
for(subj in 2:(length(accfiles_complete))) {
  next_pp <- strsplit(strsplit(accfiles_complete[subj], "_")[[1]][1], "-")[[1]][1]
  next_measurement <- strsplit(strsplit(accfiles_complete[subj], "_")[[1]][1], "-")[[1]][2]
  if(start_pp == next_pp) {
    if(measurement == next_measurement) {
      double_file <- c(double_file, next_pp)
    } else {
      multiple_days <- c(multiple_days, next_pp)
    }
  } else {
    start_pp <- next_pp
  }
}
rm(next_measurement, measurement, next_pp, start_pp, subj)

# Check other reasons for exclusion; following syntax/reasoning from GECKO cohort
remove <- c()
remove <- c(remove, which(accfiles_complete == "1168-1_31-5-2013 (2013-06-27).gt3x")) # downloaded twice
remove <- c(remove, which(accfiles_complete == "1735-1_11APR12 (2012-05-29).gt3x")) # something went wrong with data, measured again; remove first try
remove <- c(remove, which(accfiles_complete == "3074-1_15sept11 (2011-10-14).gt3x")) # data parents and actilife do not correspond, therefore other day was measured
remove <- c(remove, which(accfiles_complete == "3473-1_28juni11 (2011-08-01).gt3x")) # went in washing machine, new accelerometer was sent
remove <- c(remove, which(accfiles_complete == "3960-1_06mei11 (2011-05-24).gt3x")) # data parents and actilife do not correspond and accelerometer probably worn by someone else, therefore other day was measured
# remove is empty

# Copy raw gt3x files to be preprocessed
for (file in 1:length(spss_complete$Filename)) {
  if(!file.exists(paste0(storedir, "/rawinput/", spss_complete$Filename[file]))){
    file.copy(from = paste0(rawaccdir, "/", spss_complete$Filename[file]),
              to = paste0(storedir, "/rawinput/", spss_complete$Filename[file]), overwrite = FALSE)
  }
}

# Check if there are still double files
completefiles <- spss_complete$Filename
f_name <- c()
for (file in 1:length(completefiles)) {
  f_name <- c(f_name, strsplit(completefiles[file], "_")[[1]][1])
}
completefiles[which(duplicated(f_name))] # No double files

### STEP 2: PREPROCESS ACCELEROMETER DATA
library(GGIR)

#a) Load raw data and calculate metrics for 15-sec epochs separate for each age group, as different cut-points are needed for MVPA
GGIR::GGIR(
  datadir = "//vumc.nl/afd$/DIV10/POH/Sectie_2/2016_SB and PA pattern analysis/Team Annelinde/Physcial behavior patterns/GECKO data/subset/rawinput_5years_complete_anthro",
  outputdir = "//vumc.nl/afd$/DIV10/POH/Sectie_2/2016_SB and PA pattern analysis/Team Annelinde/Physical Behavior Patterns/GECKO data/subset/output_rawinput_5years_complete_anthro", # Save the raw data and calculated metrics here
  mode = c(1, 2, 3, 4, 5),
  #do.report = c(2, 4, 5),
  do.parallel = FALSE,
  ## Part 1 – data processing
  windowsizes = c(15, 900, 5400), # Epoch length: 15 sec, non-wear time 90min*60 = 5400
  # Metrics
  do.neishabouricounts = TRUE,
  do.hfenplus = TRUE,
  acc.metric = "NeishabouriCount_y",
  ## Part 2 – data quality and descriptive
  qlevels <- c(c(1380/1440), c(1410/1440), c(1430/1140)),
  includedaycrit = 10,
  mvpathreshold = 890,
  overwrite = FALSE, #do not overwrite data
  
  #part 5 - total volume SB, LPA MVPA use Sirard cut-points
  ignorenonwear = TRUE,
  threshold.lig = c(398),
  threshold.mod = c(890),
  threshold.vig = c(1254)
)

## c) Process accelerometer files to epoch level and store valid data only
# Note: I copied the preprocessed accelerometer data to: the newdatafolder
newdatafolder <- "//vumc.nl/afd$/DIV10/POH/Sectie_2/2016_SB and PA pattern analysis/Team Annelinde/Physcial behavior patterns/GECKO data/subset/output_rawinput_5years_complete_anthro"
files <- list.files(path = paste0(newdatafolder, "/meta/ms2.out"))
epochdir <- paste0(newdatafolder, "/epochdata")
validdatadir <- paste0(epochdir, "/validdata/")

calibration_error_na <- 0
calibration_error_bigger <- 0

for(file in 1:length(files)){
  
  if((file.exists(paste0(epochdir, "/", files[file])) & overwrit == FALSE)){
    print("epochdata already saved")
    load(paste(epochdir, files[file], sep = "/")) # Load epochdata
  } else {
    load(paste(newdatafolder, files[file], sep = "/meta/ms2.out/")) # Load RData file
    if(!is.na(as.double(SUM$summary$calib_err)) & as.double(SUM$summary$calib_err) < 0.01) { #Only load file if the post-calibration error is < 0.01 g
      # Convert timestamp to POSIX for convenience
      tsFormat = "%Y-%m-%dT%H:%M:%S%z"
      tz = ""
      IMP$metashort$timestampPOSIX = as.POSIXlt(IMP$metashort$timestamp, format = tsFormat, tz = tz)
      
      # Interpolate scores as they are in a different resolution
      shortEpochLength = IMP$windowsizes[1]
      longEpochLength = IMP$windowsizes[2]
      NlongInsideShort = longEpochLength / shortEpochLength
      scores = IMP$rout[rep(1:nrow(IMP$rout), each = NlongInsideShort),]
      colnames(scores) = c("nonwear", "clipping", "additonal_nonwear", "studyprotocol", "all") # add column names
      IMP$metashort = cbind(IMP$metashort, scores) # combine scores with the epoch level time series
      epochdata <- list(agg.epoch = IMP$metashort, day.metrics = SUM$daysummary)
      save(epochdata, file = paste0(epochdir, "/", files[file]))
    }
    else{
      if(is.na(as.double(SUM$summary$calib_err)))
        calibration_error_na <- calibration_error_na + 1
      if(as.double(SUM$summary$calib_err) >= 0.01)
        calibration_error_bigger <- calibration_error_bigger + 1
    }
  }
  # Save valid day data (at least 3 days of at least 10 hours)
  data_perday <- split(epochdata$agg.epoch, as.Date(epochdata$agg.epoch$timestampPOSIX))
  validDays <- c()
  for(day in 1:length(data_perday)){
    dataDay <- data_perday[[day]]
    
    epochsDay <- nrow(dataDay) # Number of epochs in this day
    nwEpochs <- length(c(which(dataDay$nonwear == 1), which(dataDay$additonal_nonwear == 1))) # Number of epochs labeled as non-wear in this day
    nValidEpochs <- epochsDay - nwEpochs
    validDays <- c(validDays, (nValidEpochs > (60 * hoursValidDay * epochMin)))
  }
  validData <- data_perday[validDays]
  if(length(validData) > 0){
    for(day in 1:length(validData)){
      dataDay <- validData[[day]]
      #remove non-wear data
      nonwearIndex <- c(which(dataDay$nonwear == 1), which(dataDay$additonal_nonwear == 1))
      if(length(nonwearIndex) > 0){
        validData[[day]] <- dataDay[-nonwearIndex,]
      }
    }
    if(length(validData) >= minDays){ # save data only if there are at least # valid days
      filename = paste0(strsplit(files[file], " ")[[1]][1], ".RData")
      save(validData, file = paste0(validdatadir, filename))
    }
  }
}

# Check if there are still double files
validfiles <- list.files(validdatadir, pattern = ".RData")
f_name <- c()
for (file in 1:length(validfiles)) {
  f_name <- c(f_name, strsplit(validfiles[file], "_")[[1]][1])
}

validfiles[which(duplicated(f_name))] # No double files

# Create a vector of the files that are associated with multiple measurements
f_name[which(grepl("-2", f_name))]

load(paste0(validdatadir, "1155-1_19mei11.RData")) # 5 days
load(paste0(validdatadir, "1155-2_19mei11.RData")) # 4 days
validfiles <- validfiles[-which(validfiles == "1155-2_19mei11.RData")]

load(paste0(validdatadir, "1770-1_27feb12.RData"))
load(paste0(validdatadir, "1770-2_27feb12.RData")) # Both equal number of days data, more data for first measurement
validfiles <- validfiles[-which(validfiles == "1770-2_27feb12.RData")]

load(paste0(validdatadir, "1847-1_26okt11.RData"))
load(paste0(validdatadir, "1847-2_25okt11.RData")) # Both equal number of days data, more data for first measuremenet
validfiles <- validfiles[-which(validfiles == "1847-2_25okt11.RData")]

load(paste0(validdatadir, "2774-1_19APR12.RData"))
load(paste0(validdatadir, "2774-2_19APR12.RData")) # Both equal number of days data, more data for first measuremnet
validfiles <- validfiles[-which(validfiles == "2774-2_19APR12.RData")]

load(paste0(validdatadir, "2934-1_15juni11.RData"))
load(paste0(validdatadir, "2934-2_15juni11.RData")) # Both equal number of days data, more data for the last
validfiles <- validfiles[-which(validfiles == "2934-1_15juni11.RData")]

load(paste0(validdatadir, "3497-1_7dec11.RData"))
load(paste0(validdatadir, "3497-2_7dec11.RData")) # Both equal number of days data, so keep the first
validfiles <- validfiles[-which(validfiles == "3497-2_7dec11.RData")]

load(paste0(validdatadir, "3589-1_3APR12.RData"))
load(paste0(validdatadir, "3589-2_3APR12.RData")) # Both equal number of days data, so keep the first
validfiles <- validfiles[-which(validfiles == "3589-2_3APR12.RData")]

### STEP 2: Segments the accelerometer data by fitting hidden semi-Markov models for each file
## a) Generate mhsmm data for validfiles
for(f in 1:length(validfiles)){
  cat(f)
  load(paste0(validdatadir, validfiles[f]))
  data <- validData
  filename = validfiles[f]
  
  #format data for hsmm
  N <- c()
  df <- data.frame()
  # Vectors including the data observations:
  for (day in 1:length(data)) {
    col = which(colnames(data[[day]]) == variable_name)
    df <- rbind(df, data[[day]][col])
    N <- c(N, nrow(data[[day]][col]))
  }
  colnames(df) <- c("Y")
  mhsmmdata <- list(Y = df[["Y"]], N = N)
  class(mhsmmdata) <- "hsmm.data"
  mhsmmdata$Y <- asinh(mhsmmdata$Y*1000) # Scale the observations
  mhsmmdata$Y[which(mhsmmdata$Y == 0)] <- rnorm(n = sum(mhsmmdata$Y == 0), mean = 0, sd = 0.05) # Add a little noise to zeros
  save(file = paste0(paste0(validdatadir, "mhsmmdata/"), filename), mhsmmdata)
}

## b) Fit hsmms per individual
#install.packages("DescTools")
library(mhsmm)

sequencedir <- paste0(validdatadir, "mhsmmdata/models/")
epoch <- 15
nStates <- 3 #For now 3 (if more states need to be trained look at: ~/Documents/Werk/PROGRAMMING/script/hsmm_functions.R)
qMax <- 15 # Maximum number of iterations
dMax <- (60*60) / epoch # Maximum state duration of 240 15-sec epochs, corresponds to 60 minutes - as was set in van Kuppevelt et al. (2019)
boutDurations <- c(30, 10, 5) # Note that the duration of the first bout has to be non-zero (as lambda is required to be > 0)
variable_name <- "HFENplus" #metric for training hsmm

# Load functions required for hsmm
r_files <- list.files(path = "//vumc.nl/afd$/DIV10/POH/Sectie_2/2016_SB and PA pattern analysis/Team Annelinde/Physcial behavior patterns/Physical behavior patterns/R"
                      , pattern = [file://.R$]\\.R$, full.names = TRUE)
lapply(r_files, source)

files <- list.files(path = paste0(validdatadir, "mhsmmdata/"), pattern = ".RData")

for(file in 1:length(files)){
  # Load the data from the file
  if(!file.exists(paste0(sequencedir, files[file]))){
    load(paste0(paste0(validdatadir, "mhsmmdata/"), validfiles[file]))
    # Fit the hsmm
    initParams <- derive_initialparams(mhsmmdata$Y, nStates, boutDurations, epoch)
    startmodel <- hsmmspec(
      initParams$initial,
      initParams$transition,
      initParams$emission,
      initParams$sojourn_distribution,
      dens.emis = dnorm.hsmm
    )
    hsmms <- hsmmfit(mhsmmdata$Y, startmodel, mstep = mstep.norm, maxit = qMax, M = dMax)
    # Save the best model as an .RData file
    save(hsmms, file = paste0(sequencedir, files[file]))
  }
}
