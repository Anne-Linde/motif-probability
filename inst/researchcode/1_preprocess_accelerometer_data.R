# Tabula rasa
rm(list=ls())
gc()

### STEP 1: Subset the GECKO data (for accelerometer data and anthropometry at 5 years)
### STEP 2: Preprocess the accelerometer data, by:
## a) loading the data, aggregating 15-sec epochs HFEN+, calculating time-use estimates Sirard cutpoints, using GGIR
## b) epoch level data and save valid data

# Libraries
library(haven)
library(GGIR)
library(actilifecounts)
library(anthro)

### STEP 1) Subset the GECKO data

## User settings
rawaccdir <- "//vumc.nl/afd$/DIV10/POH/Sectie_2/2016_SB and PA pattern analysis/Team Annelinde/Data/GECKO Drenthe/Accelerometer data/rawinput/gt3x"
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

length(unique(numberPP)) #1,518 participants with acceleration files

# Informed consent is missing for n = 1 "4292-1_20juli11 (2011-08-05).gt3x, check that this file is excluded after subsetting"

spssDataAcc <- spssdata[which(spssdata$GKNO %in% GECKOid),] # select the data
rm(GECKOid)

idFive <- accspss$GKNO[accspss$Age == 5] # only children under 5 years
spss_merged <- merge(spssDataAcc, accspss)
spss_merged_5years <- spss_merged[which(spss_merged$Age == 5),] # n = 619

numberPP <- c()
for (file in 1:nrow(spss_merged_5years)) {
  numberPP <- c(numberPP, strsplit(spss_merged_5years$Filename[file], "-")[[1]][1])
}
length(unique(numberPP)) #611 participants with acceleration files

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
length(multiple_days) # n=7 with multiple measurements
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

### STEP 2: PREPROCESS ACCELEROMETER DATA
#a) Load raw data and calculate metrics for 15-sec epochs separate for each age group, as different cut-points are needed for MVPA

GGIR::GGIR(
  datadir = paste0(storedir, "/rawinput"),
  outputdir = "\\vumc.nl\afd$\DIV10\POH\Sectie_2\2016_SB and PA pattern analysis\Team Annelinde\Data\GECKO Drenthe\Accelerometer data\Preprocessed\output_5years60cts0", # Save the raw data and calculated metrics here
  mode = c(1, 2, 3, 4, 5),
  #do.report = c(2, 4, 5),
  do.parallel = FALSE,
  
  ## Part 1 –  data processing
  windowsizes = c(15, 900, 3600), # Epoch length: 15 sec, non-wear time 60min*60 = 3600
  minimumFileSizeMB = 0.1,
  
  # Metrics
  do.neishabouricounts = TRUE,
  do.hfenplus     = TRUE,
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

# b) Process accelerometer files to epoch level and store valid data only
# Note: I copied the preprocessed accelerometer data to: the newdatafolder
# User input
newdatafolder = "/Users/annelindelettink/GECKO/preprocessing/manuscript/output_60min0values/meta/ms2.out" # Save the raw data and calculated metrics here
epochdir = "/Users/annelindelettink/GECKO/preprocessing/manuscript/output_60min0values/epochdata"
validdatadir = paste0(epochdir, "/validdata")
files <- list.files(path = newdatafolder, pattern = ".RData")

hoursValidDay = 10
minDays = 3
epoch = 15
epochMin <- epoch/60

calibration_error_na <- 0
calibration_error_bigger <- 0
invaliddata <- 0

for(file in 1:length(files)){
  
  if(file.exists(paste0(epochdir, files[file])) && overwrit == FALSE){
    print("epochdata already saved")
    load(paste(epochdir, files[file], sep = "/")) # Load epochdata
  } else {
    load(paste(newdatafolder, files[file], sep = "/")) # Load RData file
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
      save(validData, file = paste0(validdatadir, "/", filename))
    }else {
      invaliddata <- invaliddata + 1
    }
  }
}
