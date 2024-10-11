# Tabula rasa
rm(list=ls())
gc()

### STEP 1: Preprocesses the accelerometer data, by 
## a) Selecting the accelerometer data files for which data is complete (child is >=3 -years-old <=5 and data on BMI is available (copy this data to new folder: completedata), and saves the anthropometrics
## b) Load raw data and calculate metrics for 15-sec epochs separately for each age group, as different cut-points are needed for MVPA
## c) Process accelerometer files to epoch level and store valid data only

# Libraries
library(haven)
library(GGIR)
library(actilifecounts)
library(anthro)

## User settings
dir <- "M://Sectie_2/2016_SB and PA pattern analysis/Team Annelinde/Data/GECKO Drenthe/"
accdir <- paste0(dir, "Accelerometer data/rawinput/gt3x")
storedir <- paste0(dir, "Accelerometer data/Preprocessed")
spssanthrodir <- paste0(dir, "Deel 1")
spssaccdir <- paste0(dir, "Data opschonen")
spssdata <- haven::read_sav(paste0(spssanthrodir, "/Data_GECKO_MyLittleMoves_Deel1_2020-10-22.sav")) # Antropometrics
accspss <- haven::read_sav(paste0(spssaccdir, "/GECKONRS_ALL_Butte2014VM_RW_opgeschoond.sav"))

filedir = paste0(dir, "Accelerometer data/rawinput/csv/completedata/5years") #5-year-olds
#filedir = paste0(dir, "Accelerometer data/rawinput/csv/completedata/4years") #4-year-olds
#filedir = paste0(dir, "Accelerometer data/rawinput/csv/completedata/3years") #3-year-olds
cutpoint = 890 #MVPA Sirard cut-off point 890 for 5-year-olds, 811 for 4-year-olds, 614 for 3-year-olds

newdatafolder = "M://Sectie_2/2016_SB and PA pattern analysis/Team Annelinde/Data/GECKO Drenthe/Accelerometer data/Preprocessed/output_5years90cts0/meta/ms2.out/"
epochdir <- paste0("C://Users//P070400//OneDrive - Amsterdam UMC//Documenten/", "epochdata/")
validdatadir <- paste0(epochdir, "validdata/")
overwrit = FALSE
epoch = 15 # epoch length (sec)
epochMin = 60/epoch # number of epochs in one minute
hoursValidDay <- 10 # at least 10 hours for a valid day
minDays <- 3 # at least 3 days


## a) Select complete data only (>3 years <=5, complete data on height and weight status, gender, accelerometer file)
accfiles <- list.files(accdir)
GECKOid <- c() # Abbreviate participant id from accelerometer data to correspond with SPSS data
uniquepp <- c()
for (file in 1:length(accfiles)) {
  GECKOid <- c(GECKOid, strsplit(accfiles[file], "_")[[1]][1])
  uniquepp <- c(uniquepp, strsplit(GECKOid[file], "-")[[1]][1])
}

# Check which participants have multiple files (either double or multiple measurements)
multiple_days <- c() # measured multiple days; these data can be included
double_file <- c() # measured twice or more for 1 measurement point

start_pp <- GECKOid[1]
measurement <- strsplit(strsplit(accfiles[1], "_")[[1]][1], "-")[[1]][2]
for(subj in 2:(length(accfiles))) {
  next_pp <- strsplit(strsplit(accfiles[subj], "_")[[1]][1], "-")[[1]][1]
  next_measurement <- strsplit(strsplit(accfiles[subj], "_")[[1]][1], "-")[[1]][2]
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

# Reasons for exclusion; following syntax/reasoning from GECKO cohort
remove <- c()
remove <- c(remove, which(accfiles == "1168-1_31-5-2013 (2013-06-27).gt3x")) # downloaded twice
remove <- c(remove, which(accfiles == "1735-1_11APR12 (2012-05-29).gt3x")) # something went wrong with data, measured again; remove first try
remove <- c(remove, which(accfiles == "3074-1_15sept11 (2011-10-14).gt3x")) # data parents and actilife do not correspond, therefore other day was measured
remove <- c(remove, which(accfiles == "3473-1_28juni11 (2011-08-01).gt3x")) # went in washing machine, new accelerometer was sent
remove <- c(remove, which(accfiles == "3960-1_06mei11 (2011-05-24).gt3x")) # data parents and actilife do not correspond and accelerometer probably worn by someone else, therefore other day was measured
remove <- c(remove, which(accfiles == "4292-1_20juli11 (2011-08-05).gt3x")) # informed consent missing
if(length(remove)>0){
  accfiles <- accfiles[-remove] # Remove the (dubious) files from the file list to ensure these are not loaded
}
files <- c("1168", "1735", "3074", "3473", "3960", "4292")
double_file_corrected <- double_file[!double_file %in% files]
rm(files, double_file)

GECKOid <- c()
for (file in 1:length(accfiles)) {
  GECKOid <- c(GECKOid, strsplit(accfiles[file], "_")[[1]][1])
}

spssDataAcc <- spssdata[which(spssdata$GKNO %in% GECKOid),] # select the data
rm(spssdata)

# Check which participants age < 5 years
merged_df <- merge(accspss, spssDataAcc, by = "GKNO")
idFiveYounger <- merged_df$GKNO[merged_df$Age <= 5] # only children under 5 years

#Check which participants have data on height and weight, and save these in new data.frame
df_new <- data.frame(idFiveYounger)
Rweight <- c()
Rheight <- c()
weight <- c()
height <- c()
younger <- c()
for(id in 1:length(idFiveYounger)){
  df_new[id, 2] <- merged_df$Filename[merged_df$GKNO == idFiveYounger[id]]
  if(merged_df$Age[merged_df$GKNO == idFiveYounger[id]] == "5"){ # Lengte en gewicht op 5 jaar gemeten
    index_weight <- which(colnames(merged_df) == "GEWICHT_217")
    index_height <- which(colnames(merged_df) == "LENGTE_217")
    df_new[id, 3] <- 5
  } else if(merged_df$Age[merged_df$GKNO == idFiveYounger[id]] == "4"){ # Lengte en gewicht op 3 jaar + 9 maanden gemeten
    index_weight <- which(colnames(merged_df) == "GEWICHT_216")
    index_height <- which(colnames(merged_df) == "LENGTE_216")
    df_new[id, 3] <- 4
  } else if(merged_df$Age[merged_df$GKNO == idFiveYounger[id]] == "3"){ # Lengte en gewicht op 3 jaar
    index_weight <- which(colnames(merged_df) == "GEWICHT_215")
    index_height <- which(colnames(merged_df) == "LENGTE_215")
    df_new[id, 3] <- 3
  } else {
    younger <- c(younger, which(merged_df$GKNO == idFiveYounger[id]))
    next
  }
  
  if(is.na(merged_df[merged_df$GKNO == idFiveYounger[id], index_weight])){
    Rweight <- c(Rweight, idFiveYounger[id])
  } else{
    df_new[id, 4] <- merged_df[merged_df$GKNO == idFiveYounger[id], index_weight]
  }
  if(is.na(merged_df[merged_df$GKNO == idFiveYounger[id], index_height])){
    Rheight <- c(Rheight, idFiveYounger[id])
  } else{
    df_new[id, 5] <- merged_df[merged_df$GKNO == idFiveYounger[id], index_height]
    df_new[id, 6] <- merged_df$sex[merged_df$GKNO == idFiveYounger[id]]
  }
}
colnames(df_new) <- c("GKNO", "filename", "age", "weight_gr", "height_cm", "gender")


#Remove ids with missing values for height and weight 
Rbmi <- union(Rheight, Rweight)
Rage <- df_new$GKNO[which(is.na(df_new$age))]
Rmissings <- union(Rbmi, Rage)
remove <- which(df_new$GKNO %in% Rmissings)
if(length(Rmissings > 0)){
  df_new <- df_new[-remove,]
} else {
  df_new <- df_new
}
completePP <- df_new$GKNO

completeData <- merged_df[which(merged_df$GKNO %in% completePP), ]
df_new$gender <- factor(df_new$gender, levels = c(1, 2), labels = c("M", "F"))
write.csv(df_new, file = paste0(spssanthrodir, "/anthro_acc_measurement_minimal_format.csv"))
write_sav(completeData, paste0(spssanthrodir, "/anthro_acc_measurement_original_format.sav"))

rm(accspss, merged_df, spssDataAcc)
rm(file, id, index_height, index_weight, Rbmi, Rheight, Rweight, uniquepp, height, weight, younger, remove,GECKOid, idFiveYounger)

# Copy accelometer files with data on BMI to new folder: completedata
filestocopy <- unlist(strsplit(df_new$filename, "15sec.agd"))

for (file in 1:length(filestocopy)) {
  if(df_new$age[file] == "5"){
    newfolder = paste0(dir, "Accelerometer data/rawinput/csv/completedata/5years/")
  } else if(df_new$age[file] == "4"){
    newfolder = paste0(dir, "Accelerometer data/rawinput/csv/completedata/4years/")
  } else {
    newfolder = paste0(dir, "Accelerometer data/rawinput/csv/completedata/3years/")
  }
  if(file.exists(paste0(dir, "Accelerometer data/rawinput/csv/", filestocopy[file], ".csv"))){
    file.copy(from = paste0(dir, "Accelerometer data/rawinput/csv/", filestocopy[file], ".csv"),
              to = paste0(newfolder, filestocopy[file], ".csv"), overwrite = FALSE)
  }
}

# b) Load raw data and calculate metrics for 15-sec epochs separate for each age group, as different cut-points are needed for MVPA
GGIR::GGIR(
  datadir = filedir,
  outputdir = "C://Users//P070400//OneDrive - Amsterdam UMC//Documenten", # Save the raw data and calculated metrics here
  mode = c(1, 2),
  do.report = c(2),
  ## Part 1 –  data processing
  #windowsizes = c(15, 900, 5400), # Epoch length: 15 sec, non-wear time 90min*60 = 5400
  windowsizes = c(15, 900, 1200), # Epoch length: 15 sec, non-wear time 20min*60 = 1200
  minimumFileSizeMB = 0.1,
  # Metrics
  do.enmo = TRUE,
  do.mad = TRUE,
  do.neishabouricounts = TRUE,
  do.hfen = TRUE,
  do.hfenplus	= TRUE,
  do.parallel = FALSE,
  acc.metric = "NeishabouriCount_y",
  ## Part 2 – data quality and descriptive
  qlevels <- c(c(1380/1440), c(1410/1440), c(1430/1140)),
  includedaycrit = 10,
  mvpathreshold = cutpoint,
  overwrite = FALSE #do not overwrite data
)

## c) Process accelerometer files to epoch level and store valid data only
# Note: I copied the preprocessed accelerometer data to: the newdatafolder
files <- list.files(path = newdatafolder, pattern = ".RData")

calibration_error_na <- 0
calibration_error_bigger <- 0

for(file in 1:length(files)){
  
  if((file.exists(paste0(epochdir, files[file])) & overwrit == FALSE)){
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
      save(epochdata, file = paste0(epochdir, files[file]))
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
