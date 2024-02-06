## This script: 
# 1) Converts the old gt3x files to csv for compatibility with GGIR;
# 1) Selects the raw accelerometer data files for which data on BMI is available; 
# 3) Selects valid data only

# Tabula rasa
rm(list=ls())
gc()

# Load all scripts
my_functions_folder = "~/Documents/Work MacBook Pro Annelinde/My Little Moves (MLM)/Sequence mapping/Physical behavior patterns/motif-probability/R"
for (function_file in dir(my_functions_folder, full.names = T)) source(function_file) #load functions

# User settings
spssdir <- "/Users/annelindelettink/GECKO/deel1"
accdir <- "/Users/annelindelettink/GECKO/rawinput"


# Libraries
library(haven)
library(GGIR)

### STEP 1) Convert the old gt3x files to csv data
## TO DO: hiervoor nog de stap omzetten oude gt3x naar csv so GGIR can handle the files
list.gt3x <- list.files("/Users/annelindelettink/GECKO/rawinput", pattern = ".gt3x")
for(file in 1:length(list.gt3x)){
  gt3x_to_csv(paste0("/Users/annelindelettink/GECKO/rawinput", "/", list.gt3x[file]))
}

### STEP 2) Prepare the raw accelerometer data files for which data on BMI is available; 
spssdata <- read_sav(paste(spssdir, "Data_GECKO_MyLittleMoves_Deel1_2020-10-22.sav", sep = "/")) #antropometrics

# Select antrophometrics for which accelerometer data is available
#accFilesgt3x <- list.files(accDirGT3X, pattern = ".gt3x")
accfiles <- list.files(accdir, pattern = ".csv")
GECKOid <- c() # Abbreviate participant id from accelerometer data to correspond with SPSS data
for (file in 1:length(accfiles)) {
  GECKOid <- c(GECKOid, strsplit(accfiles[file], "_")[[1]][1])
}
spssDataAcc <- spssdata[which(spssdata$GKNO %in% GECKOid),] # select the data
rm(spssdata)

# Check which participants have complete data for BMI at age 5
Rbmi <- union(which(is.na(spssDataAcc$LENGTE_217)), which(is.na(spssDataAcc$GEWICHT_217)))
if(length(Rbmi > 0)){
  spssDataAccBMI <- spssDataAcc[-Rbmi,]
} else {
  spssDataAccBMI <- spssDataAcc
}
rm(spssDataAcc)
completePP <- spssDataAccBMI$GKNO

#Copy acc files with data on BMI to new folder: rawinput
filestocopy <- c()
for (accfile in 1:length(accfiles)) {
  index <- which(startsWith(accfiles[accfile], completePP) == TRUE)
  if(length(index) > 0){
    filestocopy <- c(filestocopy, accfiles[accfile])
  }
}

# Check which participants have multiple files (either double or multiple measurements)
multiple_days <- c() # measured multiple days; these data can be included 
double_file <- c() # measured twice or more for 1 measurement point

start_pp <- strsplit(strsplit(filestocopy[1], "_")[[1]][1], "-")[[1]][1] 
measurement <- strsplit(strsplit(filestocopy[1], "_")[[1]][1], "-")[[1]][2]
for(subj in 2:(length(filestocopy))) {
  next_pp <- strsplit(strsplit(filestocopy[subj], "_")[[1]][1], "-")[[1]][1]
  next_measurement <- strsplit(strsplit(filestocopy[subj], "_")[[1]][1], "-")[[1]][2]
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
remove <- c(remove, which(filestocopy == "1168-1_31-5-2013 (2013-06-27).csv")) # downloaded twice
remove <- c(remove, which(filestocopy == "1735-1_11APR12 (2012-05-29).csv")) # something went wrong with data, measured again; remove first try
remove <- c(remove, which(filestocopy == "3074-1_15sept11 (2011-10-14).csv")) # data parents and actilife do not correspond, therefore other day was measured
remove <- c(remove, which(filestocopy == "3473-1_28juni11 (2011-08-01).csv")) # went in washing machine, new accelerometer was sent
remove <- c(remove, which(filestocopy == "3960-1_06mei11 (2011-05-24).csv")) # data parents and actilife do not correspond and accelerometer probably worn by someone else, therefore other day was measured
remove <- c(remove, which(filestocopy == "4292-1_20juli11 (2011-08-05).csv")) # informed consent missing 

files <- c("1168", "1735", "3074", "3473", "3960", "4292")
double_file_corrected <- double_file[!double_file %in% files]
rm(files)

filestocopy <- filestocopy[-remove] # Remove the (dubious) files from the file list to ensure these are not loaded

# Copy input files to different folder and remove from current folder
newfolder <- "/Users/annelindelettink/GECKO/rawinput"
dir.create(newfolder, recursive = TRUE)
file.copy(file.path(accDirGT3X, filestocopy), newfolder)


# User settings
storedir <- "/Users/annelindelettink/Documents/Work MacBook Pro Annelinde/My Little Moves (MLM)/Sequence mapping/Physical behavior patterns/sequence-probability/data"

## Load raw data and calculate metrics for 15-sec epochs
GGIR::GGIR(
  datadir = accdir,
  outputdir = storedir, 
  mode = c(1, 2, 3),
  ## Part 1 – raw data processing 
  windowsizes = c(15, 900, 3600), # Epoch length: 15 sec
  # Metrics
  do.enmo = TRUE, 
  do.mad = TRUE,
  do.neishabouricounts = TRUE, 
  do.bfen = TRUE, 
  do.parallel = FALSE,
  ## Part 2 – data quality and descriptives
  includedaycrit = 10,
  #mvpathreshold = 100
  overwrite = FALSE #do not overwrite data
)

epochDataDir <- paste0(storeDir, "/output_rawinput/meta/ms2.out")
saveDir <- paste0(storeDir, "/epochlevelinput")
epochlevelFiles <- list.files(path = epochDataDir, pattern = ".RData")

for(pp in 1:length(epochlevelFiles)){
  load(paste(epochDataDir, epochlevelFiles[pp], sep = "/")) # Load RData file from meta/ms2.out folder
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
  IMP$metashort = cbind(IMP$metashort, scores) # combine scores with the epoch level timeseries
  epochdata <- IMP$metashort
  save(epochdata, file = paste(saveDir, epochlevelFiles[pp], sep = "/"))
}

### STEP 3) Select valid data only
storeDir <- paste0(saveDir, "/validData")

epoch = 15 # epoch length (sec)
epochMin = 60/epoch # number of epochs in one minute
hoursValidDay <- 10
minDays <- 3

files <- list.files(path = saveDir, pattern = ".RData")

for(file in 1:length(files)){
  load(paste(saveDir, files[file], sep = "/"))
  data_perday <- split(epochdata, as.Date(epochdata$timestampPOSIX))
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
      save(validData, file = paste(storeDir, filename, sep="/"))
    }
  }
}

### HIER NOG DUBBELE FILE CHECK TOEVOEGEN?
# CHECK IF THERE ARE STILL DOUBLE FILES (or were these removed when filtering the valid days?)
# #zie voorbeeldcode hieronder
# file_names <- names(GECKOvalid10h3d)
# f_name <- c()
# for (file in 1:length(file_names)) {
#   f_name <- c(f_name, strsplit(file_names[file], "_")[[1]][1])
# }
# file_names[which(duplicated(f_name))] # Double file for 2172-1 
# 
# # Search 2172-1 in validdayCriteria table shows that 2172-1_31mei2013 (2013-06-19)15sec.csv has 4 valid days, whereas, "2172-1_4sept12 (2012-09-18)15sec.csv" only 3
# remove <- which(duplicated(f_name)) # Remove the file with the least valid days
# GECKOvalid10h3d[[remove]] <- NULL
# 
# # Check if there are no more double files
# index <- c()
# for (f in 1:length(double_file_corrected)) {
#   index <- c(index, which(grepl(double_file_corrected[f], file_list_corrected)))
# }
# double_file_names <- file_list_corrected[index] # indicate the double files in the unfiltered GECKOdata1537
# 
# index <- c()
# for (f in 1:length(double_file_corrected)) {
#   index <- c(index, which(grepl(double_file_corrected[f], file_list_3d10h)))
# }
# fileCheck <- file_list_3d10h[index] 
# noMoreDoubles <- length(fileCheck) == length(unique(fileCheck))# check if there are no more doubles
# noMoreDoubles  #if TRUE then okay
# 
# # MULTIPLE MEASUREMENTS
# # Create a vector of the files that are associated with multiple measurements
# index <- c()
# for (f in 1:length(multiple_days)) {
#   index <- c(index, which(grepl(multiple_days[f], names(GECKOvalid10h3d))))
# }
# multiple_days_files <- names(GECKOvalid10h3d)[index] # indicate the file names of the multiple days
# 
# # Select the files with the most valid days, if these were equal, we selected the first measurement
# ##this was done by manually checking the validdayCriteria table
# # 1155, 1252, 1770, 1837, 1847, 2774, 3589: both measurements 4 valid days, so we removed the second measurement
# GECKOvalid10h3d$'1155-2_19mei11 (2011-06-01)15sec.csv' <- NULL
# GECKOvalid10h3d$'1252-2_2feb12 (2012-02-15)15sec.csv' <- NULL
# GECKOvalid10h3d$'1770-2_27feb12 (2012-03-07)15sec.csv' <- NULL
# GECKOvalid10h3d$'1837-2_21juni11 (2011-07-04)15sec.csv' <- NULL
# GECKOvalid10h3d$'1847-2_25okt11 (2011-11-10)15sec.csv' <- NULL
# GECKOvalid10h3d$'2774-2_19APR12 (2012-05-07)15sec.csv' <- NULL
# GECKOvalid10h3d$'3589-2_3APR12 (2012-04-17)15sec.csv' <- NULL
# #2766: both measurements 3 valid days, so we removed the second measurement
# GECKOvalid10h3d$'2766-2_15juni11 (2011-06-30)15sec.csv' <- NULL
