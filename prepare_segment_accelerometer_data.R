# Tabula rasa
rm(list=ls())
gc()

## This script: 
### STEP 1: Prepares the accelerometer data, by 
## a) Selecting the accelerometer data files for which data is complete (child is >=3 -years-old <=5 and data on BMI is available (copy this data to new folder: completedata)
## b) Load raw data and calculate metrics for 15-sec epochs separate for each age group, as different cut-points are needed for MVPA
## c) Process accelerometer files to epoch level and store valid data only

### STEP 2: Segment the accelerometer data, by 
## a) Generating the mhsmm format for the HFENplus variable
## b) Fitting hidden semi-Markov models for each file
## c) Copy the trained hsmms for the least and most active subjects to a folder

# Libraries
library(haven)
library(GGIR)
library(actilifecounts)
library(mhsmm)

## User settings
# Step 1:
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
hoursValidDay <- 10
minDays <- 3

# Step 2:
sequencedir <- "C://Users//P070400//OneDrive - Amsterdam UMC//Documenten/sequences/"
nStates <- 3 #For now 3, if we train multiple states look at: ~/Documents/Werk/PROGRAMMING/script/hsmm_functions.R
qMax <- 15 # Maximum number of iterations
dMax <- (60*60) / epoch  # Maximum state duration of 240 15-sec epochs, corresponds to 60 minutes - as was set in van Kuppevelt et al. (2019)
boutDurations <- c(30, 10, 5) # Note that the duration of the first bout has to be non-zero (as lambda is required to be > 0)
variable_name <- "HFENplus" #metric for training hsmm

### STEP 1) Prepare the accelerometer data: 
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
write.csv(df_new, file = paste0(spssanthrodir, "anthro_acc_measurement_minimal_format.csv"))
write_sav(completeData, paste0(spssanthrodir, "anthro_acc_measurement_original_format.sav"))

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

# Check if there are still double files
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
    #mhsmmdata <- formatMhsmm(validData, variable_name = "ENMO")
    #mhsmmdata <- formatMhsmm(validData,  "ENMO")
    #acceleration.descriptives[f, 1]$id <- filename
    #mhsmmdata <- formatMhsmm(validData, variablename = "NeishabouriCount_vm")
    mhsmmdata$Y <- asinh(mhsmmdata$Y*1000) # Scale the observations 
    mhsmmdata$Y[which(mhsmmdata$Y == 0)] <- rnorm(n = sum(mhsmmdata$Y == 0), mean = 0, sd = 0.05) # Add a little noise to zeros
    save(file = paste0(paste0(validdatadir, "mhsmmdata/"), filename), mhsmmdata) 
}

## b) Fit Hsmms
acceleration.descriptives <- data.frame()
files <- list.files(path = paste0(validdatadir, "mhsmmdata/"), pattern = ".RData")

for(f in 1:length(files)){
  cat(f)
  load(paste0(paste0(validdatadir, "mhsmmdata/"), files[f]))
  filename <- files[f]

  acceleration.descriptives[f,"id"] <- strsplit(filename, "_")[[1]][1]
  acceleration.descriptives[f,"min"] <- min(mhsmmdata$Y)
  acceleration.descriptives[f,"max"] <- max(mhsmmdata$Y)
  acceleration.descriptives[f,"Q1"] <-  quantile(mhsmmdata$Y, 0.25)
  acceleration.descriptives[f,"Median"] <-  quantile(mhsmmdata$Y, 0.5)
  acceleration.descriptives[f,"Q3"] <-  quantile(mhsmmdata$Y, 0.75)
  
  if(!file.exists(paste(sequencedir, filename, sep="/"))){
    # Fit the hsmm
    initParams <- derive_initialparams(mhsmmdata$Y, nStates, boutDurations, epoch)
    startmodel <- hsmmspec(initParams$initial, initParams$transition, initParams$emission, 
                           initParams$sojourn_distribution, dens.emis = dnorm.hsmm)
    hsmms <- hsmmfit(mhsmmdata$Y, startmodel, mstep = mstep.norm, maxit = qMax, M = dMax) 
    save(file = paste(sequencedir, filename, sep="/"), hsmms)
  } else{
    load(paste(sequencedir, filename, sep="/"))
  }
  
  # # Plot state space and sojourn distributions of the states
  colPalette <- scales::hue_pal()(hsmms$J) # Set the color palette depending on the number of states
  # png(paste0(sequencedir, "/plots/" , filename, "_hsmm_state_information.png"))
  # visualizeHSMM(hsmms, colPalette) # Save the plot with the model information
  # dev.off()
  # png(paste0(sequencedir, "/plots/" , filename, "hsmm_sojourn.png"))
  # plotSoj(hsmms, colPalette) # Save the plot with the sojourn densities
  # dev.off()
  # 
  # overlay data with estimated hidden states
  stateMeans <- hsmms$model$parms.emission$mu
  stateQ975s <- hsmms$model$parms.emission$mu + 2*sqrt(hsmms$model$parms.emission$sigma)
  stateQ025s <- hsmms$model$parms.emission$mu - 2*sqrt(hsmms$model$parms.emission$sigma)
  stateSummaries <- cbind(stateMeans, stateQ025s, stateQ975s)
  rownames(stateSummaries) <- paste("state", c(1:nStates), sep="")
  colnames(stateSummaries) <- c("mean", "min2sd", "plus2sd")
  
  # check residuals
  hist(mhsmmdata$Y - stateMeans[hsmms$yhat], n=100, col="blue", border="lightblue")
  qqnorm(mhsmmdata$Y - stateMeans[hsmms$yhat], pch=20)
  
  # find segments for nicer plotting
  segments <- numeric()
  init     <- 1
  for (j in 2:length(hsmms$yhat)){
    if (hsmms$yhat[j] != hsmms$yhat[j-1]){
      segments <- rbind(segments, c(init, j-1, hsmms$yhat[j-1]))
      init <- j
    }
  }
  png(paste0(sequencedir, "plots/" , filename, "overlay_data_states.png"))
  plot(mhsmmdata$Y, pch=".", cex=2)
  for (k in 1:nrow(segments)){
    ids <- segments[k,1]:segments[k,2]
    lines(stateMeans[hsmms$yhat[ids]] ~ ids, col=colPalette[segments[k,3]], lwd=2)
  }
  dev.off()
  
  # extract transition probabilities
  transitions <- round(hsmms$model$transition, 3)
  rownames(transitions) <- paste("state", c(1:nStates), sep="")
  colnames(transitions) <- paste("state", c(1:nStates), sep="")
  
  # calculate stationary distribution
  statDist <- transitions
  for (u in 1:1000){
    statDist <- statDist %*% transitions
  }
  statDist <- statDist[1,]
  
  # time spent in each state
  stateTimes <- matrix(ncol=6, nrow=nStates)
  rownames(stateTimes) <- paste("state", c(1:nStates), sep="")
  colnames(stateTimes) <- c("obs.time (%)", "min.length", "1st.quartile", "median", "3rd.quartile", "max.length")
  for (k in 1:nStates){
    stateTimes[k,1] <- sum(hsmms$yhat == k)
  }
  stateTimes[,1]  <- round(100*stateTimes[,1] / sum(stateTimes[,1]), 3)
  timeSpent   <- segments[,2] - segments[,1] + 1
  for (k in 1:nStates){
    stateTimes[k,-1] <- quantile(timeSpent[which(segments[,3]==k)], c(0, 0.25, 0.5, 0.75, 1))
  }
  
  # print summary
  sink(paste0(sequencedir, "summary/" , filename, "summaryStates_profile.txt"))
  cat(paste("number of identified states: ", nStates, sep=""), "\n")
  cat("\n")
  cat(paste("acceleration statistics of states: ", sep=""), "\n")
  print(round(stateSummaries, 3), quote=FALSE)
  cat("\n")
  cat("\n")
  cat(paste("state transition probabilities: ", sep=""), "\n")
  print(round(transitions, 3), quote=FALSE)
  cat("\n")
  cat("\n")
  cat(paste("number of state changes: ", nrow(segments), sep=""), "\n")
  cat("x\n")
  cat(paste("time spent in each state: ", sep=""), "\n")
  print(stateTimes, quote=FALSE)
  sink()
}
write.csv(acceleration.descriptives, file = paste0(sequencedir, "mhsmmdata/descriptives.csv"))

## c) Make a selection of the files of the most (Q4) and least (Q1) active participants and copy the fitted hsmm to a folder
# using total time spent in Daily MVPA (Sirard Cut-point for Y-axis): MVPA_E15S_T890_NeishabouriCount_y_0-24hr 
ggir_summary <- read.csv("C://Users//P070400//OneDrive - Amsterdam UMC//Documenten/output_5years90cts0/results/part2_daysummary.csv", sep = ",")

pp <- c()
for(file in 1:length(files)){
  pp <- c(pp, strsplit(files[file], ".RData")[[1]][1])
}

MVPAdaily <- c()
ppDaily <- c()
for(row in 1:nrow(ggir_summary)){
  if(strsplit(ggir_summary$ID[row], " ")[[1]][1] %in% pp){
    if(ggir_summary$N.valid.hours[row] >= 10){
      MVPAdaily <- c(MVPAdaily, ggir_summary$MVPA_E15S_T890_NeishabouriCount_y_0.24hr[row])
      ppDaily <- c(ppDaily, strsplit(ggir_summary$ID[row], " ")[[1]][1])
    }
  }
}

Q1 <- ifelse(MVPAdaily <= quantile(MVPAdaily)[2], 1, 0)
Q4 <- ifelse(MVPAdaily >= quantile(MVPAdaily)[4], 1, 0)

MVPAdf <- as.data.frame(cbind(ppDaily, MVPAdaily, Q1, Q4))
ppQ4 <- c() # If a participant is in Q4 during all measurement days
ppQ1 <- c() # If a participant is in Q1 during all measurement days
pp.active <- c() # If a participant is in Q4 during more than half of the measurement days
pp.sedentary <- c() # If a participant is in Q1 during more than half of the measurement days

for(subj in 1:length(unique(MVPAdf$ppDaily))){
  tmp <- MVPAdf[which(MVPAdf$ppDaily == unique(MVPAdf$ppDaily)[subj]),]
  if(sum(as.numeric(tmp$Q4)) == nrow(tmp)){
    ppQ4 <- c(ppQ4, unique(MVPAdf$ppDaily)[subj])
    if(!file.exists(paste0(sequencedir, "most active/", unique(MVPAdf$ppDaily)[subj], ".RData"))){
      file.copy(from = paste0(sequencedir, unique(MVPAdf$ppDaily)[subj], ".Rdata"),
                to = paste0(sequencedir, "most active/", unique(MVPAdf$ppDaily)[subj], ".RData"), overwrite = FALSE)
    }
  }
  if(sum(as.numeric(tmp$Q1)) == nrow(tmp)){
    ppQ1 <- c(ppQ1, unique(MVPAdf$ppDaily)[subj])
    if(!file.exists(paste0(sequencedir, "least active/", unique(MVPAdf$ppDaily)[subj], ".RData"))){
      file.copy(from = paste0(sequencedir, unique(MVPAdf$ppDaily)[subj], ".Rdata"),
                to = paste0(sequencedir, "least active/", unique(MVPAdf$ppDaily)[subj], ".RData"), overwrite = FALSE)
    }
  }
  if(sum(as.numeric(tmp$Q4)) >= nrow(tmp)/2){
    pp.active <- c(pp.active, unique(MVPAdf$ppDaily)[subj])
    if(!file.exists(paste0(sequencedir, "active/", unique(MVPAdf$ppDaily)[subj], ".RData"))){
      file.copy(from = paste0(sequencedir, unique(MVPAdf$ppDaily)[subj], ".Rdata"),
                to = paste0(sequencedir, "active/", unique(MVPAdf$ppDaily)[subj], ".RData"), overwrite = FALSE)
    }
  }
  if(sum(as.numeric(tmp$Q1)) >= nrow(tmp)/2){
    pp.sedentary <- c(pp.sedentary, unique(MVPAdf$ppDaily)[subj])
    if(!file.exists(paste0(sequencedir, "sedentary/", unique(MVPAdf$ppDaily)[subj], ".RData"))){
      file.copy(from = paste0(sequencedir, unique(MVPAdf$ppDaily)[subj], ".Rdata"),
                to = paste0(sequencedir, "sedentary/", unique(MVPAdf$ppDaily)[subj], ".RData"), overwrite = FALSE)
    }
  }
}

