## This script: 
### STEP 1: Prepares the accelerometer data, by 
## a) Converting old gt3x files to csv data
## b) Loading raw data and calculating metrics for 15-sec epochs (using GGIR)
## c) Selecting the accelerometer data files for which data on BMI is available (copy this data to new folder: GECKOaccBMI)
## d) Process these files to epoch level and store valid data only, remove multiple diles per partcipants with least valid days

### STEP 2: Segments the accelerometer data by fitting hidden semi-Markov models for each file

# Tabula rasa
rm(list=ls())
gc()

# Load all scripts
my_functions_folder = "~/Documents/Work MacBook Pro Annelinde/My Little Moves (MLM)/Sequence mapping/Physical behavior patterns/motif-probability/R"
for (function_file in dir(my_functions_folder, full.names = T)) source(function_file) #load functions

# User settings
spssdir <- "/Users/annelindelettink/GECKO/deel1"
accdir <- "/Users/annelindelettink/GECKO/rawinput" # Save the .csv raw acceleration files here (STEP 1a)
storedir <- "/Users/annelindelettink/Documents/Work MacBook Pro Annelinde/My Little Moves (MLM)/Sequence mapping/Physical behavior patterns/data"

# Libraries
library(haven)
library(GGIR)

### STEP 1) Prepare the accelerometer data: 
## a) Convert the old gt3x files to csv data
# list.gt3x <- list.files(accdir, pattern = ".gt3x")
# for(file in 1:length(list.gt3x)){
#   gt3x_to_csv(paste0(accdir "/", list.gt3x[file]))
# }

## b) Load raw data and calculate metrics for 15-sec epochs
# GGIR::GGIR(
#   datadir = accdir,
#   outputdir = storedir, # Save the raw data and calculated metrics here 
#   mode = c(1, 2, 3),
#   ## Part 1 – raw data processing 
#   windowsizes = c(15, 900, 3600), # Epoch length: 15 sec
#   # Metrics
#   do.enmo = TRUE, 
#   do.mad = TRUE,
#   do.neishabouricounts = TRUE, 
#   do.bfen = TRUE, 
#   do.parallel = FALSE,
#   ## Part 2 – data quality and descriptives
#   includedaycrit = 10,
#   #mvpathreshold = 100
#   overwrite = FALSE #do not overwrite data
# )

## c) Copy accelerometer data for the participants with available data for BMI to new folder
newfolder <- paste0(storedir, "/GECKOaccBMI/")
spssdata <- read_sav(paste(spssdir, "Data_GECKO_MyLittleMoves_Deel1_2020-10-22.sav", sep = "/")) # Antropometrics
batch <- c(paste0(storedir, "/output_rawinput"), paste0(storedir, "/output_rawinput1"), paste0(storedir, "/output_rawinput2")) # Folders that hold all GGIR processed GECKO data
accfiles <- c()
for(folder in 1:length(batch)){
  accfiles <- c(accfiles, list.files(paste0(batch[folder], "/meta/ms2.out/"), pattern = ".RData")) # List accelerometer files
}

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

Rsex <- which(is.na(spssDataAcc$sex)) # complete data on gender?
if(length(Rsex > 0)){
  spssDataAccSEX <- spssDataAcc[-Rsex,]
} else {
  spssDataAccSEX <- spssDataAcc
}
rm(spssDataAcc)
completePP <- spssDataAccBMI$GKNO

# Copy accelometer files with data on BMI to new folder: GECKOaccBMI
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
remove <- c(remove, which(filestocopy == "1168-1_31-5-2013 (2013-06-27).csv.RData")) # downloaded twice
remove <- c(remove, which(filestocopy == "1735-1_11APR12 (2012-05-29).csv.RData")) # something went wrong with data, measured again; remove first try
remove <- c(remove, which(filestocopy == "3074-1_15sept11 (2011-10-14).csv.RData")) # data parents and actilife do not correspond, therefore other day was measured
remove <- c(remove, which(filestocopy == "3473-1_28juni11 (2011-08-01).csv.RData")) # went in washing machine, new accelerometer was sent
remove <- c(remove, which(filestocopy == "3960-1_06mei11 (2011-05-24).csv.RData")) # data parents and actilife do not correspond and accelerometer probably worn by someone else, therefore other day was measured
remove <- c(remove, which(filestocopy == "4292-1_20juli11 (2011-08-05).csv.RData")) # informed consent missing
if(length(remove)>0){
  filestocopy <- filestocopy[-remove] # Remove the (dubious) files from the file list to ensure these are not loaded
}
files <- c("1168", "1735", "3074", "3473", "3960", "4292")
double_file_corrected <- double_file[!double_file %in% files]
rm(files)

for (file in 1:length(filestocopy)) {
  # Iterate over each path
  for (b in 1:length(batch)) {
    if(file.exists(paste0(batch[b], "/meta/ms2.out/", filestocopy[file]))){
      file.copy(from = paste0(batch[b], "/meta/ms2.out/", filestocopy[file]),
      to = paste0(newfolder, filestocopy[file]), overwrite = FALSE)
    }
  }
}

## d) Process accelerometer files to epoch level and store valid data only
epochdir <- paste0(newfolder, "epochdata/")
validdatadir <- paste0(newfolder, "validdata/")
overwrit = FALSE
epoch = 15 # epoch length (sec)
epochMin = 60/epoch # number of epochs in one minute
hoursValidDay <- 10
minDays <- 3

files <- list.files(path = newfolder, pattern = ".RData")

for(file in 1:length(files)){
  
  if(file.exists(paste0(epochdir, files[file])) & overwrit == FALSE){
    print("epochdata already saved")
    load(paste(epochdir, files[file], sep = "/")) # Load epochdata
  } else {
    load(paste(newfolder, files[file], sep = "/")) # Load RData file
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
      save(validData, file = paste(validdatadir, filename, sep="/"))
    }
  }
}

# Check if there are still double files
validfiles <- list.files(validdatadir, pattern = ".RData")
f_name <- c()
for (file in 1:length(validfiles)) {
   f_name <- c(f_name, strsplit(validfiles[file], "_")[[1]][1])
}

validfiles[which(duplicated(f_name))] # Double files for 2094-1, 3496-1
load(paste0(validdatadir, "2094-1_13mrt12.RData"))
load(paste0(validdatadir, "2094-1_24APR12.RData")) # more valid days, so keep this one
validfiles <- validfiles[-which(validfiles == "2094-1_13mrt12.RData")]

# Create a vector of the files that are associated with multiple measurements
f_name[which(grepl("-2", f_name))]

load(paste0(validdatadir, "1155-1_19mei11.RData"))
load(paste0(validdatadir, "1155-2_19mei11.RData")) # Both equal number of days data, so keep the first
validfiles <- validfiles[-which(validfiles == "1155-2_19mei11.RData")]

load(paste0(validdatadir, "1770-1_27feb12.RData"))
load(paste0(validdatadir, "1770-2_27feb12.RData")) # Both equal number of days data, so keep the first
validfiles <- validfiles[-which(validfiles == "1770-2_27feb12.RData")]

load(paste0(validdatadir, "1837-1_21juni11.RData"))
load(paste0(validdatadir, "1837-2_21juni11.RData")) # Both equal number of days data, so keep the first
validfiles <- validfiles[-which(validfiles == "1837-2_21juni11.RData")]

load(paste0(validdatadir, "1847-1_26okt11.RData"))
load(paste0(validdatadir, "1847-2_25okt11.RData")) # Both equal number of days data, so keep the first
validfiles <- validfiles[-which(validfiles == "1847-2_25okt11.RData")]

load(paste0(validdatadir, "2774-1_19APR12.RData"))
load(paste0(validdatadir, "2774-2_19APR12.RData")) # Both equal number of days data, so keep the first
validfiles <- validfiles[-which(validfiles == "2774-2_19APR12.RData")]

load(paste0(validdatadir, "3497-1_7dec11.RData"))
load(paste0(validdatadir, "3497-2_7dec11.RData")) # Both equal number of days data, so keep the first
validfiles <- validfiles[-which(validfiles == "3497-2_7dec11.RData")]

load(paste0(validdatadir, "3589-1_3APR12.RData"))
load(paste0(validdatadir, "3589-2_3APR12.RData")) # Both equal number of days data, so keep the first
validfiles <- validfiles[-which(validfiles == "3589-2_3APR12.RData")]

### STEP 2: Segments the accelerometer data by fitting hidden semi-Markov models for each file
# Load all packages
packages = c("mhsmm")
lapply(packages, FUN = function(X) {
  tmp = unlist(strsplit(X, "/"))
  X = tmp[length(tmp)]
  do.call("require", list(X)) 
})

length(validfiles)
sequencedir <- paste0(storedir, "/sequences")


# General parameters
epoch <- 15 # epoch length (sec)
nStates <- 4 #For now 4, if we train multiple states look at: ~/Documents/Werk/PROGRAMMING/script/hsmm_functions.R
qMax <- 15 # Maximum number of iterations
dMax <- (60*60) / epoch  # Maximum state duration of 240 15-sec epochs, corresponds to 60 minutes - as was set in van Kuppevelt et al. (2019)
boutDurations <- c(30, 10, 5, 1) # Note that the duration of the first bout has to be non-zero (as lambda is required to be > 0)

acceleration.descriptives <- data.frame()

for(f in 1:length(validfiles)){
  cat(f)
  load(paste0(validdatadir, validfiles[f]))
  filename = validfiles[f]
  
  if(!file.exists(paste0(paste0(sequencedir, "/mhsmmdata/")))){
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
    acceleration.descriptives[f]$id <- 
    #mhsmmdata <- formatMhsmm(validData, variablename = "NeishabouriCount_vm")
    mhsmmdata$Y <- asinh(mhsmmdata$Y*1000) # Scale the observations 
    mhsmmdata$Y[which(mhsmmdata$Y == 0)] <- rnorm(n = sum(mhsmmdata$Y == 0), mean = 0, sd = 0.05) # Add a little noise to zeros
    save(file = paste0(paste0(sequencedir, "/mhsmmdata/"), filename), mhsmmdata)
  } else{
    load(paste0(paste0(sequencedir, "/mhsmmdata/"), filename))
  }
  
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
  png(paste0(sequencedir, "/plots/" , filename, "overlay_data_states.png"))
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
  sink(paste0(sequencedir, "/summary/" , filename, "summaryStates_profile.txt"))
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
write.csv(acceleration.descriptives, file = paste0(sequencedir, "/mhsmmdata/descriptives.csv"))

