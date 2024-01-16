## TO DO: rewrite and code documentation

# Tabula rasa
rm(list=ls())
gc()

# user input required:
dataDir <- "/Users/annelindelettink/Documents/Work MacBook Pro Annelinde/My Little Moves (MLM)/Sequence mapping/Physical behavior patterns/sequence-probability"
epochDataDir <- paste0(dataDir, "/data/epochlevelinput/validData")
sequenceDir <- paste0(dataDir, "/data/sequences")

# Load all packages
packages = c("mhsmm")
lapply(packages, FUN = function(X) {
  tmp = unlist(strsplit(X, "/"))
  X = tmp[length(tmp)]
  do.call("require", list(X)) 
})

# Load all scripts
my_functions_folder = "~/Documents/Work MacBook Pro Annelinde/My Little Moves (MLM)/Sequence mapping/Physical behavior patterns/sequence-probability/functions"
for (function_file in dir(my_functions_folder, full.names = T)) source(function_file) #load functions

#----------------------------------------------------------------

# For each participant, segment the data using a hidden semi-Markov model
files <- list.files(path = epochDataDir) #List data files

# General parameters
epoch <- 15 # epoch length (sec)
nStates <- 4 #For now 4, if we train multiple states look at: ~/Documents/Werk/PROGRAMMING/script/hsmm_functions.R
qMax <- 15 # Maximum number of iterations
dMax <- (60*60) / epoch  # Maximum state duration of 240 15-sec epochs, corresponds to 60 minutes - as was set in van Kuppevelt et al. (2019)
boutDurations <- c(30, 10, 5, 1) # Note that the duration of the first bout has to be non-zero (as lambda is required to be > 0)

# Load file
files <- list.files(epochDataDir)

for(f in 1:length(files)){
  cat(f)
  load(paste0(epochDataDir, "/", files[f]))
  filename = files[f]
  
  #format data for hsmm
  mhsmmdata <- formatMhsmm(validData, variablename = "NeishabouriCount_vm")
  mhsmmdata$Y <- asinh(mhsmmdata$Y*1000) # Scale the observations 
  mhsmmdata$Y[which(mhsmmdata$Y == 0)] <- rnorm(n = sum(mhsmmdata$Y == 0), mean = 0, sd = 0.05) # Add a little noise to zeros
  save(file = paste0(paste0(sequenceDir, "/mhsmmdata/"), filename), mhsmmdata)
  
  # Fit the hsmm
  initParams <- deriveInitialParameters(mhsmmdata$Y, nStates, boutDurations, epoch)
  startmodel <- hsmmspec(initParams$initial, initParams$transition, initParams$emission, 
                         initParams$sojourn_duration$sojourn_poisson, dens.emis = dnorm.hsmm)
  hsmms <- hsmmfit(mhsmmdata$Y, startmodel, mstep = mstep.norm, maxit = qMax, M = dMax) 
  save(file = paste(sequenceDir, filename, sep="/"), hsmms)
  
  # Plot state space and sojourn distributions of the states
  colPalette <- scales::hue_pal()(hsmms$J) # Set the color palette depending on the number of states
  png(paste0(sequenceDir, "/plots/" , filename, "_hsmm_state_information.png"))
  visualizeHSMM(hsmms, colPalette) # Save the plot with the model information
  dev.off()
  png(paste0(sequenceDir, "/plots/" , filename, "hsmm_sojourn.png"))
  plotSoj(hsmms, colPalette) # Save the plot with the sojourn densities
  dev.off()
  
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
  png(paste0(sequenceDir, "/plots/" , filename, "overlay_data_states.png"))
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
  sink(paste0(sequenceDir, "/summary/" , filename, "summaryStates_profile.txt"))
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
  cat("\n")
  cat(paste("time spent in each state: ", sep=""), "\n")
  print(stateTimes, quote=FALSE)
  sink()
}
