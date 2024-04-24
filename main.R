# Tabula rasa
rm(list=ls())
gc()

# user input required:
datadir <- "/Users/annelindelettink/Documents/Work MacBook Pro Annelinde/My Little Moves (MLM)/Sequence mapping/Physical behavior patterns"
sequencedir <- paste0(datadir, "/data/sequences")
spssdir <- "/Users/annelindelettink/GECKO/deel1"

# Load packages
library(mhsmm)
library(haven)
library(anthro)

# Load all scripts
# This can be replaced by the package install later on...
my_functions_folder = "/Users/annelindelettink/Documents/Work MacBook Pro Annelinde/My Little Moves (MLM)/Sequence mapping/Physical behavior patterns/motif-probability/R/"
for (function_file in dir(my_functions_folder, full.names = TRUE)) source(function_file) #load functions

### 1 ) Define motif (i.e. a sequence of events, each characterized by their acceleration range and length)
# Check hypotheses (proof of concept): 
#descriptives <- read.csv(paste0(sequencedir, "//descriptives.csv"), sep = ",")[,-1] # acceleration descriptive statistics

# Define epoch length
epoch_length <- 15
n_epochs_min <- 60 / epoch_length

# Function to create motif data frame
create_motif <- function(Amin, Amax, length, upper_bound) {
  return(data.frame(Amin = Amin, Amax = Amax, length = length, upper_bound = upper_bound))
}

# Function to create motifs for PA bouts
create_PA_motifs <- function(bout_lengths, upper_bounds, n_epochs_min) {
  motifs <- list()
  for (i in seq_along(bout_lengths)) {
    motifs[[paste0("motif.PA_", bout_lengths[i], "min")]] <- create_motif(Amin.PA, Amax.PA, bout_lengths[i] * n_epochs_min, upper_bounds[i])
  }
  return(motifs)
}

# Function to create motifs for SB bouts
create_SB_motifs <- function(bout_lengths, upper_bounds, n_epochs_min) {
  motifs <- list()
  for (i in seq_along(bout_lengths)) {
    motifs[[paste0("motif.SB_", bout_lengths[i], "min")]] <- create_motif(Amin.SB, Amax.SB, bout_lengths[i] * n_epochs_min, upper_bounds[i])
  }
  return(motifs)
}

# Define acceleration ranges for PA and SB
Amin.PA <- 2.5
Amax.PA <- 7
Amin.SB <- 0
Amax.SB <- 0.5

# Define lengths and upper bounds for PA and SB bouts
PA_bout_lengths <- c(0.5, 1, 5, 10, 15)
PA_upper_bounds <- c(1, 5, 10, 15, Inf) # Upper bounds for PA bout lengths
SB_bout_lengths <- c(10, 20, 30)
SB_upper_bounds <- c(20, 30, Inf) # Upper bounds for SB bout lengths

# Create motifs for PA and SB bouts
PA_motifs <- create_PA_motifs(PA_bout_lengths, PA_upper_bounds)
SB_motifs <- create_SB_motifs(SB_bout_lengths, SB_upper_bounds)

# Combine motifs for SB bouts interrupted and uninterrupted
motif.SB_30min.interrupted <- do.call(rbind, replicate(3, SB_motifs[["motif.SB_10min"]], simplify = FALSE))

# Print motifs
print(PA_motifs)
print(SB_motifs)
print(motif.SB_30min.interrupted)

### 2) Apply forward algorithm for motif probability calculation

# Define function to calculate probabilities for motifs
calculate_probabilities <- function(filelist, sequencedir, motif_definitions) {
  # Create empty data frame to store probabilities
  probabilities <- data.frame(id = character())
  for (pp in 1:length(filelist)) {
    # Load the fitted hsmm
    load(paste0(sequencedir, "/", filelist[pp]))

    # Extract participant ID
    participant_id <- strsplit(filelist[pp], "_")[[1]][1]
    # Initialize list to store probabilities for each motif
    motif_probabilities <- list()
    for (motif_name in names(motif_definitions)) {
      motif_probabilities[[motif_name]] <- motif_probability(motif_definitions[[motif_name]], hsmms)
    }
    
    # Append calculated probabilities to the data frame
    probabilities <- rbind(probabilities, c(id = participant_id, as.numeric(unlist(motif_probabilities))))
  }
  
  return(probabilities)
}


### 3) Link zBMI 
spssdata <- read_sav(paste(spssdir, "Data_GECKO_MyLittleMoves_Deel1_2020-10-22.sav", sep = "/")) # Antropometrics
spss_acc <- spssdata[which(spssdata$GKNO %in% probabilities.heuristic$id),]
rm(spssdata)
# Calculate zBMI based on WHO standards
spss_acc$sex_new <- ifelse(spss_acc$sex==1, "M", "F")
#TO DO: if age is missing (spss_acc$measure_age_217) impute mean age check volgende regel!
spss_acc$measure_age_217[which(is.na(spss_acc$measure_age_217))] <- mean(spss_acc$measure_age_217, na.rm = TRUE)
antro <- anthro_zscores(sex = spss_acc$sex_new, age = spss_acc$measure_age_217, weight = (spss_acc$GEWICHT_217/1000), lenhei = spss_acc$LENGTE_217)
spss_acc$bmi <- (spss_acc$GEWICHT_217/1000)/((spss_acc$LENGTE_217/100)^2)
library(dplyr)
antro <- antro %>%
  mutate(overweight = factor(case_when(
    zbmi > 2 ~ 2,
    zbmi > 1 & zbmi <= 2 ~ 1,
    zbmi < -1 ~ -1,
    TRUE ~ 0
  )))

waist <- c()
bmi <- c()
zBMI <- c()
SDQ <- c()
gender <- c()
for (pp in 1:nrow(probabilities.heuristic)){
  row <- which(spss_acc$GKNO == probabilities.heuristic$id[pp])
  gender <- c(gender, spss_acc$sex[row])
  zBMI <- c(zBMI, antro$zbmi[row])
  bmi <- c(bmi, spss_acc$bmi[row])
  waist <- c(waist, spss_acc$BUIKOM_217_QC[row])
  SDQ <- c(SDQ, spss_acc$SDQ_HYP_217[row])
}
probabilities.heuristic <- cbind(probabilities.heuristic, gender, zBMI, bmi, waist, SDQ)
probabilities.heuristic$gender <- ifelse(probabilities.heuristic$gender == 1, "M", ifelse(probabilities.heuristic$gender == 2, "F", "NA"))


#colnames(probabilities.heuristic) <- c("id", "prob_long","prob_alternating", "zBMI", "waist", "WHR")
save(probabilities.heuristic, file = paste0(datadir, "/data/motif-probabilities.RData"))

# Correlations with zBMI to test hypotheses
# higher probability of longer uninterrupted activity is negatively related to BMI
cor.test(probabilities.heuristic$zBMI, exp(as.numeric(probabilities.heuristic$PA.05min)), na.action = "omit", method = "kendall")
cor.test(probabilities.heuristic$zBMI, exp(as.numeric(probabilities.heuristic$PA.1min)), na.action = "omit", method = "kendall")
cor.test(probabilities.heuristic$zBMI, exp(as.numeric(probabilities.heuristic$PA.5min)), na.action = "omit", method = "kendall")
cor.test(probabilities.heuristic$zBMI, exp(as.numeric(probabilities.heuristic$PA.10min)), na.action = "omit", method = "kendall")
cor.test(probabilities.heuristic$zBMI, exp(as.numeric(probabilities.heuristic$PA.15min)), na.action = "omit", method = "kendall")

# higher probability of longer uninterrupted sedentary behavior is positively related to BMI -- NOT TRUE
cor.test(probabilities.heuristic$zBMI, exp(as.numeric(probabilities.heuristic$SB.10min)), na.action = "omit", method = "kendall")
cor.test(probabilities.heuristic$zBMI, exp(as.numeric(probabilities.heuristic$SB.20min)), na.action = "omit", method = "kendall")
cor.test(probabilities.heuristic$zBMI, exp(as.numeric(probabilities.heuristic$SB.30min)), na.action = "omit", method = "kendall")

# higher probability of longer uninterrupted activity different for boys vs girls?
cor.test(probabilities.heuristic$zBMI, exp(as.numeric(probabilities.heuristic$PA.1min)), na.action = "omit")
cor.test(probabilities.heuristic$zBMI, exp(as.numeric(probabilities.heuristic$PA.5min)), na.action = "omit")
cor.test(probabilities.heuristic$zBMI, exp(as.numeric(probabilities.heuristic$PA.10min)), na.action = "omit")
cor.test(probabilities.heuristic$zBMI, exp(as.numeric(probabilities.heuristic$PA.15min)), na.action = "omit")

# Differences between boys and girls?
boys <- subset(probabilities.heuristic, gender == "M")
girls <- subset(probabilities.heuristic, gender == "F")

# longer uninterrupted activity 
wilcox.test(girls$PA.1min, boys$PA.1min, na.action = "omit")
wilcox.test(girls$PA.5min, boys$PA.5min, na.action = "omit")
wilcox.test(girls$PA.10min, boys$PA.10min, na.action = "omit")

wilcox.test(girls$SB.10min, boys$SB.10min)
wilcox.test(girls$SB.20min, boys$SB.20min)
wilcox.test(girls$SB.30min, boys$SB.30min)

wilcox.test(girls$SB.30min.interrupted, boys$SB.30min.interrupted)

cor.test(probabilities.heuristic$zBMI, exp(as.numeric(probabilities.heuristic$PA.5min)), na.action = "omit")
cor.test(probabilities.heuristic$zBMI, exp(as.numeric(probabilities.heuristic$PA.10min)), na.action = "omit")
cor.test(probabilities.heuristic$zBMI, exp(as.numeric(probabilities.heuristic$PA.15min)), na.action = "omit")



t.test(girls$PA.10min, boys$PA.10min, na.action = "omit", paired = FALSE)

t.test(girls$SB.30min.interrupted, boys$SB.30min.interrupted, na.action = "omit", paired = FALSE)


cor.test(probabilities.heuristic$zBMI, exp(as.numeric(probabilities.heuristic$SB.30min.interrupted)), na.action = "omit")


lm(probabilities.heuristic$zBMI ~ exp(as.numeric(probabilities.heuristic$SB.30min.interrupted)), na.action = "omit")


# Heatmaps
probabilities$prob_long[which(probabilities$prob_long == "-Inf")] <- 0
probabilities$prob_alternating[which(probabilities$prob_alternating == "-Inf")] <- 0
probabilities$prob_long <- as.numeric(probabilities$prob_long)
probabilities$prob_alternating <- as.numeric(probabilities$prob_alternating)



plot(prob_long/prob_alternating)

# heatmap <- ggplot2::ggplot(probabilities, ggplot2::aes(prob_long, prob_alternating, fill = zBMI)) +
#   ggplot2::geom_tile() 
scatter_plot <- ggplot2::ggplot(probabilities, ggplot2::aes(zBMI, (prob_long/prob_alternating))) +
  ggplot2::geom_point(color = "blue")

scatter_plot_long <- ggplot2::ggplot(probabilities, ggplot2::aes(zBMI, prob_long)) +
  ggplot2::geom_point(color = "blue")

scatter_plot_alternating <- ggplot2::ggplot(probabilities, ggplot2::aes(zBMI, prob_alternating)) +
  ggplot2::geom_point(color = "blue")


scatter_plot_long_whr <- ggplot2::ggplot(probabilities, ggplot2::aes(WHR, prob_long)) +
  ggplot2::geom_point(color = "blue")

scatter_plot_alternating_whr <- ggplot2::ggplot(probabilities, ggplot2::aes(WHR, prob_alternating)) +
  ggplot2::geom_point(color = "blue")



hist(prob_long/prob_alternating)

spss_acc$GEWICHT_215/(spss_acc$GEWICHT_215^2)

