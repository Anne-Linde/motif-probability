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

# Motif length: length from minutes to epoch level, lambda was based on 15-sec epoch time series
epoch_length = 15
n_epochs_min = 60/epoch_length 
# higher probability of longer uninterrupted activity is negatively related to BMI
L.PA_bouts_1_uninterrupted <- c(1) * n_epochs_min
L.PA_bouts_5_uninterrupted <- c(5) * n_epochs_min
L.PA_bouts_10_uninterrupted <- c(10) * n_epochs_min
L.PA_bouts_15_uninterrupted <- c(15) * n_epochs_min

# higher probability of longer uninterrupted sedentary behavior is positively related to BMI
L.SB_bouts_10_uninterrupted <- c(10) * n_epochs_min
L.SB_bouts_20_uninterrupted <- c(20) * n_epochs_min
L.SB_bouts_30_uninterrupted <- c(30) * n_epochs_min

# Motif: acceleration
Amin.PA = 2.5 #(mean(descriptives$Q3)+ mean(descriptives$Median))/2# Q3
Amax.PA = 7 #max(descriptives$max) 

Amin.SB = 0
Amax.SB = 0.5 #mean(descriptives$Q1)# Q1

# Motif definition
motif.PA_0.5min <- data.frame(Amin.PA, Amax.PA, 2)
colnames(motif.PA_0.5min) <- c("Amin", "Amax", "length")

motif.PA_1min <- data.frame(Amin.PA, Amax.PA, L.PA_bouts_1_uninterrupted)
colnames(motif.PA_1min) <- c("Amin", "Amax", "length")
motif.PA_5min <- data.frame(Amin.PA, Amax.PA, L.PA_bouts_5_uninterrupted)
colnames(motif.PA_5min) <- c("Amin", "Amax", "length")
motif.PA_10min <- data.frame(Amin.PA, Amax.PA, L.PA_bouts_10_uninterrupted)
colnames(motif.PA_10min) <- c("Amin", "Amax", "length")
motif.PA_15min <- data.frame(Amin.PA, Amax.PA, L.PA_bouts_15_uninterrupted)
colnames(motif.PA_15min) <- c("Amin", "Amax", "length")

motif.SB_10min <- data.frame(Amin.SB, Amax.SB, L.SB_bouts_10_uninterrupted)
colnames(motif.SB_10min) <- c("Amin", "Amax", "length")
motif.SB_20min <- data.frame(Amin.SB, Amax.SB, L.SB_bouts_20_uninterrupted)
colnames(motif.SB_20min) <- c("Amin", "Amax", "length")
motif.SB_30min <- data.frame(Amin.SB, Amax.SB, L.SB_bouts_30_uninterrupted)
colnames(motif.SB_30min) <- c("Amin", "Amax", "length")

motif.SB_30min.interrupted <- rbind(motif.SB_10min, motif.SB_10min, motif.SB_10min)

### 2) Apply forward algorithm for motif probability calculation

# For each participant load the fitted hsmm
filelist <- list.files(sequencedir, pattern = ".RData")

prob_PA_05min <- c()
prob_PA_5min <- c()
prob_PA_10min <- c()
prob_PA_15min <- c()

prob_SB_10min <- c()
prob_SB_20min <- c()
prob_SB_30min <- c()
prob_SB_30min.interrupted <- c()

id <- c()

for(pp in 1:length(filelist)){
  load(paste0(sequencedir, "/", filelist[pp]))
  id <- c(id, strsplit(filelist[pp], "_")[[1]][1])
  prob_PA_05min <- c(prob_PA_05min, motif_probability(motif.PA_0.5min, hsmms))
  prob_PA_5min <- c(prob_PA_5min, motif_probability(motif.PA_5min, hsmms))
  prob_PA_10min <- c(prob_PA_10min, motif_probability(motif.PA_10min, hsmms))
  prob_PA_15min <- c(prob_PA_15min, motif_probability(motif.PA_15min, hsmms))
  prob_SB_10min <- c(prob_SB_10min, motif_probability(motif.SB_10min, hsmms))
  prob_SB_20min <- c(prob_SB_20min, motif_probability(motif.SB_20min, hsmms))
  prob_SB_30min <- c(prob_SB_30min, motif_probability(motif.SB_30min, hsmms))
  prob_SB_30min.interrupted <- c(prob_SB_30min.interrupted, motif_probability(motif.SB_30min.interrupted, hsmms))
}

probabilities.heuristic <- as.data.frame(cbind(id, prob_PA_05min, prob_PA_5min, prob_PA_10min, prob_PA_15min, prob_SB_10min, prob_SB_20min, prob_SB_30min, prob_SB_30min.interrupted))
probabilities.heuristic <- cbind(probabilities.heuristic$id, dplyr::mutate_all(probabilities.heuristic[,-1], as.numeric))
colnames(probabilities.heuristic) <- c("id", "PA.05min", "PA.5min", "PA.10min", "PA.15min", "SB.10min", "SB.20min", "SB.30min", "SB.30min.interrupted")






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
cor.test(probabilities.heuristic$zBMI, exp(as.numeric(probabilities.heuristic$PA.1min)), na.action = "omit")
cor.test(probabilities.heuristic$zBMI, exp(as.numeric(probabilities.heuristic$PA.5min)), na.action = "omit")
cor.test(probabilities.heuristic$zBMI, exp(as.numeric(probabilities.heuristic$PA.10min)), na.action = "omit")
cor.test(probabilities.heuristic$zBMI, exp(as.numeric(probabilities.heuristic$PA.15min)), na.action = "omit")

# higher probability of longer uninterrupted sedentary behavior is positively related to BMI -- NOT TRUE
cor.test(probabilities.heuristic$zBMI, exp(as.numeric(probabilities.heuristic$SB.10min)), na.action = "omit")
cor.test(probabilities.heuristic$zBMI, exp(as.numeric(probabilities.heuristic$SB.20min)), na.action = "omit")
cor.test(probabilities.heuristic$zBMI, exp(as.numeric(probabilities.heuristic$SB.30min)), na.action = "omit")

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

