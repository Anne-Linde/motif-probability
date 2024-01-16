## TO DO: check welke funcites nog van belang zijn/ worden gebruikt

# Function to randomly selects a fraction from the data to train the model (trainSet) and to fit the model on the data (validationSet)
validationSetApproach <- function(data, trainingsFract) {
  numberPP <- nrow(data)
  trainIndex <- sample(1:numberPP, trainingsFract * numberPP)
  trainSet <- data[trainIndex, ]
  validationSet <- data[-trainIndex, ]
  return(list(trainSet = trainSet, validationSet = validationSet))
}

# Function to put the data into a hsmm.data object (required format for training hsmm with the mhsmm package) incorporating non-wear time
formatMhsmm <- function(data) {
  nSequences <- nrow(data)
  nObservations <- rep(ncol(data), nSequences)
  object <- list(x = unlist(data), N = nObservations)
  class(object) <- "hsmm.data"
  return(object)
}

# Function to train the HSMMs with the varying number of states, given model specifications. The function saves all trained models
train.mhsmm <- function(trainData, minStates, maxStates, maxIter, maxDur, bouts) {
  train <- formatMhsmm(trainData) # Put the training data into a hsmm.data object
  #train <- formatMhsmmOmitNA(trainData) # Put the training data into a hsmm.data object
  states <- minStates : maxStates # Vector of the possible number of states
  # Train HSMMs with varying number of states
  bics <- c()
  hsmms <- list()
  for(n in 1:length(states)) { # Loop over each number of states in the state vector
    # Initial model parameters
    init <- rep(1 / states[n], states[n]) # State probabilities
    tm0 <- matrix(runif(states[n] ^ 2), nrow = states[n]) # Transition matrix
    diag(tm0) <- rep(0, states[n]) # Diagonal of the transition matrix must be 0
    for (j in 1:states[n]) {
      if (minStates > 2 & n != 1)
        tm0[j, ] <- tm0[j, ] / sum(tm0[j, ])
    }
    if (any(is.na(train$x)))
      m0 <- seq(quantile(train$x, 0.05, na.rm = TRUE), 
                quantile(train$x, 0.95, na.rm = TRUE), length.out = states[n])
    else
      m0 <- seq(quantile(train$x, 0.05), quantile(train$x, 0.95), 
                length.out = states[n])
    odpar0 <- list(mu = m0, sigma = 10 / (1 : states[n] + 10)) # Observation distribution (emission matrix): Gaussian
    soj_dur <- bouts * 60 / epoch # Number of epochs for bout durations
    soj0 <- list(lambda = seq(soj_dur[1], soj_dur[4], length.out = states[n]), 
                 shift = rep(1, states[n]), type = "poisson") # Run length distribution: Poisson
    
    # Model specification
    startmodel <- hsmmspec(init, tm0, odpar0, soj0, dens.emis = dnorm.hsmm)
    
    # EM algorithm to fit a HSMM to the training data
    hsmm <- hsmmfit(train, startmodel, mstep = mstep.norm, M = dMax, 
                    maxit = qMax)
    #assign(paste("model", toString(states[n]), sep = "_"), hsmm) # Save the fitted model
    hsmms[[n]] <- hsmm
    
    # Decide if early stopping is necessary is BIC of current model lower than BIC of previous model
    bic <- extractBIC.mhsmm(hsmm) # Extract BIC of current model % increase of 2 N^D_(js)
    if (n == 1)
      bics <- c(bics, bic)
    else {
      if ((abs(bic - bics[n - 1])/bic) * 100 > 5) # If more complex model does not change much (5%) stop training the model
        bics <- c(bics, bic)
      else
        break
    }
  }
  return(hsmms)
}

# Decide on optimal number of states based on BIC criterium and return the optimal model
nStatesOptimal <- function(hsmmList) {
  BICs <- c()
  for (s in 1:length(hsmmList)) {
    model <- hsmmList[[s]] # Get fitted model information
    bic <- extractBIC.mhsmm(model) # Calculate BIC
    BICs <- c(BICs, bic)
  }
  finalModel <- hsmmList[[which.min(BICs)]] # Position of min BIC value in vector is the position of the model with the optimal number of states
  return(finalModel)
}

extractBIC.mhsmm <- function(hsmm) {
  maxLogl <- max(hsmm$loglik, na.rm = TRUE) # Get maximum loglikelihood value of the model's iterations
  nParameters <- length(unlist(hsmm$model)) - sum(length(unlist(hsmm$model$d)), 
                                                  length(unlist(hsmm$model$D)), length(unlist(hsmm$model$dens.emission))) # Get number model parameters
  nObservations <- hsmm$estep_variables$nsequences # Get number observations
  bic <- (nParameters * log(nObservations)) - (2 * log(maxLogl)) # Calculate BIC
  return(bic)
}

# Function to train the HSMM given model specifications of the optimal model. The function returns the most likely sequences given the observed (validation) data
fit.mhsmm <- function(validationData, finalModel) {
  # Optimal model parameter specification
  sp <- finalModel$model$init
  tm <- finalModel$model$transition
  odpar <- list(mu = finalModel$model$parms.emission$mu, sigma = 
                  finalModel$model$parms.emission$sigma)
  soj <- list(lambda = finalModel$model$sojourn$lambda, shift = 
                finalModel$model$sojourn$shift, type = finalModel$model$sojourn$type)
  
  # Model specification
  model <- hsmmspec(sp, tm, odpar, soj, dens.emis = dnorm.hsmm)
  
  # Information aquisition
  validation <- formatMhsmm(validationData) # Put the validation data into a hsmm.data object
  yhat <- predict(model, validation, method = "viterbi") # Infer hidden states: Determine the most likely sequence given the observed (validation) data
  return(yhat)
}

### Functions for plotting ###
# Function that returns a graph object of the hsmm
visualizeHSMM   <- function(hsmmModel, colPalette, save = FALSE, path = NULL) {
  require(igraph)
  if(save)
    png(file = paste(path, "hsmm_state_information.png"))
  
  g <- make_ring(hsmmModel$J, directed = TRUE, mutual = TRUE, circular = TRUE) # Make a directed circle graph with the optimal number of states
  # Vertex specifications
  vertex_attr(g, "label") <- paste(rep("state", finalModel$J), 
                                   c(1:finalModel$J), sep = " ")
  vertex_attr(g, "color") <- colPalette
  vertex_attr(g, "size") <- 75
  
  # Edge specifications
  edges <- as.data.frame(get.edgelist(g)) # Put the edges of the graph in a data.frame object
  edgeTrans <- c()
  edgeCol <- c()
  for(e in 1:length(E(g))) {
    edgeTrans <- c(edgeTrans, hsmmModel$model$transition[edges$V1[e], 
                                                         edges$V2[e]]) # Get the transitions corresponding to the graph edges and put in a vector
    edgeCol <- c(edgeCol, edges$V2[e]) # Set the color the edge to the vertex (outgoing)
  }
  edge_attr(g, "label") <- round(edgeTrans, digits = 4) # Set the edge labels to the corresponding transition value, rounded to 4 digits
  edge.start <- ends(g, es = E(g), names = F)[, 1] # Set start point of edges to get the correct color
  edge_attr(g, "color") <- V(g)$color[edge.start] # Set edge color attributes
  edgeCol <- rgb(t(col2rgb(V(g)$color[edge.start]) / 255), alpha = .7) # Set edge label color attributes somewhat lighter for readability
  plot(g, edge.curved = 0.4, edge.arrow.size = .7, edge.label.color = edgeCol)
  if(save) {
    while(!is.null(dev.list()))
      dev.off()
  }
}

# Function to plot the sojourn densities
plotSoj <- function(model, color, save = FALSE, path = NULL) {
  tmp = model$model$d
  if(save)
    png(file = paste(path, "hsmm_sojourn.png"))
  plot(1:nrow(tmp), tmp[, 1], type = "l", ylab = "d(u)", xlab = "u", 
       ylim = range(tmp), col = color, lwd = 2)
  for (i in 2:model$J) {
    lines(tmp[, i], type = "l", col = color[i], lwd = 2)
  }
  legend("topright", legend = paste("state" , c(1:model$J)), 
         col = color, lty = 1, lwd = 2)
  if(save) {
    while(!is.null(dev.list()))
      dev.off()
  }
}

# Function to put the predicted data into a data.frame object (required format for plotting accelerometer data with state information)
toDF <- function(hsmmdataPredicted, validationSet) {
  pp <- c() # Make a vector to link the data with the subject
  time <- c() # Make a vector to link the time steps
  for (s in 1:length(hsmmdataPredicted$N)) {
    pp <- c(pp, rep(rownames(validationSet)[s], hsmmdataPredicted$N[s]))
    time  <- c(time, c(1:hsmmdataPredicted$x$N[s]))
  }
  # Make data.frame object of the data required to plot
  data <- cbind(time, pp, hsmmdataPredicted$s, hsmmdataPredicted$x$x)
  colnames(data) <- c("time", "pp", "state", "acc")
  data <- as.data.frame(data)
  return(data)
}

# Function to save or plot the accelerometer data for each subject separately with or without states (requires data.frame object as input); adapted from J. O'Connell plot.hsmm.data() mhsmm package
plot.hsmm.data2 <- function(df, col, addStates = TRUE, save = FALSE, path = NULL) {
  nPp <- as.numeric(as.vector(unique(df$pp)))
  for (s in 1:length(nPp)) { # Plot for each subject
    pp <- nPp[s] # Current subject
    if(save) {
      filepath <- paste(paste("states_pp", toString(pp), sep = ""), ".png", 
                        sep = "")
      png(file = paste(path, filepath))
    }
    plot(ts(as.vector(df$acc[df$pp == pp])), ylab = paste("pp", pp, 
                                                          sep = " ")) #The time series
    if (addStates) { # Add the states
      states <- as.vector(as.integer(df$state[df$pp == pp]))
      addStates2(states, col)
    }
    if(save) {
      while(!is.null(dev.list()))
        dev.off()
    }
  }
}

# Function to plot state heatmaps for each subject (requires data.frame object as input)
plotHeat <- function(df, color, save = FALSE, path = NULL) {
  require(ggplot2)
  if(save) {
    filepath <- paste("heatmap_states", ".png", sep = "")
    png(file = paste(path, filepath))
  }
  ggplot(df, aes(x = time, y = pp, fill = as.factor(state))) + geom_tile() +
    guides(col = guide_legend(title = "States")) + ylab("subject") +
    scale_fill_manual(values = color) +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = 
                         element_blank())
  if(save) {
    while(!is.null(dev.list()))
      dev.off()
  }
}

# Function to add the states to the time series plot adopted from J. O'Connell (only the color adapted!)
addStates2 <- function (states, cols, x = NULL, ybot = axTicks(2)[1], 
                        ytop = ybot + (axTicks(2)[2] - axTicks(2)[1]) / 5, dy  = ytop - ybot, 
                        greyscale = FALSE, leg = NA, J = length(unique(states)), time.scale = 1, 
                        shiftx = 0) {
  draw.it  <- function(hats, ybot, ytop, cols, greyscale) {
    ##cat("ybot", ybot, "ytop", ytop, "\n")
    for (ii in 1:length(hats$state)) {
      if (greyscale) {
        rect(xleft = hats$intervals[ii], ybottom = ybot,
             xright = hats$intervals[ii + 1], ytop = ytop,
             col = cols[hats$state[ii]], border = 1)
      }
      else {
        rect(xleft = hats$intervals[ii], ybottom = ybot,
             xright = hats$intervals[ii + 1], ytop = ytop,
             col = cols[hats$state[ii]], border = cols[hats$state[ii]])
      }
    }
  }
  if (is.null(states)) {
    states <- x
    if (!is.list(states))
      states <- list(states)
    x <- seq_along(states[[1]])
  }
  else {
    if (!is.list(states))
      states <- list(states)
    if (is.null(x))
      x <- seq_along(states[[1]])
  }
  x <- as.numeric(x)
  rr <- range(x)
  J = length(unique(states))
  if (greyscale)
    cols <- c("#FFFFFF", "#F0F0F0", "#D9D9D9", "#BDBDBD", "#969696", 
              "#737373", "#525252", "#252525")
  else
    cols <- rgb(t(col2rgb(col)/255), alpha = .95) # Color adapted here! (AL)
  st.list <- lapply (states, function(st) {
    runs = rle(st)
    cs  <- cumsum(c(0, runs$lengths))
    hats <- list(intervals=rr[1]+ diff(rr)*cs/max(cs), states=runs$values)
    hats
  })
  for (ii in seq_along(st.list)) {
    draw.it (st.list[[ii]], ybot, ytop, cols, greyscale)
    ybot <- ytop + .2 * dy
    ytop <- ybot + dy
  }
  if (any(!is.na(leg)))
    legend("topleft", legend = leg, fill = cols, bg = "white")
}
