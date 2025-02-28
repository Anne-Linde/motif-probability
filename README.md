# motif-probability
The code in this repository corresponds to our motif probability methodology for accelerometer data.

The probability of a specific physical behavior pattern in the data can be calculated using the forward algorithm of the hidden semi-Markov model (HSMM). This algorithm requires 2 inputs:

1. User-defined motif(s)
The motif probability is defined as the probability of a specific sequence of bouts, each characterized by acceleration ranges = [Amin, Amax] and their corresponding duration in minutes (L). 

For example the hypothetical motif:

| Sequence of bouts | Amin | Amax | duration |
|----------|----------|----------|----------|
| Bout 1    | 3   | 4  | 30   |
| Bout 2    | 0   | 1.5   | 30 |

This hypothetical motif can be defined as follows using the function defineMotif.R:
motif <- defineMotif(Amin = c(3, 0), Amax = c(4, 1.5), duration = 30, nEpochsMin = 60/15)

Note that the nEpochsMin is required to determine the number of epochs corresponding to the duration and epoch length of the accelerometer data used to fit the hsmm

2. Parameters of a fitted HSMM
For data reduction purposes an HSMM is fitted on the high dimensional accelerometer data using the mhsmm package.
The parameters of this fitted HSMM (transition matrix, and state parameters (i.e., Gaussian acceleration distribution and Poisson duration distribution parameters)) then lend themselves for the calculation of this motif probability.

**Motif probability calculation**
The motifProbability function calculates the probabilities of the user-defined motifs for each participant using the estimated HSMM parameters. The process is as follows:

1. Calculation of Bout Probabilities:
For each bout in the motif and each state: 
  a. Acceleration probability: the probability of observing the acceleration values given the underlying state, using the acceleration distribution. 
  b. Duration probability: the probability of observing the bout duration given the underlying state, using the, using the duration distribution. 
  c. Accounting for the transition probability: probability of transitioning from the previous state to the current state. 
  d. Sum these probabilities over all possible states to get the forward probability for each bout in the motif.
2. Combine the Bout Probabilities

Note that the function calculateMotifProbabilities function can be used for the calculation of multiple motifs.
