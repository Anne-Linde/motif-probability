hsmmDir <- "/Users/annelindelettink/GECKO/preprocessing/manuscript/subset/output_rawinput_5years_complete_anthro/epochdata/validdata/mhsmmdata/models"
load(paste0("/Users/annelindelettink/GECKO/preprocessing/manuscript/subset/analyses", "/merged_data_rerun.RData"))

## Calculate descriptive statistics
merged_data$zbmi <- as.numeric(merged_data$zbmi)
merged_data$measure_age <- as.numeric(merged_data$measure_age)

funs <- list(sd=sd, mean=mean)
descriptives <- as.data.frame(sapply(funs, function(x) sapply(merged_data, x, na.rm=T)))
table(merged_data$sex) # 1 = boys, 0 = girls


# Calculate average state means, sd, and durations
files <- list.files(hsmmDir)

means_state1 <- c()
vars_state1 <- c()
means_state2 <- c()
vars_state2 <- c()
means_state3 <- c()
vars_state3 <- c()
lambda_state1 <- c()
lambda_state2 <- c()
lambda_state3 <- c()

for(file in 1:length(files)){
  load(paste0(hsmmDir, "/", files[file]))
  means <- sort(hsmms$model$parms.emission$mu, decreasing = FALSE)
  means_state1 <- c(means_state1, means[1])
  means_state2 <- c(means_state2, means[2])
  means_state3 <- c(means_state3, means[3])
  index = which(means == hsmms$model$parms.emission$mu)
  vars <- hsmms$model$parms.emission$sigma
  
  vars_state1 <- c(vars_state1, vars[index[1]])
  vars_state2 <- c(vars_state2, vars[index[2]])
  vars_state3 <- c(vars_state3, vars[index[3]])
  lambdas <- hsmms$model$sojourn$lambda / (60/15)
  lambda_state1 <- c(lambda_state1, lambdas[index[1]])
  lambda_state2 <- c(lambda_state2, lambdas[index[2]])
  lambda_state3 <- c(lambda_state3, lambdas[index[3]])
}

state_df <- as.data.frame(cbind(means_state1, vars_state1, means_state2, vars_state2, means_state3, vars_state3,
                  lambda_state1, lambda_state2, lambda_state3))
summary(state_df)
## Figure 1: graph (rest of the figure was created in Canva)
# Plot and save images final hsmm and the state sojourn
library(lubridate)
pp = 15
load(paste0(hsmmDir, "/",list.files(hsmmDir)[pp]))
load(paste0("/Users/annelindelettink/GECKO/preprocessing/manuscript/subset/output_rawinput_5years_complete_anthro/epochdata/validdata/", "/",list.files("/Users/annelindelettink/GECKO/preprocessing/manuscript/subset/output_rawinput_5years_complete_anthro/epochdata/validdata/")[pp]))
load(paste0("/Users/annelindelettink/GECKO/preprocessing/manuscript/subset/output_rawinput_5years_complete_anthro/epochdata/validdata/mhsmmdata", "/",list.files("/Users/annelindelettink/GECKO/preprocessing/manuscript/subset/output_rawinput_5years_complete_anthro/epochdata/validdata/")[pp]))

palette <- viridis::viridis(length(unique(data_plot$State)), option = "D")

data_plot <- data.frame(
  Time = validData[[1]]$timestamp,
  Observed = mhsmmdata$Y[1:mhsmmdata$N[1]],
  State = as.factor(hsmms$yhat[1:mhsmmdata$N[1]]))

data_plot$Time <- as.POSIXct(data_plot$Time, format = "%Y-%m-%dT%H:%M:%S", tz = "Europe/Amsterdam")

# Initialize data for the rectangles
categories <- unique(data_plot$State)

xmin <- c()
xmax <- c()
ymin <- c()
ymax <- c()
col = c()

for (t in 1:(length(data_plot$Time)-1)) {
  xmin <- c(xmin, as.POSIXct(data_plot$Time[t], origin = "1970-01-01", tz = "Europe/Amsterdam"))
  xmax <- c(xmax, as.POSIXct(data_plot$Time[t + 1], origin = "1970-01-01", tz = "Europe/Amsterdam"))
  ymin = c(ymin, 0)
  ymax = c(ymax, max(data_plot$Observed, na.rm = TRUE))
  col <- c(col, palette[which(categories == data_plot$State[t])])
}


df_rect <- as.data.frame(list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, col = col))
df_rect$xmin <- as.POSIXct(df_rect$xmin, origin = "1970-01-01", tz = "Europe/Amsterdam")
df_rect$xmax <- as.POSIXct(df_rect$xmax, origin = "1970-01-01", tz = "Europe/Amsterdam")


# Create the plot
g <- ggplot() +
  # Add state rectangles
  geom_rect(
    data = df_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = col),
    #alpha = 0.35
  ) +
  # Add hip data line
  geom_line(data = data_plot, aes(x = Time, y = Observed), color = "black") +
  theme_classic() +
  scale_fill_identity(name = "States", labels = categories, guide = "legend") +
  scale_y_continuous(name = "Acceleration") +
  xlab("Time (hh:mm)")






# Define the states' mu (mean) and sigma (variance)
mu <- hsmms$model$parms.emission$mu #means
sigma_sd <- sqrt(hsmms$model$parms.emission$sigma)  # sd

# Create a sequence of values for the x-axis (observed data range)
x_vals <- seq(min(mu - 3*sigma_sd), max(mu + 3*sigma_sd), length.out = 100)

# Plot the Gaussian distributions for each state
library(ggplot2)

plot_data <- data.frame(x = rep(x_vals, 3),
                        y = c(dnorm(x_vals, mean = mu[1], sd = sigma_sd[1]),
                              dnorm(x_vals, mean = mu[2], sd = sigma_sd[2]),
                              dnorm(x_vals, mean = mu[3], sd = sigma_sd[3])),
                        State = rep(c(1, 2, 3), each = length(x_vals)))

ggplot(plot_data, aes(x = x, y = y, color = factor(State))) +
  scale_color_manual(values = palette) +  # Map states to your custom palette
  
  geom_line() +
  labs(title = "Emission Distributions for Each State",
       x = "Observed Value",
       y = "Density",
       color = "State") +
  theme_minimal()



# Example: Poisson durations modeled by exponential distributions
# lambda is the rate parameter, i.e., the expected number of events per unit time.

# Define lambda for each state
lambda <-  hsmms$model$sojourn$lambda + 1

# Create a sequence of possible durations (x-values)
x_vals <- seq(0, 1, length.out = 100)

# Create the probability density for each state's exponential distribution
plot_data <- data.frame(x = rep(x_vals, 3),
                        y = c(dexp(x_vals, rate = lambda[1]),
                              dexp(x_vals, rate = lambda[2]),
                              dexp(x_vals, rate = lambda[3])),
                        State = rep(c(1, 2, 3), each = length(x_vals)))

# Create the plot using ggplot2
library(ggplot2)

ggplot(plot_data, aes(x = x, y = y, color = factor(State))) +
  geom_line(size = 1.2) +  # Line thickness
  scale_color_manual(values = palette) +  # Map states to your custom palette
  labs(title = "Exponential Distributions for Poisson Durations",
       x = "Duration (Time between events)",
       y = "Density",
       color = "State") +
  theme_minimal()

