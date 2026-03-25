rm(list=ls())
library(readxl)
path_X <- "data/Thorax_X_101.xlsx"
path_Y <- "data/Thorax_Y_101.xlsx"
path_Z <- "data/Thorax_Z_101.xlsx"
SheetsNames <- excel_sheets(path_X)
X_list <- lapply(SheetsNames, function(sheet) read_excel(path_X, sheet = sheet))
Y_list <- lapply(SheetsNames, function(sheet) read_excel(path_Y, sheet = sheet))
Z_list <- lapply(SheetsNames, function(sheet) read_excel(path_Z, sheet = sheet))
names(X_list) <- SheetsNames
names(Y_list) <- SheetsNames
names(Z_list) <- SheetsNames

row_means_X_male <- X_list[["Male"]]
row_means_X_female <- X_list[["Female"]]

row_means_Y_male <- Y_list[["Male"]]
row_means_Y_female <- Y_list[["Female"]]

row_means_Z_male <- Z_list[["Male"]]
row_means_Z_female <- Z_list[["Female"]]




mean_male <- rowMeans(row_means_X_male, na.rm = TRUE)
mean_female <- rowMeans(row_means_X_female, na.rm = TRUE)

mean_male_Y <- rowMeans(row_means_Y_male, na.rm = TRUE)
mean_female_Y <- rowMeans(row_means_Y_female, na.rm = TRUE)

mean_male_Z <- rowMeans(row_means_Z_male, na.rm = TRUE)
mean_female_Z <- rowMeans(row_means_Z_female, na.rm = TRUE)



df_plot <- data.frame(
  time = seq_len(length(mean_male)),
  Male = mean_male,
  Female = mean_female
)

df_plot_Y <- data.frame(
  time = seq_len(length(mean_male_Y)),
  Male = mean_male_Y,
  Female = mean_female_Y
)

df_plot_Z <- data.frame(
  time = seq_len(length(mean_male_Z)),
  Male = mean_male_Z,
  Female = mean_female_Z
)

library(ggplot2)

ggplot(df_plot, aes(x = time)) +
  geom_line(aes(y = Male, color = "Male"), linewidth = 1) +
  geom_line(aes(y = Female, color = "Female"), linewidth = 1) +
  labs(
    title = "Mean Thorax X Over Time",
    x = "Time",
    y = "Mean Thorax X",
    color = "Sex"
  ) +
  theme_minimal() + ylim(-4, 16)

# Note: -Male and -Female
ggplot(df_plot_Y, aes(x = time)) +
  geom_line(aes(y = Male, color = "Male"), linewidth = 1) +
  geom_line(aes(y = Female, color = "Female"), linewidth = 1) +
  labs(
    title = "Mean Thorax Y Over Time",
    x = "Time",
    y = "Mean Thorax Y",
    color = "Sex"
  ) +
  theme_minimal() + ylim(-5, 4)


ggplot(df_plot_Z, aes(x = time)) +
  geom_line(aes(y = Male, color = "Male"), linewidth = 1) +
  geom_line(aes(y = Female, color = "Female"), linewidth = 1) +
  labs(
    title = "Mean Thorax Z Over Time",
    x = "Time",
    y = "Mean Thorax Z",
    color = "Sex"
  ) +
  theme_minimal() + ylim(-8, 6)



#####Attempt 1#####
#First attempt at tau_fun, using pointwise squared difference and mean of that as
#scale. Not great, but shows the idea of plotting the scaled difference and horizontal lines
#for different tau.
tau_fun <- function(y1, y2, taus = seq(0.1, 1, by = 0.1)) {
  
  
  
  dif <- (y1 - y2)^2
  mean_total <- mean(dif)
  
  scaled_dif <- dif / mean_total
  
  # Set up colors automatically
  colors <- rainbow(length(taus))
  
  # Plot main curve
  plot(0:(length(y1)-1), scaled_dif, type = "l",
       main = "Scaled Difference",
       ylab = "Scaled Difference",
       xlab = "x",
       ylim = c(0, max(taus , scaled_dif) * 1.1))
  
  # Add horizontal lines
  for (i in seq_along(taus)) {
    abline(h = taus[i] ,
           col = colors[i],
           lwd = 2,
           lty = 2)
  }
  
  # Add legend
  legend("topright",
         legend = paste("tau =", taus),
         col = colors,
         lty = 2,
         lwd = 2)
  
}


tau_fun(mean_male, mean_female)
tau_fun(mean_male_Y, mean_female_Y)
tau_fun(mean_male_Z, mean_female_Z)

zero_signal <- rep(0, 101)
zeroone_signal <- rep(0.1, 101)
tau_fun(zero_signal, zeroone_signal)

four_signal <- rep(4, 101)
fourone_signal <- rep(4.1, 101)
tau_fun(four_signal, fourone_signal) 

oh_signal <- rep(100, 101)
oh_p_one_signal <- rep(100.1, 101)
tau_fun(oh_signal, oh_p_one_signal) 

oh_signal <- rep(100, 101)
ohtw_signal <- rep(120, 101)
tau_fun(oh_signal, ohtw_signal) 

source("R/utilities.R")
# effect trajectory nonzero range percentages
cent_ranges <- centered_ranges(seq(5, 100, by = 5))
p_range <- which(cent_ranges[, "Percentage"] == 10)
signal1 <- square_pulse(start_end_pulse = c(cent_ranges$Start[p_range],
                                            cent_ranges$End[p_range]),
                                            start_height = 0, pulse_height = 1)
signal2 <- square_pulse(start_end_pulse = c(cent_ranges$Start[p_range],
                                            cent_ranges$End[p_range]),
                        start_height = 0.00001, pulse_height = 1)

tau_fun(signal1, zero_signal)
tau_fun(signal2, zero_signal)
######simulated signals#####


tau_fun <- function(y1, y2, taus = seq(0.1, 1, by = 0.1)) {
  
  rel_dif <- abs(y1 - y2) / sqrt(mean((y1^2 + y2^2)/2)) #best
  # rel_dif <- (y1 - y2)^2 / sqrt(mean((y1^2 + y2^2)/2)) #alt
  
  colors <- rainbow(length(taus))
  
  plot(rel_dif, type = "l",
       main = "Scale-normalized difference",
       ylab = "Relative difference",
       xlab = "x",
       ylim = c(0, max(taus , rel_dif) * 1.1) * 1.1)
  
  for (i in seq_along(taus)) {
    abline(h = taus[i],
           col = colors[i],
           lwd = 2,
           lty = 2)
  }
  
  legend("topright",
         legend = paste("tau =", taus),
         col = colors,
         lty = 2,
         lwd = 2)
}




tau_fun(mean_male, mean_female)
tau_fun(mean_male_Y, mean_female_Y)
tau_fun(mean_male_Z, mean_female_Z)

zero_signal <- rep(0, 101)
zeroone_signal <- rep(0.1, 101)
tau_fun(zero_signal, zeroone_signal) 

four_signal <- rep(4, 101)
fourone_signal <- rep(4.1, 101)
tau_fun(four_signal, fourone_signal) 

oh_signal <- rep(100, 101)
oh_p_one_signal <- rep(100.1, 101)
tau_fun(oh_signal, oh_p_one_signal) 

oh_signal <- rep(100, 101)
ohtw_signal <- rep(120, 101)
tau_fun(oh_signal, ohtw_signal) 

source("R/utilities.R")
# effect trajectory nonzero range percentages
cent_ranges <- centered_ranges(seq(5, 100, by = 5))
p_range <- which(cent_ranges[, "Percentage"] == 10)
signal1 <- square_pulse(start_end_pulse = c(cent_ranges$Start[p_range],
                                            cent_ranges$End[p_range]),
                        start_height = 0, pulse_height = 1)
signal2 <- square_pulse(start_end_pulse = c(cent_ranges$Start[p_range],
                                            cent_ranges$End[p_range]),
                        start_height = 0.00001, pulse_height = 1)

tau_fun(signal1, zero_signal)
tau_fun(signal2, zero_signal)



##########Effect size########
tau_fun_sd <- function(m1, m2, sd1, sd2, taus = seq(0.1, 1, by = 0.1)) {
  
  pooled_sd <- sqrt((sd1^2 + sd2^2) / 2)
  effect <- abs(m1 - m2) / pooled_sd
  
  colors <- rainbow(length(taus))
  
  plot(effect, type = "l",
       main = "Variance-adjusted functional effect size",
       ylab = "Absolute Cohen's d(x)",
       xlab = "x",
       ylim = c(0, max(taus , effect) * 1.1))
  
  for (i in seq_along(taus)) {
    abline(h = taus[i],
           col = colors[i],
           lwd = 2,
           lty = 2)
  }
  
  legend("topright",
         legend = paste("tau =", taus),
         col = colors,
         lty = 2,
         lwd = 2)
}

D1prime <- function(m1, m2, sd1, sd2, taus = seq(0.1, 1, by = 0.1)) {
  
  pooled_sd <- sqrt((sd1^2 + sd2^2) / 2)
  effect <- abs(m1 - m2) / pooled_sd
  
  x <- seq_along(effect)
  colors <- rainbow(length(taus))
  
  # Empty plot (no axes, no labels)
  plot(NA, NA,
       xlim = range(x),
       ylim = c(0.5, length(taus) + 0.5),
       xlab = "",
       ylab = "",
       bty = "n")
  
  # For each tau create horizontal segments
  for (i in seq_along(taus)) {
    
    tau <- taus[i]
    idx <- which(effect > tau)
    
    if (length(idx) > 0) {
      
      # Split into consecutive runs
      runs <- split(idx, cumsum(c(1, diff(idx) != 1)))
      
      for (r in runs) {
        segments(min(r), i, max(r), i,
                 col = colors[i],
                 lwd = 6)
      }
    }
  }
}

sd_male <- apply(row_means_X_male, 1, sd)
sd_female <- apply(row_means_X_female, 1, sd)

sd_male_Y <- apply(row_means_Y_male, 1, sd)
sd_female_Y <- apply(row_means_Y_female, 1, sd)

sd_male_Z <- apply(row_means_Z_male, 1, sd)
sd_female_Z <- apply(row_means_Z_female, 1, sd)

tau_fun_sd(mean_male, mean_female, sd1=sd_male , sd2=sd_female)
tau_fun_sd(mean_male_Y, mean_female_Y, sd1=sd_male_Y , sd2=sd_female_Y)
tau_fun_sd(mean_male_Z, mean_female_Z, sd1=sd_male_Z , sd2=sd_female_Z)

zero_signal <- rep(0, 101)
zeroone_signal <- rep(0.1, 101)
tau_fun_sd(zero_signal, zeroone_signal, sd1=rep(0.1,101), sd2=rep(0.1,101))
tau_fun_sd(zero_signal, zeroone_signal, sd1=rep(1,101), sd2=rep(1,101))

four_signal <- rep(4, 101)
fourone_signal <- rep(4.1, 101)
tau_fun_sd(four_signal, fourone_signal, sd1=rep(0.2,101), sd2=rep(0.2,101))
tau_fun_sd(four_signal, fourone_signal, sd1=rep(2,101), sd2=rep(2,101))


oh_signal <- rep(100, 101)
oh_p_one_signal <- rep(100.1, 101)
tau_fun_sd(oh_signal, oh_p_one_signal, sd1=rep(10,101), sd2=rep(10,101)) 
tau_fun_sd(oh_signal, oh_p_one_signal, sd1=rep(0.2,101), sd2=rep(0.2,101))

oh_signal <- rep(100, 101)
ohtw_signal <- rep(120, 101)
tau_fun_sd(oh_signal, ohtw_signal, sd1=rep(2,101), sd2=rep(2,101)) 
tau_fun_sd(oh_signal, ohtw_signal, sd1=rep(20,101), sd2=rep(20,101))
tau_fun_sd(oh_signal, ohtw_signal, sd1=rep(40,101), sd2=rep(40,101)) 

source("R/utilities.R")
# effect trajectory nonzero range percentages
cent_ranges <- centered_ranges(seq(5, 100, by = 5))
p_range <- which(cent_ranges[, "Percentage"] == 10)
signal1 <- square_pulse(start_end_pulse = c(cent_ranges$Start[p_range],
                                            cent_ranges$End[p_range]),
                        start_height = 0, pulse_height = 1)
signal2 <- square_pulse(start_end_pulse = c(cent_ranges$Start[p_range],
                                            cent_ranges$End[p_range]),
                        start_height = 0.00001, pulse_height = 1)

tau_fun_sd(signal1, zero_signal, sd1=rep(1,101), sd2=rep(1,101))
tau_fun_sd(signal1, zero_signal, sd1=rep(0.5,101), sd2=rep(0.5,101))
tau_fun_sd(signal1, zero_signal, sd1=rep(4,101), sd2=rep(4,101))

tau_fun_sd(signal2, zero_signal, sd1=rep(1,101), sd2=rep(1,101))
tau_fun_sd(signal2, zero_signal, sd1=rep(0.5,101), sd2=rep(0.5,101))
tau_fun_sd(signal2, zero_signal, sd1=rep(4,101), sd2=rep(4,101))

