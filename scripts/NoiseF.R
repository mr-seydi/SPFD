set.seed(123)
library(tidyverse)
library(MASS)
library(ggplot2)
# Vectorized squared-exponential covariance
sq_exp_cov <- function(x, lambda, alpha, nugget = 1e-8) {
  # distance matrix (1D)
  d <- abs(outer(x, x, "-"))
  if (lambda == 0) {
    return(alpha^2 * diag(length(x)) + diag(nugget, length(x)))
  }
  K <- alpha^2 * exp(-(d^2) / (2 * lambda^2))
  K + diag(nugget, length(x))
}

# Corrected exponential covariance (vectorized)
exp_cov <- function(x, lambda, alpha, nugget = 1e-8) {
  d <- abs(outer(x, x, "-"))
  if (lambda == 0) {
    return(alpha^2 * diag(length(x)) + diag(nugget, length(x)))
  }
  K <- alpha^2 * exp(-d / lambda)
  K + diag(nugget, length(x))
}



GP <- function(x, alpha, lambda, n_samples = 10, cov_type = c("exp", "sqexp")) {
  cov_type <- match.arg(cov_type)
  if (cov_type == "sqexp") {
    K <- sq_exp_cov(x, lambda, alpha)
  } else {
    K <- exp_cov(x, lambda, alpha)
  }
  
  mat <- MASS::mvrnorm(
    n = n_samples,
    mu = rep(0, length(x)),
    Sigma = K
  )
  
  # Force mat to be a matrix even if n_samples = 1
  if(n_samples == 1) mat <- matrix(mat, nrow = 1)
  # Assign names (X1, X2, ...) to the rows of mat (which become columns after t())
  rownames(mat) <- paste0("X", seq_len(nrow(mat)))
  
  # transpose so columns are samples, then tidy
  t(mat) |>
    data.frame() |>
    mutate(alpha = alpha, lambda = lambda, x = x) |>
    pivot_longer(c(-x, -alpha, -lambda), names_to = "sample", values_to = "y")
}


grd <- expand.grid(alpha = c(0.5), lambda = c(0,0.01,0.05,0.1, 0.2, 0.5))
samples <- map2(grd$alpha, grd$lambda, ~ GP(seq(0, 1, length.out = 800), .x, .y, n_samples = 10, cov_type = "sqexp")) |> bind_rows()

ggplot(samples, aes(x = x, y = y, group = sample, colour = sample)) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~ alpha + lambda) +
  theme_minimal()

#######################################################

grd <- expand.grid(alpha = c(0.5), lambda = c(0,0.01,0.05,0.1, 0.2, 0.5))
samples <- map2(grd$alpha, grd$lambda, ~ GP(seq(0, 1, length.out = 101), .x, .y, n_samples = 10, cov_type = "sqexp")) |> bind_rows()

ggplot(samples, aes(x = x, y = y, group = sample, colour = sample)) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~ alpha + lambda) +
  theme_minimal()

######################################################
S <- 10000
SD_0.05 <-  matrix(0, nrow=101, ncol=S)
SD_0.25 <-  matrix(0, nrow=101, ncol=S)
SD_0.50 <-  matrix(0, nrow=101, ncol=S)
for (i in 1:S){
  samples_0.05 <- GP(seq(0, 1, length.out = 101), alpha = 0.5, lambda = 0.05, n_samples = 10, cov_type = "sqexp")
  samples_0.25 <- GP(seq(0, 1, length.out = 101), alpha = 0.5, lambda = 0.25, n_samples = 10, cov_type = "sqexp")
  samples_0.50 <- GP(seq(0, 1, length.out = 101), alpha = 0.5, lambda = 0.50, n_samples = 10, cov_type = "sqexp")
  
  sd_0.05 <- samples_0.05 |>
    group_by(x) |>
    summarise(sd_y = sd(y), .groups = "drop") |>
    pull(sd_y)
  
  sd_0.25 <- samples_0.25 |>
    group_by(x) |>
    summarise(sd_y = sd(y), .groups = "drop") |>
    pull(sd_y)
  
  sd_0.50 <- samples_0.50 |>
    group_by(x) |>
    summarise(sd_y = sd(y), .groups = "drop") |>
    pull(sd_y)
  
  SD_0.05[,i] <- sd_0.05
  SD_0.25[,i] <- sd_0.25
  SD_0.50[,i] <- sd_0.50
}


library(tidyverse)

sd_to_long <- function(SD) {
  as_tibble(SD) |>
    mutate(x_id = row_number()) |>
    pivot_longer(-x_id, values_to = "sd")
}

df_005 <- sd_to_long(SD_0.05)
df_025 <- sd_to_long(SD_0.25)
df_050 <- sd_to_long(SD_0.50)


# plot the point-wise boxplots of SDs, each column is a different sample
df_all <- bind_rows(
  mutate(df_005, lambda = "0.05"),
  mutate(df_025, lambda = "0.25"),
  mutate(df_050, lambda = "0.50")
)

plot1 <- ggplot(df_all, aes(x = factor(x_id), y = sd, fill = lambda)) +
  geom_boxplot(outlier.size = 0.2) +
  facet_wrap(~ lambda, ncol = 1) +
  labs(x = "x index", y = "SD", fill = "Lambda") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_fill_brewer(palette = "Set2")



####################################
source("R/utilities.R")
#Noise_generator (Sample_size, Continuum_size, Noise_mu, Noise_sig, Noise_fwhm)

S <- 10000
SD_0.05A <-  matrix(0, nrow=101, ncol=S)
SD_0.25A <-  matrix(0, nrow=101, ncol=S)
SD_0.50A <-  matrix(0, nrow=101, ncol=S)
for (i in 1:S){
  GP(seq(0, 1, length.out = 101), alpha = 0.5, lambda = 0.05, n_samples = 10, cov_type = "sqexp")
  SD_0.05A[,i] <- Noise_generator(Sample_size = 10, Continuum_size = 101, Noise_mu = 0, Noise_sig = 0.5, Noise_fwhm = 5)$SD
  SD_0.25A[,i] <- Noise_generator(Sample_size = 10, Continuum_size = 101, Noise_mu = 0, Noise_sig = 0.5, Noise_fwhm = 25)$SD
  SD_0.50A[,i] <- Noise_generator(Sample_size = 10, Continuum_size = 101, Noise_mu = 0, Noise_sig = 0.5, Noise_fwhm = 50)$SD
}


df_005 <- sd_to_long(SD_0.05A)
df_025 <- sd_to_long(SD_0.25A)
df_050 <- sd_to_long(SD_0.50A)


# plot the point-wise boxplots of SDs, each column is a different sample
df_all <- bind_rows(
  mutate(df_005, lambda = "0.05"),
  mutate(df_025, lambda = "0.25"),
  mutate(df_050, lambda = "0.50")
)

plot2 <- ggplot(df_all, aes(x = factor(x_id), y = sd, fill = lambda)) +
  geom_boxplot(outlier.size = 0.2) +
  facet_wrap(~ lambda, ncol = 1) +
  labs(x = "x index", y = "SD", fill = "Lambda") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_fill_brewer(palette = "Set2")

mean(SD_0.05)
mean(SD_0.25)
mean(SD_0.50)

mean(SD_0.05A)
mean(SD_0.25A)
mean(SD_0.50A)


#two column plot1 and plot2
library(gridExtra)
grid.arrange(plot1, plot2, ncol = 2)










samples |>
  group_by(alpha, lambda, x) |>
  summarise(sd_y = sd(y), .groups = "drop")

library(ggplot2)

# 1. Calculate the pointwise SD
sd_summary <- samples |>
  group_by(alpha, lambda, x) |>
  summarise(sd_y = sd(y), .groups = "drop")

# 2. Plot the results
ggplot(sd_summary, aes(x = x, y = sd_y)) +
  geom_line(color = "steelblue", linewidth = 1) +
  facet_grid(alpha ~ lambda, labeller = label_both) +
  labs(
    title = "Pointwise Standard Deviation of Gaussian Process Samples",
    subtitle = "Grouped by Alpha (rows) and Lambda (columns)",
    x = "Input Variable (x)",
    y = "Standard Deviation (σ)"
  ) +
  theme_minimal()

