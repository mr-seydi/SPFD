rm(list = ls())
#rows are sample sizes, columns are MC iterations
n_samples <- 5:50

X_ttest <- readRDS(paste0("outputs/TestResults/X_ttest", ".rds"))
Y_ttest <- readRDS(paste0("outputs/TestResults/Y_ttest", ".rds"))
Z_ttest <- readRDS(paste0("outputs/TestResults/Z_ttest", ".rds"))

X_wilcox <- readRDS(paste0("outputs/TestResults/X_wilcox", ".rds"))
Y_wilcox <- readRDS(paste0("outputs/TestResults/Y_wilcox", ".rds"))
Z_wilcox <- readRDS(paste0("outputs/TestResults/Z_wilcox", ".rds"))

# Calculate power for each sample size and test
power_ttest_X <- rowMeans(X_ttest < 0.05)
power_ttest_Y <- rowMeans(Y_ttest < 0.05)
power_ttest_Z <- rowMeans(Z_ttest < 0.05)

power_wilcox_X <- rowMeans(X_wilcox < 0.05)
power_wilcox_Y <- rowMeans(Y_wilcox < 0.05)
power_wilcox_Z <- rowMeans(Z_wilcox < 0.05)

# Create a data frame for plotting
# power_df <- data.frame(
#   SampleSize = rep(n_samples, 6),
#   Power = c(power_ttest_X, power_ttest_Y, power_ttest_Z,
#             power_wilcox_X, power_wilcox_Y, power_wilcox_Z),
#   Test = rep(c("t-test", "Wilcoxon"), each = length(n_samples
#   ) * 3),
#   Axis = rep(c("X", "Y", "Z"), times = 2 *
#   length(n_samples))
# )


library(readxl)
source("R/utilities.R")

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

n_X_male <- ncol(row_means_X_male)
n_X_female <- ncol(row_means_X_female)

n_Y_male <- ncol(row_means_Y_male)
n_Y_female <- ncol(row_means_Y_female)

n_Z_male <- ncol(row_means_Z_male)
n_Z_female <- ncol(row_means_Z_female)



mean_male_X <- rowMeans(row_means_X_male, na.rm = TRUE)
mean_female_X <- rowMeans(row_means_X_female, na.rm = TRUE)

mean_male_Y <- rowMeans(row_means_Y_male, na.rm = TRUE)
mean_female_Y <- rowMeans(row_means_Y_female, na.rm = TRUE)

mean_male_Z <- rowMeans(row_means_Z_male, na.rm = TRUE)
mean_female_Z <- rowMeans(row_means_Z_female, na.rm = TRUE)


sd_male_X <- apply(row_means_X_male, 1, sd)
sd_female_X <- apply(row_means_X_female, 1, sd)

sd_male_Y <- apply(row_means_Y_male, 1, sd)
sd_female_Y <- apply(row_means_Y_female, 1, sd)

sd_male_Z <- apply(row_means_Z_male, 1, sd)
sd_female_Z <- apply(row_means_Z_female, 1, sd)



D1_X <- which(mean_male_X != mean_female_X)
D1_Y <- which(mean_male_Y != mean_female_Y)
D1_Z <- which(mean_male_Z != mean_female_Z)

taus = seq(0.1, 1, by = 0.1)
D1p_X <- compute_D1p(mean_male_X, mean_female_X, sd_male_X, sd_female_X,
                     n_X_male, n_X_female, taus)
D1p_X_ind <- make_indices(D1p_X)
D1p_Y <- compute_D1p(mean_male_Y, mean_female_Y, sd_male_Y, sd_female_Y,
                     n_Y_male, n_Y_female, taus)
D1p_Y_ind <- make_indices(D1p_Y)
D1p_Z <- compute_D1p(mean_male_Z, mean_female_Z, sd_male_Z, sd_female_Z,
                     n_Z_male, n_Z_female, taus)
D1p_Z_ind <- make_indices(D1p_Z)


results <- list()

for (n in n_samples) {
  
  cat("Processing sample size:", n, "\n")
  
  X_TWT <- readRDS(paste0("outputs/TestResults/X_TWT_n_", n, ".rds"))
  Y_TWT <- readRDS(paste0("outputs/TestResults/Y_TWT_n_", n, ".rds"))
  Z_TWT <- readRDS(paste0("outputs/TestResults/Z_TWT_n_", n, ".rds"))
  
  X_IATSE <- readRDS(paste0("outputs/TestResults/X_IATSE_n_", n, ".rds"))
  Y_IATSE <- readRDS(paste0("outputs/TestResults/Y_IATSE_n_", n, ".rds"))
  Z_IATSE <- readRDS(paste0("outputs/TestResults/Z_IATSE_n_", n, ".rds"))
  
  # --- Omnibus ---
  omnibus <- list(
    TWT = list(
      X = omnibus_power(X_TWT, alpha = 0.05),
      Y = omnibus_power(Y_TWT, alpha = 0.05),
      Z = omnibus_power(Z_TWT, alpha = 0.05)
    ),
    IATSE = list(
      X = omnibus_power(X_IATSE, alpha = 0.05),
      Y = omnibus_power(Y_IATSE, alpha = 0.05),
      Z = omnibus_power(Z_IATSE, alpha = 0.05)
    )
  )
  
  # --- Sensitivity ---
  sensitivity_res <- list(
    TWT = list(
      X = sensitivity(X_TWT, D1_X, alpha = 0.05),
      Y = sensitivity(Y_TWT, D1_Y, alpha = 0.05),
      Z = sensitivity(Z_TWT, D1_Z, alpha = 0.05)
    ),
    IATSE = list(
      X = sensitivity(X_IATSE, D1_X, alpha = 0.05),
      Y = sensitivity(Y_IATSE, D1_Y, alpha = 0.05),
      Z = sensitivity(Z_IATSE, D1_Z, alpha = 0.05)
    )
  )
  
  # --- ROI ---
  ROI <- list(
    TWT = list(
      X = lapply(D1p_X_ind, function(ind) ROI_power(X_TWT, ind, alpha = 0.05)),
      Y = lapply(D1p_Y_ind, function(ind) ROI_power(Y_TWT, ind, alpha = 0.05)),
      Z = lapply(D1p_Z_ind, function(ind) ROI_power(Z_TWT, ind, alpha = 0.05))
    ),
    IATSE = list(
      X = lapply(D1p_X_ind, function(ind) ROI_power(X_IATSE, ind, alpha = 0.05)),
      Y = lapply(D1p_Y_ind, function(ind) ROI_power(Y_IATSE, ind, alpha = 0.05)),
      Z = lapply(D1p_Z_ind, function(ind) ROI_power(Z_IATSE, ind, alpha = 0.05))
    )
  )
  
  
  # --- ROI ---
  FPower <- list(
    TWT = list(
      X = lapply(D1p_X_ind, function(ind) sensitivity(X_TWT, ind, alpha = 0.05)),
      Y = lapply(D1p_Y_ind, function(ind) sensitivity(Y_TWT, ind, alpha = 0.05)),
      Z = lapply(D1p_Z_ind, function(ind) sensitivity(Z_TWT, ind, alpha = 0.05))
    ),
    IATSE = list(
      X = lapply(D1p_X_ind, function(ind) sensitivity(X_IATSE, ind, alpha = 0.05)),
      Y = lapply(D1p_Y_ind, function(ind) sensitivity(Y_IATSE, ind, alpha = 0.05)),
      Z = lapply(D1p_Z_ind, function(ind) sensitivity(Z_IATSE, ind, alpha = 0.05))
    )
  )
  
  # --- Store everything ---
  results[[as.character(n)]] <- list(
    omnibus = omnibus,
    sensitivity = sensitivity_res,
    ROI = ROI,
    FPower = FPower
  )
}

#print(results$`20`)



##################Make data frame for plotting#################
library(dplyr)
library(tidyr)
library(purrr)

roi_df <- map_dfr(names(results), function(n) {
  
  res_n <- results[[n]]
  
  map_dfr(names(res_n$ROI), function(method) {        # TWT / IATSE
    
    map_dfr(names(res_n$ROI[[method]]), function(var) {  # X / Y / Z
      
      roi_list <- res_n$ROI[[method]][[var]]
      
      tibble(
        n = as.numeric(n),
        method = method,
        variable = var,
        tau = names(roi_list),
        value = unlist(roi_list)
      )
      
    })
  })
})


FPower_df <- map_dfr(names(results), function(n) {
  
  res_n <- results[[n]]
  
  map_dfr(names(res_n$FPower), function(method) {        # TWT / IATSE
    
    map_dfr(names(res_n$FPower[[method]]), function(var) {  # X / Y / Z
      
      FP_list <- res_n$FPower[[method]][[var]]
      
      tibble(
        n = as.numeric(n),
        method = method,
        variable = var,
        tau = names(FP_list),
        value = unlist(FP_list)
      )
      
    })
  })
})


library(ggplot2)

ggplot(roi_df, aes(x = n, y = value, color = tau, group = tau)) + geom_line() +
  facet_grid(method ~ variable) + theme_minimal() + ylim(0, 1)

ggplot(roi_df, aes(x = n, y = value, color = tau, group = tau)) +
  geom_line(position = position_jitter(width = 0.005, height = 0.005)) +
  facet_grid(method ~ variable) +
  theme_minimal()

ggplot(FPower_df, aes(x = n, y = value, color = tau, group = tau)) +
  geom_line() +
  facet_grid(method ~ variable) +
  theme_minimal() +
  ylim(0, 1)

omnibus_df <- map_dfr(names(results), function(n) {
  
  res_n <- results[[n]]
  
  map_dfr(names(res_n$omnibus), function(method) {
    
    map_dfr(names(res_n$omnibus[[method]]), function(var) {
      
      tibble(
        n = as.numeric(n),
        method = method,
        variable = var,
        value = res_n$omnibus[[method]][[var]]
      )
      
    })
  })
})

sensitivity_df <- map_dfr(names(results), function(n) {
  
  res_n <- results[[n]]
  
  map_dfr(names(res_n$sensitivity), function(method) {
    
    map_dfr(names(res_n$sensitivity[[method]]), function(var) {
      
      tibble(
        n = as.numeric(n),
        method = method,
        variable = var,
        value = res_n$sensitivity[[method]][[var]]
      )
      
    })
  })
})


omnibus_df <- omnibus_df %>%
  mutate(metric = "omnibus")

sensitivity_df <- sensitivity_df %>%
  mutate(metric = "sensitivity")
combined_df <- bind_rows(omnibus_df, sensitivity_df)


ggplot(combined_df, aes(x = n, y = value, color = metric)) +
  geom_line() +
  facet_grid(method ~ variable) +
  theme_minimal() +
  ylim(0, 1)

power_df <- data.frame(
  n = rep(n_samples, times = 6),  # 6 vectors: t-test X/Y/Z + Wilcox X/Y/Z
  value = c(power_ttest_X, power_ttest_Y, power_ttest_Z,
            power_wilcox_X, power_wilcox_Y, power_wilcox_Z),
  method = rep(c("t-test", "Wilcoxon"), each = length(n_samples) * 3),
  variable = rep(rep(c("X", "Y", "Z"), each = length(n_samples)), times = 2),
  metric = "power"
)

matplot(n_samples, power_ttest_X, type = "l", col = "blue", ylim = c(0, 1), xlab = "Sample Size", ylab = "Power", main = "Power of t-test")

ggplot(power_df, aes(x = n, y = value, color = method)) +
  geom_line() +
  facet_grid(metric ~ variable) +
  theme_minimal() +
  ylim(0, 1)
