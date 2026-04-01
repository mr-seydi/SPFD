rm(list = ls())
library(fda)
source("R/TWT2.R")
source("R/utilities.R")

set.seed(26032026)

n_samples <- 5:50
MC_size <- 10000
con_size <- 101

for (n in n_samples) {
  
  X_array <- matrix(NA, nrow = con_size, ncol = MC_size)
  Y_array <- matrix(NA, nrow = con_size, ncol = MC_size)
  Z_array <- matrix(NA, nrow = con_size, ncol = MC_size)
  
  X_male <- readRDS(paste0("data/GeneratedData/X_male_n_", n, ".rds"))
  X_female <- readRDS(paste0("data/GeneratedData/X_female_n_", n, ".rds"))
  Y_male <- readRDS(paste0("data/GeneratedData/Y_male_n_", n, ".rds"))
  Y_female <- readRDS(paste0("data/GeneratedData/Y_female_n_", n, ".rds"))
  Z_male <- readRDS(paste0("data/GeneratedData/Z_male_n_", n, ".rds"))
  Z_female <- readRDS(paste0("data/GeneratedData/Z_female_n_", n, ".rds"))
  
  col_idx <- 1
  
  for (mc in 1:MC_size) {
    
    if (mc %% 500 == 0){
      print(paste("P-value for n =", n, "MC iteration:", mc))
    }
    
    cols <- col_idx:(col_idx + n - 1)
    
    X_array[,mc] <- TWT(t(X_male[, cols]), t(X_female[, cols]))
    Y_array[,mc] <- TWT(t(Y_male[, cols]), t(Y_female[, cols]))
    Z_array[,mc] <- TWT(t(Z_male[, cols]), t(Z_female[, cols]))
    
    
    col_idx <- col_idx + n
  }
  
  saveRDS(X_array, paste0("outputs/TestResults/X_TWT_n_", n, ".rds"), compress = "xz")
  saveRDS(Y_array, paste0("outputs/TestResults/Y_TWT_n_", n, ".rds"), compress = "xz")
  saveRDS(Z_array, paste0("outputs/TestResults/Z_TWT_n_", n, ".rds"), compress = "xz")
  
  # Cleanup AFTER finishing n
  rm(X_male, X_female, Y_male, Y_female, Z_male, Z_female)
  gc()
}

