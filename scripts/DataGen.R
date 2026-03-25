rm(list = ls())
print(getwd())
print(.libPaths())

if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}

print(.libPaths())

source("R/utilities.R")
source("R/ThoraxDataParameters.R")


results <- compute_metrics(
  "data/Thorax_X_101.xlsx",
  "data/Thorax_Y_101.xlsx",
  "data/Thorax_Z_101.xlsx"
)

mean_X_male <- apply(results$row_means$X_male, 1, mean)
mean_X_female <- apply(results$row_means$X_female, 1, mean)
mean_Y_male <- apply(results$row_means$Y_male, 1, mean)
mean_Y_female <- apply(results$row_means$Y_female, 1, mean)
mean_Z_male <- apply(results$row_means$Z_male, 1, mean)
mean_Z_female <- apply(results$row_means$Z_female, 1, mean)
fwhm_X <- results$fwhm$X
fwhm_Y <- results$fwhm$Y
fwhm_Z <- results$fwhm$Z
sd_X <- results$sd$X
sd_Y <- results$sd$Y
sd_Z <- results$sd$Z

set.seed(24032026)

n_samples <- 5:50
MC_size <- 10000

con_size <- length(sd_X)

for (n in n_samples) {
  
  X_array_male <- matrix(NA, nrow = con_size, ncol = n * MC_size)
  X_array_female <- matrix(NA, nrow = con_size, ncol = n * MC_size)
  Y_array_male <- matrix(NA, nrow = con_size, ncol = n * MC_size)
  Y_array_female <- matrix(NA, nrow = con_size, ncol = n * MC_size)
  Z_array_male <- matrix(NA, nrow = con_size, ncol = n * MC_size)
  Z_array_female <- matrix(NA, nrow = con_size, ncol = n * MC_size)
  
  col_idx <- 1
  
  for (mc in 1:MC_size) {
    
    # if (mc %% 500 == 0){
    #   print(paste("Generating data for n =", n, "MC iteration:", mc))
    # }
    
    # Generate noise per MC iteration
    noise_X_male <- Noise_generator(n, con_size, 0, sd_X, fwhm_X)$noise
    noise_X_female <- Noise_generator(n, con_size, 0, sd_X, fwhm_X)$noise
    
    noise_Y_male <- Noise_generator(n, con_size, 0, sd_Y, fwhm_Y)$noise
    noise_Y_female <- Noise_generator(n, con_size, 0, sd_Y, fwhm_Y)$noise
    
    noise_Z_male <- Noise_generator(n, con_size, 0, sd_Z, fwhm_Z)$noise
    noise_Z_female <- Noise_generator(n, con_size, 0, sd_Z, fwhm_Z)$noise
    
    # Generate data
    X_male <- data_generator(signal = mean_X_male, noise = noise_X_male)
    X_female <- data_generator(signal = mean_X_female, noise = noise_X_female)
    
    Y_male <- data_generator(signal = mean_Y_male, noise = noise_Y_male)
    Y_female <- data_generator(signal = mean_Y_female, noise = noise_Y_female)
    
    Z_male <- data_generator(signal = mean_Z_male, noise = noise_Z_male)
    Z_female <- data_generator(signal = mean_Z_female, noise = noise_Z_female)
    
    
    cols <- col_idx:(col_idx + n - 1)
    
    X_array_male[, cols] <- X_male
    X_array_female[, cols] <- X_female
    Y_array_male[, cols] <- Y_male
    Y_array_female[, cols] <- Y_female
    Z_array_male[, cols] <- Z_male
    Z_array_female[, cols] <- Z_female
    
    
    col_idx <- col_idx + n
    
    #MEMORY CLEANUP every 50 iterations
    if (mc %% 50 == 0) {
      rm(
        noise_X_male, noise_X_female,
        noise_Y_male, noise_Y_female,
        noise_Z_male, noise_Z_female,
        X_male, X_female,
        Y_male, Y_female,
        Z_male, Z_female
      )
      gc()
    }
  }
  
  saveRDS(X_array_male, paste0("data/GeneratedData/X_male_n_", n, ".rds"), compress = "xz")
  saveRDS(X_array_female, paste0("data/GeneratedData/X_female_n_", n, ".rds"), compress = "xz")
  saveRDS(Y_array_male, paste0("data/GeneratedData/Y_male_n_", n, ".rds"), compress = "xz")
  saveRDS(Y_array_female, paste0("data/GeneratedData/Y_female_n_", n, ".rds"), compress = "xz")
  saveRDS(Z_array_male, paste0("data/GeneratedData/Z_male_n_", n, ".rds"), compress = "xz")
  saveRDS(Z_array_female, paste0("data/GeneratedData/Z_female_n_", n, ".rds"), compress = "xz")
}




