rm(list=ls())
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

# tau_values = seq(0.1, 1, by = 0.1)

# compute_D1p(mean_male_X, mean_female_X, sd_male_X, sd_female_X,
#               n_X_male, n_X_female, tau_values)
# compute_D1p(mean_male_Y, mean_female_Y, sd_male_Y, sd_female_Y,
#               n_Y_male, n_Y_female, tau_values)
# compute_D1p(mean_male_Z, mean_female_Z, sd_male_Z, sd_female_Z,
#               n_Z_male, n_Z_female, tau_values)

POI_X <- which.max(SAE_function(mean_male_X, mean_female_X, sd_male_X, sd_female_X,
                                n_X_male, n_X_female))
POI_Y <- which.max(SAE_function(mean_male_Y, mean_female_Y, sd_male_Y, sd_female_Y,
                                n_Y_male, n_Y_female))
POI_Z <- which.max(SAE_function(mean_male_Z, mean_female_Z, sd_male_Z, sd_female_Z,
                                n_Z_male, n_Z_female))

n_samples <- 5:50
MC_size <- 10000
con_size <- 101
X_array_ttest <- matrix(NA, nrow = length(n_samples), ncol = MC_size)
Y_array_ttest <- matrix(NA, nrow = length(n_samples), ncol = MC_size)
Z_array_ttest <- matrix(NA, nrow = length(n_samples), ncol = MC_size)

X_array_wilcox <- matrix(NA, nrow = length(n_samples), ncol = MC_size)
Y_array_wilcox <- matrix(NA, nrow = length(n_samples), ncol = MC_size)
Z_array_wilcox <- matrix(NA, nrow = length(n_samples), ncol = MC_size)
for (n in n_samples) {
  
  X_male <- readRDS(paste0("data/GeneratedData/X_male_n_", n, ".rds"))[POI_X,]
  X_female <- readRDS(paste0("data/GeneratedData/X_female_n_", n, ".rds"))[POI_X,]
  Y_male <- readRDS(paste0("data/GeneratedData/Y_male_n_", n, ".rds"))[POI_Y,]
  Y_female <- readRDS(paste0("data/GeneratedData/Y_female_n_", n, ".rds"))[POI_Y,]
  Z_male <- readRDS(paste0("data/GeneratedData/Z_male_n_", n, ".rds"))[POI_Z,]
  Z_female <- readRDS(paste0("data/GeneratedData/Z_female_n_", n, ".rds"))[POI_Z,]
  
  col_idx <- 1
  sample_counter <- which(n_samples == n)
  
  for (mc in 1:MC_size) {
    
    if (mc %% 5000 == 0){
      print(paste("P-value for n =", n, "MC iteration:", mc))
    }
    
    cols <- col_idx:(col_idx + n - 1)
    
    X_array_ttest[sample_counter, mc] <- t.test(X_male[cols], X_female[cols], var.equal = TRUE)$p.value
    Y_array_ttest[sample_counter, mc] <- t.test(Y_male[cols], Y_female[cols], var.equal = TRUE)$p.value
    Z_array_ttest[sample_counter, mc] <- t.test(Z_male[cols], Z_female[cols], var.equal = TRUE)$p.value
    
    X_array_wilcox[sample_counter, mc] <- wilcox.test(X_male[cols], X_female[cols])$p.value
    Y_array_wilcox[sample_counter, mc] <- wilcox.test(Y_male[cols], Y_female[cols])$p.value
    Z_array_wilcox[sample_counter, mc] <- wilcox.test(Z_male[cols], Z_female[cols])$p.value
    
    col_idx <- col_idx + n
  }
}
#check
length(n_samples)
dim(X_array_wilcox)


#write
#rows are sample sizes, columns are MC iterations
saveRDS(X_array_ttest, paste0("outputs/TestResults/X_ttest", ".rds"), compress = "xz")
saveRDS(Y_array_ttest, paste0("outputs/TestResults/Y_ttest", ".rds"), compress = "xz")
saveRDS(Z_array_ttest, paste0("outputs/TestResults/Z_ttest", ".rds"), compress = "xz")

saveRDS(X_array_wilcox, paste0("outputs/TestResults/X_wilcox", ".rds"), compress = "xz")
saveRDS(Y_array_wilcox, paste0("outputs/TestResults/Y_wilcox", ".rds"), compress = "xz")
saveRDS(Z_array_wilcox, paste0("outputs/TestResults/Z_wilcox", ".rds"), compress = "xz")
