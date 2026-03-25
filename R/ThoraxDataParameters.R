library(readxl)
source("R/utilities.R")

compute_metrics <- function(path_X, path_Y, path_Z) {
  
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
  
  resid_X <- residuals_data(t(row_means_X_male), t(row_means_X_female))
  resid_Y <- residuals_data(t(row_means_Y_male), t(row_means_Y_female))
  resid_Z <- residuals_data(t(row_means_Z_male), t(row_means_Z_female))
  
  fwhm_estimate_X <- estimate_fwhm(resid_X)
  fwhm_estimate_Y <- estimate_fwhm(resid_Y)
  fwhm_estimate_Z <- estimate_fwhm(resid_Z)
  
  sd_X <- pooled_sd(row_means_X_male, row_means_X_female)
  sd_Y <- pooled_sd(row_means_Y_male, row_means_Y_female)
  sd_Z <- pooled_sd(row_means_Z_male, row_means_Z_female)
  
  return(list(
    row_means = list(
      X_male = row_means_X_male,
      X_female = row_means_X_female,
      Y_male = row_means_Y_male,
      Y_female = row_means_Y_female,
      Z_male = row_means_Z_male,
      Z_female = row_means_Z_female
    ),
    fwhm = list(
      X = fwhm_estimate_X,
      Y = fwhm_estimate_Y,
      Z = fwhm_estimate_Z
    ),
    sd = list(
      X = sd_X,
      Y = sd_Y,
      Z = sd_Z
    )
  ))
}


