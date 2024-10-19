# Imputation Wrapper ------------------------------------------------------
require(missForest)
require(impute)
require(magrittr)
require(imputeLCMD)
# source('MVI_global.R')

# Packages installed with BiocdManager: impute, pcaMethods

RF_wrapper <- function(data, ...) {
  result <- missForest(data, ...)[[1]]
  return (result)
}

kNN_wrapper <- function(data, ...) {
  result <- data %>% data.matrix %>% impute.knn(., ...) %>% extract2(1)
  return(result)
}

SVD_wrapper <- function(data, K = 5) {
  data_sc_res <- scale_recover(data, method = 'scale')
  data_sc <- data_sc_res[[1]]
  data_sc_param <- data_sc_res[[2]]
  result <- data_sc %>% impute.wrapper.SVD(., K = K) %>% 
    scale_recover(., method = 'recover', param_df = data_sc_param) %>% extract2(1)
  return(result)
}

Mean_wrapper <- function(data) {
  result <- data
  result[] <- lapply(result, function(x) {
    x[is.na(x)] <- mean(x, na.rm = T)
    x
  })
  return(result)
}

Median_wrapper <- function(data) {
  result <- data
  result[] <- lapply(result, function(x) {
    x[is.na(x)] <- median(x, na.rm = T)
    x
  })
  return(result)
}


Min_0.1_wrapper <- function(data) {
  result <- data
  result[] <- lapply(result, function(x) {
    x[is.na(x)] <- min(x, na.rm = T)/10
    x
  })
  return(result)
}

Min_0.2_wrapper <- function(data) {
  result <- data
  result[] <- lapply(result, function(x) {
    x[is.na(x)] <- min(x, na.rm = T)/5
    x
  })
  return(result)
}

Min_0.4_wrapper <- function(data) {
  result <- data
  result[] <- lapply(result, function(x) {
    x[is.na(x)] <- min(x, na.rm = T)*0.4
    x
  })
  return(result)
}

Min_0.5_wrapper <- function(data) {
  result <- data
  result[] <- lapply(result, function(x) {
    x[is.na(x)] <- min(x, na.rm = T)/2
    x
  })
  return(result)
}


QRILC_wrapper <- function(data, ...) {
  result <- data %>% log %>% impute.QRILC(., ...) %>% extract2(1) %>% exp
  return(result)
}


# NRMSE Calculation Function by Column (MSE for the whole column)
nrmse_by_column <- function(imputed_df, full_df) {
  # Ensure that the datasets have the same dimensions
  if (!all(dim(imputed_df) == dim(full_df))) {
    stop("The dimensions of the imputed and full datasets do not match.")
  }
  
  # Initialize a vector to store NRMSE for each column
  nrmse_values <- numeric(ncol(imputed_df))
  names(nrmse_values) <- colnames(imputed_df)
  
  # Loop through each column
  for (col in 1:ncol(imputed_df)) {
    # Calculate the mean squared error for the entire column (all values)
    mse <- mean((imputed_df[, col] - full_df[, col])^2, na.rm = TRUE)
    
    # Calculate the variance of the entire column
    variance_whole <- var(full_df[, col], na.rm = TRUE)
    
    # If the variance is zero (all values in the column are the same), set variance to 1 to avoid division by zero
    if (variance_whole == 0) {
      variance_whole <- 1
    }
    
    # Calculate NRMSE (normalized by variance of the whole column)
    nrmse_values[col] <- sqrt(mse / variance_whole)
  }
  
  return(nrmse_values)
}

# Define the Weighted NRMSE Function
weighted_nrmse_by_column <- function(imputed_df, full_df, missing_df) {
  # Initialize a vector to store weighted NRMSE for each column
  weighted_nrmse_values <- numeric(ncol(imputed_df))
  names(weighted_nrmse_values) <- colnames(imputed_df)
  
  for (col in 1:ncol(imputed_df)) {
    # Get the proportion of missing values for this column
    missing_prop <- sum(is.na(missing_df[, col])) / nrow(missing_df)
    
    # Calculate the mean squared error for the entire column (including observed + imputed)
    mse <- mean((imputed_df[, col] - full_df[, col])^2, na.rm = TRUE)
    
    # Calculate the variance of the full column (non-imputed)
    variance_whole <- var(full_df[, col], na.rm = TRUE)
    
    # If variance is zero, set it to 1 to avoid division by zero
    if (variance_whole == 0) variance_whole <- 1
    
    # Calculate Weighted NRMSE
    nrmse <- sqrt(mse / variance_whole)
    weighted_nrmse_values[col] <- nrmse * missing_prop  # Weight NRMSE by missing proportion
  }
  
  return(weighted_nrmse_values)
}
