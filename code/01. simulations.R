#----------------------------------------------------------------------------------
# Libraries 
#----------------------------------------------------------------------------------

# df manipulation:
library(readxl)
library(dplyr)
library(tidyr)
library(reshape2)
library(patchwork)  # For combining multiple ggplots


# Plotting:
library(ggplot2)
library(Amelia) # For missing value maps visualization

# Wrappers:
source('Impute_wrapper.R')

# Load cleaned df
V5_unique_species <- readRDS("V5_unique_species.rds")

#----------------------------------------------------------------------------------
# Simulation of missing values MCAR (PER COLUMN) to a complete df
#----------------------------------------------------------------------------------
# Function to introduce missing data (MCAR) into a dataset (returns MCAR df and indices) per column
MCAR_generate_per_column <- function(data, missing_prop = 0.05) {
  
  # Copy the dataset to avoid modifying the original
  data_res <- data
  
  # Initialize a list to store indices of missing values for each column
  missing_idx_per_column <- list()
  
  # Loop over each column in the dataset
  for (col in 1:ncol(data_res)) {
    
    # Get the indices of non-NA values in the current column
    non_na_idx <- which(!is.na(data_res[, col]))
    
    # Randomly select a proportion of these indices to set as missing
    num_missing <- round(length(non_na_idx) * missing_prop)
    
    if (num_missing > 0) {
      random_idx <- sample(non_na_idx, num_missing, replace = FALSE)
      
      # Set selected indices to NA in the current column
      data_res[random_idx, col] <- NA
      
      # Store the indices of the missing values for the current column
      missing_idx_per_column[[paste0("col_", col)]] <- random_idx
    }
  }
  
  return(list(data_MCAR = data_res, MCAR_idx = missing_idx_per_column))
}

# Run the function with missing proportions from 0.05 to 0.4 in increments of 0.05, repeated num_repeats times
missing_props <- seq(0.05, 0.45, by=0.1)
num_repeats <- 10 
mcar_datasets <- list()  # List to store datasets
mcar_indices <- list()   # List to store indices

# Loop over each missing proportion
for (prop in missing_props) {
  for (rep in 1:num_repeats) {
    # Run MCAR_generate_per_column for each repetition
    MCAR_result <- MCAR_generate_per_column(V5_unique_species, missing_prop = prop)
    
    # Extract the data with missing values
    data_MCAR <- MCAR_result$data_MCAR
    dataset_name <- paste0("MCAR_", prop, "_rep_", rep)
    mcar_datasets[[dataset_name]] <- data_MCAR  # Store the result in a named list
    
    # Extract and store the missing value indices
    MCAR_idx <- MCAR_result$MCAR_idx
    index_name <- paste0("MCAR_", prop, "_rep_", rep, "_idx")
    mcar_indices[[index_name]] <- MCAR_idx  # Store the indices in a named list
  }
}



rm(MCAR_idx, data_MCAR, MCAR_result)

# Set up the plot area to have 4 panels (2 rows x 2 columns)
par(mfrow = c(2, 2))  # 2 rows, 2 columns

# Plot 4 datasets in the same figure
for (i in seq(6,40,10)) {
  name <- names(mcar_datasets)[i]
  cat("Plotting missing map for: ", name, "\n")
  plot_missing_map(mcar_datasets[[name]], name)
}
saveRDS(mcar_datasets,file='mcar_datasets.rds')

# Reset the plot area back to 1 panel 
par(mfrow = c(1, 1))

#----------------------------------------------------------------------------------
# Simulation of missing values MNAR to a complete df for prop_miss PER COLUMN
#----------------------------------------------------------------------------------
### We can vary up to 40% of columns where there will be missing
# It is important to choose several different columns (run the experiment several times to not be biased by the data in 1 column) 
# MNAR generation function
MNAR_generate <- function(data, var_prop = 0.05) {
  # var_prop: A proportion used to determine the cutoff for introducing missing values.
  
  data_mis <- data
  
  # Loop through all columns
  for (i in 1:ncol(data)) {
    # Get column data
    current_column_data <- data[, i]
    
    # Calculate a cutoff value for the current column using the current var_prop.
    # Values will be ordered increasingly and we set the cutoff where the proportion value falls.
    cutoff <- quantile(current_column_data, var_prop, na.rm = TRUE)
    
    # Set the values under the cutoff as NA
    current_column_data[current_column_data < cutoff] <- NA
    data_mis[, i] <- current_column_data
  }
  
  # Return the dataframe with missing values
  return(data_mis)
}

# Define the different var_prop values
var_prop_values <- c(0.05, 0.1, 0.2, 0.3, 0.4)

# Create a list to store the dataframes with missing values
MNAR_datasets <- list()

# Loop through each var_prop value
for (prop in var_prop_values) {
  # Generate the dataset with missing values for the current var_prop
  MNAR_data <- MNAR_generate(V5_unique_species, var_prop = prop)
  
  # Store the result in the list with a descriptive name
  MNAR_datasets[[paste0("MNAR_var_prop_", prop)]] <- MNAR_data
}

saveRDS(MNAR_datasets,file='MNAR_datasets.rds')


# The MNAR_datasets list now contains 5 dataframes, each corresponding to a different var_prop value

# Plotting some of the generated datasets (optional)
par(mfrow = c(2, 2))  # Set up the plot area for 2 rows and 2 columns

for (i in seq(2,5)) {
  name <- names(MNAR_datasets)[i]
  cat("Plotting missing map for: ", name, "\n")
  plot_missing_map(MNAR_datasets[[name]], name)
}

# Reset the plot area to 1 panel
par(mfrow = c(1, 1))

#----------------------------------------------------------------------------------
# Checking distributions of datasets
#----------------------------------------------------------------------------------
# MCAR Define the missing proportions to plot
proportions <- c(0.2, 0.4)

# Loop through the missing proportions
for (prop in proportions) {
  
  # Work only with the first repetition for each proportion (rep_1)
  dataset_name <- paste0("MCAR_", prop, "_rep_1")
  
  # Extract the corresponding dataset with missing values
  missing_data <- mcar_datasets[[dataset_name]]
  
  # Find the columns in the dataset where there are missing values
  cols_with_missing <- which(colSums(is.na(missing_data)) > 0)
  
  # Ensure we have at least 3 columns with missing data
  if (length(cols_with_missing) < 3) {
    stop(paste("Not enough columns with missing data for proportion", prop))
  }
  
  # Randomly select 3 columns from those that contain missing values
  selected_columns <- sample(cols_with_missing, 3)
  
  # Extract the full dataset with no missing values (V5_unique_species)
  full_data_subset <- V5_unique_species[, selected_columns, drop = FALSE]
  
  # Subset the missing data to these selected columns
  missing_data_subset <- missing_data[, selected_columns, drop = FALSE]
  
  # Melt and label the full dataset (no missing values)
  melted_full <- melt(as.data.frame(full_data_subset), variable.name = "Variable", value.name = "Value")
  melted_full$Type <- "Complete"
  
  # Melt and label the missing dataset (with MCAR missingness)
  melted_missing <- melt(as.data.frame(missing_data_subset), variable.name = "Variable", value.name = "Value")
  melted_missing$Type <- "MCAR"
  
  # Combine both datasets
  combined_data <- rbind(melted_full, melted_missing)
  
  # Add the proportion and repetition for easy labeling
  combined_data$Prop <- prop
  combined_data$Rep <- 1
  
  # Print a density plot for this proportion
  print(
    ggplot(combined_data, aes(x = Value, fill = Type, color = Type)) +
      geom_density(alpha = 0.5) +
      facet_grid(Variable ~ ., scales = "free") +  # Facet by Variable, 3 panels
      labs(title = paste("Density Comparison: Complete vs MCAR Simulated Data (Proportion:", prop, ")"),
           x = "Value", y = "Density") +
      theme_minimal() +
      theme(legend.position = "right")
  )
}

# Plot histograms to confirm
par(mfrow = c(2, 2)) 
mcar_0.2 <- mcar_datasets[["MCAR_0.2_rep_1"]]
mcar_0.4 <- mcar_datasets[["MCAR_0.4_rep_1"]]
plot_missing_proportions(mcar_0.2)
plot_missing_proportions(mcar_0.4)

# MNAR
# Loop through the missing proportions
for (i in 1:5) {
  
  # Extract the corresponding dataset with missing values
  missing_data <- MNAR_datasets[[i]]
  
  # Find the columns in the dataset where there are missing values
  cols_with_missing <- which(colSums(is.na(missing_data)) > 0)
  
  # Randomly select 3 columns from those that contain missing values
  selected_columns <- sample(cols_with_missing, 3)
  
  # Extract the full dataset with no missing values (V5_unique_species)
  full_data_subset <- V5_unique_species[, selected_columns, drop = FALSE]
  
  # Subset the missing data to these selected columns
  missing_data_subset <- missing_data[, selected_columns, drop = FALSE]
  
  # Melt and label the missing dataset (with MNAR missingness)
  melted_missing <- melt(as.data.frame(missing_data_subset), variable.name = "Variable", value.name = "Value")
  melted_missing$Type <- "MNAR"
  
  # Melt and label the full dataset (no missing values)
  melted_full <- melt(as.data.frame(full_data_subset), variable.name = "Variable", value.name = "Value")
  melted_full$Type <- "Complete"
  
  # Combine both datasets
  combined_data <- rbind(melted_missing, melted_full)
  
  # Add the proportion and repetition for easy labeling
  combined_data$Prop <- prop
  combined_data$Rep <- 1
  
  # Print a density plot for this proportion
  print(
    ggplot(combined_data, aes(x = Value, fill = Type, color = Type)) +
      geom_density(alpha = 0.5) +
      facet_grid(Variable ~ ., scales = "free") +  # Facet by Variable, 3 panels
      labs(title = paste("Density Comparison: Complete vs Simulated Missing Data", names(MNAR_datasets)[i]),
           x = "Value", y = "Density") +
      theme_minimal() +
      theme(legend.position = "right")
  )
}

par(mfrow = c(2, 2)) 
mnar_0.2 <- MNAR_datasets[[3]]
mnar_0.4 <- MNAR_datasets[[5]]
plot_missing_proportions(mnar_0.2)
plot_missing_proportions(mnar_0.4)


