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

#----------------------------------------------------------------------------------
# Functions for Plotting
#----------------------------------------------------------------------------------
# Plotting missing values map
plot_missing_map <- function(df, df_name, ...) {
  # Create a dynamic title using the passed dataset name
  main_title <- sprintf("Missing values %s", df_name)
  
  # Plot the missing values matrix for raw data
  Amelia::missmap(df, main = main_title,
                  col = c("red", "black"), legend = FALSE,
                  rank.order = FALSE, axes = TRUE, ... )
}

# Histogram of proportions of Missing values
plot_missing_proportions <- function(df, bin_width = 0.1) {
  
  df_name <- deparse(substitute(df))
  
  # Calculate the proportion of missing values for each column
  total_rows <- nrow(df)
  missing_proportions <- sapply(df, function(x) sum(is.na(x)) / total_rows)
  
  # Create a data frame for ggplot
  df_missing <- data.frame(proportion_missing = missing_proportions)
  
  # Plot using ggplot2 with binning
  ggplot(df_missing, aes(x = proportion_missing)) +
    geom_histogram(binwidth = bin_width, fill = '#00BFC4', color = "white",alpha=0.5) +
    labs(title = paste0("Columns count by proportion of missing values in ",df_name),
         x = "Proportion of Missing Values", 
         y = "Number of Columns") +
    theme_minimal() +
    scale_x_continuous(breaks = seq(0, 1, by = bin_width))  # Adjust x-axis for better readability
}

# Example usage:
plot_missing_proportions(v0_raw_data, bin_width = 0.1)

#----------------------------------------------------------------------------------
# v0. Read the data and create data frame
#----------------------------------------------------------------------------------
# Read the data
v0_raw_data <- read_excel("Data/Master_Project_NIST_MissingData.xlsx", sheet = "Data")
# Convert to a regular data frame
v0_raw_data <- as.data.frame(v0_raw_data)
# Set the first column as row names
rownames(v0_raw_data) <- v0_raw_data[[1]]
# Remove the first column from the data frame (name)
v0_raw_data <- v0_raw_data[-1]
# First row is the molecule species
# Second row is the subspecies
# First column is the name of the sample / repetition 

# Extract the first row as species names
subspecies_names <- v0_raw_data[1, ]
# Remove the first row from the data frame
v0_raw_data <- v0_raw_data[-1, ]

# Convert all columns to numeric. 
v0_raw_data <- v0_raw_data %>% mutate(across(everything(), as.numeric))

plot_missing_map(v0_raw_data,'Raw data')

#----------------------------------------------------------------------------------
# v1. Delete columns with repeated integer values (standards)
#----------------------------------------------------------------------------------
# Function to check if a column has repeated integers
is_repeated_integer <- function(column) {
  
  # Convert the column to numeric values, omitting any NA (missing) values.
  column_values <- na.omit(as.numeric(column))
  
  # If all values in column were NA or non-numeric, length will be 0 
  # so column should not be eliminated in this step. 
  if (length(column_values) == 0) return(FALSE)
  
  # Check if values are exact (they don't change if rounding)
  return(all(column_values == round(column_values)))
}

# Identify the columns to remove
temp_columns_to_remove <- sapply(v0_raw_data, is_repeated_integer)

# Remove the identified columns
v1_data_wo_standards <- v0_raw_data[, !temp_columns_to_remove]
plot_missing_map(v1_data_wo_standards,'v1_data_wo_standards')

# Remove temporary variables
rm(temp_columns_to_remove, is_repeated_integer)

#----------------------------------------------------------------------------------
# v2. Eliminate species(molecules) where all rows are NA for sub-species. 
#----------------------------------------------------------------------------------
# Function to eliminate columns where all rows are NA in the whole dataframe
eliminate_na_columns <- function(df) {
  # Identify columns where not all values are NA
  keep_columns <- which(colSums(!is.na(df)) > 0)
  
  # Subset the dataframe to keep only the columns where there is at least one data point available
  df_filtered <- df[, keep_columns]
  
  return(df_filtered)
}

# Apply the function to the dataset
v2_data_deleted_NA_species <- eliminate_na_columns(v1_data_wo_standards)

plot_missing_map(v2_data_deleted_NA_species,'v2_data_deleted_NA_species')

#----------------------------------------------------------------------------------
# Extract names of molecules (species)
#----------------------------------------------------------------------------------
# Function to extract molecule names from column names (e.g., extract "DAG 32:1" from "DAG 32:1...5")
extract_molecule_name <- function(col_name) {
  return(sub("\\.\\.\\..*$", "", col_name))
  # The pattern is a regular expression that matches this: {'...'}{anything except a new line}{end}
  # In R, you need to double the backslash (\\) because backslashes also need to be escaped in R string literals. 
  # The asterisk (*) means zero or more occurrences of the preceding element. 
  # Together, .* matches any sequence of characters (including an empty sequence) following the three dots.
  # $: This matches the end of the string.
  
}

# Apply the extraction function to all column names
molecule_names <- sapply(names(v2_data_deleted_NA_species), extract_molecule_name)

#---------------------------------------------------------------------------------
# v3. If applies, fill NAs with existing values of other sub-species of the same molecule 
#----------------------------------------------------------------------------------
# ---->>>Confirm if values can be just filled up from other species (why are they missing??) 
# ---->>>The algorithm is made so that values are not filled in if there is nothing similar within two species (even if they belong to same molecule)

# Function to fill NAs with values from the same molecule's columns:
process_molecule_columns <- function(df, molecule) {
  tryCatch({
    # Subset the columns that belong to the same molecule
    molecule_columns <- which(molecule_names == molecule)
    
    # Skip processing if there is only one column
    if (length(molecule_columns) <= 1) {
      message(sprintf("Skipping molecule: %s (only one column)", molecule))
      return(df)  # Return df without changes
    }
    
    molecule_data <- df[, molecule_columns]  
    
    # Get the number of non-NA values for each column
    non_na_counts <- sapply(molecule_data, function(col) sum(!is.na(col)))
    
    # Sort columns by the number of non-NA values (descending)
    sorted_columns <- order(non_na_counts, decreasing = TRUE) # This gives a vector of indexes 
    
    shortest_col_index <- sorted_columns[length(sorted_columns)]
    short_column <- molecule_data[, shortest_col_index]
    skip_columns <- 0
    
    # Only try imputation if there are still NAs in the column and if we have not already tried all the columns
    while (any(is.na(short_column)) && skip_columns < length(sorted_columns))  {
      
      # We will iterate through the rest of the longer columns, starting by the longest
      for (i in 1:(length(sorted_columns) - 1)) {
        long_column <- molecule_data[,sorted_columns[i]]
        
        # Check if non-NA values are the same in both columns
        non_na_short_col <- short_column[!is.na(short_column)]
        non_na_long_col <- long_column[!is.na(short_column)]
        
        # Sometimes, the selected indexes are NA in the longer column
        if (any(is.na(non_na_long_col))) {
          non_na_indices <- !is.na(non_na_long_col)
          non_na_long_col <- non_na_long_col[non_na_indices]
          non_na_short_col <- non_na_short_col[non_na_indices]
        }
        
        # Only enter conditional if there are still non NA values to compare
        if (length(non_na_long_col) > 0 && all(non_na_short_col == non_na_long_col)) {
          # Fill NAs in the shorter column with values from the long column
          na_indices <- which(is.na(short_column))
          short_column[na_indices] <- long_column[na_indices]
          molecule_data[, shortest_col_index] <- short_column
        }
        
        # Break the loop if no more NAs in the short_column
        if (!any(is.na(short_column))) {
          break
        }
      }
      
      ## After comparing to all the 'longer' columns, we check if we were able to fill any missing value:
      # Get the number of non-NA values for each column
      non_na_counts <- sapply(molecule_data, function(col) sum(!is.na(col)))
      
      # Sort columns by the number of non-NA values (descending)
      sorted_columns <- order(non_na_counts, decreasing = TRUE) # This gives a vector of indexes 
      
      # Check if we were able to make any changes in the short column:
      shortest_col_index <- sorted_columns[length(sorted_columns) - skip_columns]
      new_short_column <- molecule_data[, shortest_col_index]
      # If nothing changed, we can skip this column for the next iteration
      if (identical(short_column, new_short_column)) {
        skip_columns <- skip_columns + 1
        shortest_col_index <- sorted_columns[length(sorted_columns) - skip_columns]
        short_column <- molecule_data[, shortest_col_index]
      } else {
        short_column <- new_short_column
      }
    }
    
    df[, molecule_columns] <- molecule_data
    return(df)
    
  }, error = function(e) {
    message(sprintf("Error encountered with molecule: %s", molecule))
    message(e)
    return(df)  # Return the original df in case of error
  })
}

# Iterate over all test molecules
v3_data_filled_by_mol <- v2_data_deleted_NA_species  # Initialize data_filled
for (molecule in unique(molecule_names)) {  
  v3_data_filled_by_mol <- process_molecule_columns(v3_data_filled_by_mol, molecule)
}
plot_missing_map(v3_data_filled_by_mol,'v3_data_filled_by_mol')

#---------> FX CASE TEST - process_molecule_columns()

# test_molecule <- unique(molecule_names)[35] #35 will not fill up columns, test 25 for full data
# 
# before_test_df <- v2_data_deleted_NA_species[,(which(molecule_names == test_molecule))]
# after_test_df <- process_molecule_columns(v2_data_deleted_NA_species,test_molecule)
# after_test_df <- after_test_df[,(which(molecule_names == test_molecule))]

#----------------------------------------------------------------------------------
# v4. FINAL CLEANUP - Eliminate columns(subspecies) where any rows are NA (optn for thres!=0)
#----------------------------------------------------------------------------------
# Function to eliminate columns where the number of non-NA values is less than or equal to a threshold
eliminate_na_species <- function(df, threshold = 0) {
  # Identify rows where the number of non-NA values is greater than the threshold
  keep_species <- which(colSums(is.na(df)) <= threshold)
  
  # Subset the dataframe to keep only the rows that meet the threshold condition
  df_filtered <- df[,keep_species ]
  
  return(df_filtered)
}

# Apply the function to the dataset
V4_data_only_nonNA <- eliminate_na_species(v3_data_filled_by_mol)
plot_missing_map(V4_data_only_nonNA,'V4_data_only_nonNA')


# keep subspecies with the most information: 
keep_unique_columns <- function(df) {
  #This function works by treating the columns as lists, checking for duplicated columns, 
  #and removing the duplicates while keeping the first occurrence of each unique column. 
  #The result will be a data frame with only unique columns.
  df_unique <- df[, !duplicated(as.list(df))]
  return(df_unique)
}

V5_unique_species <- keep_unique_columns(V4_data_only_nonNA)

# Get the current column names of both dataframes
# Check which columns in the first row are not equal to 'MS1'
columns_to_keep <- subspecies_names[1, ] != 'MS1'
subspecies_colnames <- colnames(subspecies_names)[columns_to_keep]

unique_columns_colnames <- colnames(V5_unique_species)

# Find the intersection of the column names in both dataframes
matching_columns <- intersect(subspecies_colnames, unique_columns_colnames)

# For the columns that are in both dataframes, rename them in V5_unique_species
# using the first row of subspecies_names
for (col in matching_columns) {
  # Get the new name from the first row of subspecies_names
  new_name <- subspecies_names[1, col]
  
  # Rename the corresponding column in V5_unique_species
  colnames(V5_unique_species)[colnames(V5_unique_species) == col] <- new_name
}

rm(subspecies_colnames,unique_columns_colnames,matching_columns,columns_to_keep)

# Save dataframe as CSV
write.csv(V5_unique_species, file = "full_cleaned_df.csv", row.names = FALSE)
# Save dataframe as RDS
saveRDS(V5_unique_species, file = "V5_unique_species.rds")

#----------------------------------------------------------------------------------
# Simulation of missing values MCAR to a complete df(based on Wei 2018)
#----------------------------------------------------------------------------------
# Function to introduce missing data (MCAR) into a dataset (returns MCAR df and indices)
MCAR_generate <- function(data, missing_prop = 0.05) {
  
  # Get a list of all the indices in the data frame (exclude NA values)
  all_idx <- which(!is.na(data), arr.ind = TRUE)
  
  # Random sampling of indices based on desired proportion
  random_idx <- sample(1:nrow(all_idx), round(nrow(all_idx) * missing_prop), replace = FALSE)
  selected_idx <- all_idx[random_idx, ]
  
  # Introduce missing values to the dataframe
  data_res <- data
  for (i in seq_len(nrow(selected_idx))) {
    row <- selected_idx[i, 1]
    col <- selected_idx[i, 2]
    data_res[row, col] <- NA
  }
  
  return(list(data_MCAR = data_res, MCAR_idx = selected_idx))
}
# Generate datasets with missing proportions from 0.05 to 0.4 in increments of 0.05, repeated num_repeats times
missing_props <- seq(0.05, 0.4, by = 0.05)
num_repeats <- 5  # Set the number of repetitions
mcar_datasets <- list()  # List to store datasets
mcar_indices <- list()   # List to store indices

# Loop over each missing proportion
for (prop in missing_props) {
  for (rep in 1:num_repeats) {
    # Run MCAR_generate for each repetition
    MCAR_result <- MCAR_generate(V5_unique_species, missing_prop = prop)
    
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

rm(MCAR_idx,data_MCAR,MCAR_result)

# Set up the plot area to have 4 panels (2 rows x 2 columns)
par(mfrow = c(2, 2))  # 2 rows, 2 columns


# Plot 4 datasets in the same figure
for (i in 1:8) {
  name <- names(mcar_datasets)[i]
  cat("Plotting missing map for: ", name, "\n")
  plot_missing_map(mcar_datasets[[name]], name)
}

# Reset the plot area back to 1 panel 
par(mfrow = c(1, 1))
#----------------------------------------------------------------------------------
# Simulation of missing values MNAR to a complete df (based on Wei 2018)
#----------------------------------------------------------------------------------
### We can vary up to 40% of columns where there will be missing
# It is important to choose several different columns (run the experiment several times to not be biased by the data in 1 column) 

# MNAR generation function
MNAR_generate <- function(data, missing_var_prop = 0.1, var_prop = seq(.05, .4, .05)){
  # missing_var_prop: Proportion of variables (columns) to introduce missing values in. 
  # var_prop: A sequence of proportions used to determine the cutoff for introducing missing values.
  
  data_mis <- data
  
  # Randomly select the column indices for the proportion of columns specified
  var_mis_list <- sample(1:ncol(data), round(ncol(data)*missing_var_prop)) 
  
  # Loop through selected columns and introduce missing variables
  for (i in 1:length(var_mis_list)) {
    # Get column data
    variable_index <- var_mis_list[i]
    current_column_data <- data[, variable_index]
    
    # Calculate a cutoff value for the current column using a randomly selected proportion from var_prop.
    # Values will be ordered increasingly and we set the cutoff where the proportion value falls.
    cutoff <- quantile(current_column_data, sample(var_prop, 1))
    
    # Set the values under the cutoff as NA
    current_column_data[current_column_data < cutoff] <- NA
    data_mis[, variable_index] <- current_column_data
  }
  
  # Get the indices for the introduced NA vals 
  mis_idx_df <- which(is.na(data_mis), arr.ind = TRUE)
  return(list(data_mis = data_mis, mis_idx_df = mis_idx_df))
}

# Define parameters
missing_var_props <- seq(0.1, 0.4, by = 0.1)  # Proportion of variables (columns) for missing data
num_repeats <- 20  # Number of repetitions per missing_var_prop
var_prop <- seq(.05, .4, .05)  # Proportion sequence to determine the cutoff

# Lists to store the generated datasets and indices for each experiment
datasets_MNAR <- list()
indices_MNAR <- list()

# Loop through each missing_var_prop and repeat the process 5 times
for (missing_var_prop in missing_var_props) {
  for (rep in 1:num_repeats) {
    # Generate MNAR dataset
    MNAR_result <- MNAR_generate(V5_unique_species, missing_var_prop = missing_var_prop, var_prop = var_prop)
    
    # Extract the data with missing values and the indices
    data_mis <- MNAR_result$data_mis
    mis_idx_df <- MNAR_result$mis_idx_df
    
    # Store the results in the lists, with names indicating the proportion and repetition
    datasets_MNAR[[paste0("MNAR_", missing_var_prop, "_rep_", rep)]] <- data_mis
    indices_MNAR[[paste0("MNAR_idx_", missing_var_prop, "_rep_", rep)]] <- mis_idx_df
  }
}

# Plotting some of the generated datasets (optional)
par(mfrow = c(2, 2))  # Set up the plot area for 2 rows and 2 columns

for (i in seq(1,20,5)) {
  name <- names(datasets_MNAR)[i]
  cat("Plotting missing map for: ", name, "\n")
  plot_missing_map(datasets_MNAR[[name]], name)
}

# Reset the plot area to 1 panel
par(mfrow = c(1, 1))
# Check if 30 percent still makes sense to do imputation
# NRMSE based on percentage of missing to identify the 'cutoff' where it still makes sense  
# (for example make a graph of NRMSE vs MISSING PROPORTION FOR EACH IMPUTATION METHOD )
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

# MNAR
# Define the missing proportions to plot for MNAR
proportions <- c(0.2, 0.4)

# Loop through the missing proportions
for (prop in proportions) {
  
  # Work only with the first repetition for each proportion (rep_1)
  dataset_name <- paste0("MNAR_", prop, "_rep_1")
  
  # Extract the corresponding dataset with missing values
  missing_data <- datasets_MNAR[[dataset_name]]
  
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
      labs(title = paste("Density Comparison: Complete vs MNAR Simulated Data (Proportion:", prop, ")"),
           x = "Value", y = "Density") +
      theme_minimal() +
      theme(legend.position = "right")
  )
}