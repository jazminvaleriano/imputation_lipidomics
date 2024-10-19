source('Impute_wrapper.R')
library(dplyr)      # For data manipulation (group_by, summarise, etc.)
library(reshape2)   # For reshaping the data (melt)
library(ggplot2)    # For visualizations (boxplots, line plots)
library(patchwork)  # For combining plots
library(vegan)      # For PROCRUSTES analysis

##### IMPUTATION MCAR------------------------------------
# List of imputation functions
MCAR_imputation_methods <- list(
  Mean = Mean_wrapper,
  Median = Median_wrapper,
  RF = RF_wrapper,
  kNN = kNN_wrapper
)

# Empty list to store the imputed dataframes for each method and dataset
MCAR_imputed_dataframes <- list()

# Loop over each dataset in MCAR_datasets
for (dataset_name in names(mcar_datasets)) {
  # Extract the dataset with missing values
  missing_df <- mcar_datasets[[dataset_name]]
  
  # Create a list to store the imputed dataframes for this dataset
  dataset_imputed <- list()
  
  # Loop over each imputation method
  for (method_name in names(MCAR_imputation_methods)) {
    # Apply the imputation method and store the imputed dataframe
    imputed_df <- MCAR_imputation_methods[[method_name]](missing_df)
    
    # Store the imputed dataframe under the method and dataset
    dataset_imputed[[method_name]] <- imputed_df
  }
  print(paste0("Done processing ",dataset_name))
  # Store the imputed dataframes for this dataset
  MCAR_imputed_dataframes[[dataset_name]] <- dataset_imputed
}

saveRDS(MCAR_imputed_dataframes,file = 'MCAR_imputed_dataframes.rds')

##### IMPUTATION MNAR ----------------------------------
# List of imputation functions for MNAR
MNAR_imputation_methods <- list(
  Min_0.5 = Min_0.5_wrapper,
  RF = RF_wrapper,
  kNN = kNN_wrapper,
  QRILC = QRILC_wrapper
)

# Empty list to store the imputed dataframes for each method and dataset
MNAR_imputed_dataframes <- list()

# Loop over each dataset in MNAR_datasets
for (dataset_name in names(MNAR_datasets)) {
  # Extract the dataset with missing values
  missing_df <- MNAR_datasets[[dataset_name]]
  
  # Create a list to store the imputed dataframes for this dataset
  dataset_imputed <- list()
  print(paste("Processing dataset", dataset_name))
  
  # Loop over each imputation method
  for (method_name in names(MNAR_imputation_methods)) {
    # Apply the imputation method
    imputed_df <- MNAR_imputation_methods[[method_name]](missing_df)
    
    # Store the imputed dataframe under the method and dataset
    dataset_imputed[[method_name]] <- imputed_df
  }
  
  # Store the imputed dataframes for this dataset
  MNAR_imputed_dataframes[[dataset_name]] <- dataset_imputed
}

saveRDS(MNAR_imputed_dataframes,file = 'MNAR_imputed_dataframes_v2.rds')

##### EVALUATION of MCAR imputation methods (weighted NRMSE) --------------
# List of imputation methods
MCAR_imputed_dataframes <- readRDS('MCAR_imputed_dataframes.rds')
V5_unique_species <- readRDS('V5_unique_species.rds')
MCAR_imputation_methods <- list("Mean", "Median", "RF", "kNN")
mcar_dataframes <- readRDS('mcar_datasets.rds')

# Empty list to store NRMSE results for each dataset and method
MCAR_nrmse_results <- list()

# Loop over each dataset in MCAR_datasets
for (dataset_name in names(mcar_dataframes)) {
  # Extract the dataset with missing values
  missing_df <-mcar_dataframes[[dataset_name]]
  
  # The full dataset is always V5_unique_species
  full_df <- V5_unique_species
  
  # Create a list to store the NRMSE results for this dataset
  dataset_nrmse <- list()
  
  # Loop over each imputation method (including kNN)
  for (method_name in MCAR_imputation_methods) {
    # Retrieve the pre-imputed dataframe
    imputed_df <- MCAR_imputed_dataframes[[dataset_name]][[method_name]]
    
    # Calculate the Weighted NRMSE for each column in this imputation method
    nrmse_values <- weighted_nrmse_by_column(imputed_df, full_df, missing_df)
    
    # Store Weighted NRMSE for each column under the method and dataset
    dataset_nrmse[[method_name]] <- nrmse_values
  }
  
  # Store the Weighted NRMSE results for this dataset
  MCAR_nrmse_results[[dataset_name]] <- dataset_nrmse
}

# Convert the Weighted NRMSE results into a long format for plotting
mcar_nrmse_df <- do.call(rbind, lapply(names(MCAR_nrmse_results), function(dataset_name) {
  # Create a data frame for each dataset, including the Weighted NRMSE values per column for each imputation method
  data.frame(
    Dataset = dataset_name,
    Column = rep(1:ncol(V5_unique_species), length(MCAR_imputation_methods)),
    kNN = MCAR_nrmse_results[[dataset_name]]$kNN,
    Mean = MCAR_nrmse_results[[dataset_name]]$Mean,
    Median = MCAR_nrmse_results[[dataset_name]]$Median,
    RF = MCAR_nrmse_results[[dataset_name]]$RF
  )
}))

# Extract the MCAR proportion from the dataset names (e.g., "MCAR_0.1_rep_1" -> 0.1)
mcar_nrmse_df$MCAR_Prop <- as.numeric(gsub("MCAR_([0-9.]+).*", "\\1", mcar_nrmse_df$Dataset))

# Reshape the data for plotting (long format), keeping the column information
mcar_nrmse_melted <- reshape2::melt(mcar_nrmse_df, id.vars = c("Dataset", "MCAR_Prop", "Column"), 
                                    variable.name = "Imputation_Method", value.name = "NRMSE")

# Plot boxplots showing the variance across columns for each method and MCAR proportion
ggplot(mcar_nrmse_melted, aes(x = as.factor(MCAR_Prop), y = NRMSE, fill = Imputation_Method)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +  # Boxplot with transparency to show variance across columns
  labs(title = "Weighted NRMSE Boxplots Across Columns for Imputation Methods (MCAR)",
       x = "MCAR Proportion (%)", y = "Weighted NRMSE", fill = "Imputation Method") +
  theme_minimal() +
  theme(legend.position = "right")

# Calculate mean Weighted NRMSE for each method and MCAR proportion (to plot the trend lines)
mcar_nrmse_means <- mcar_nrmse_melted %>%
  dplyr::group_by(MCAR_Prop, Imputation_Method) %>%
  dplyr::summarise(mean_nrmse = mean(NRMSE, na.rm = TRUE))

# Create the combined plot (boxplots + mean trend lines)
ggplot(mcar_nrmse_melted, aes(x = as.factor(MCAR_Prop), y = NRMSE, fill = Imputation_Method)) +
  # Create the boxplot to show variance across columns
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +  
  # Overlay with a line plot showing the mean trends
  geom_line(data = mcar_nrmse_means, aes(x = as.factor(MCAR_Prop), y = mean_nrmse, 
                                         group = Imputation_Method, color = Imputation_Method), size = 1) +
  geom_point(data = mcar_nrmse_means, aes(x = as.factor(MCAR_Prop), y = mean_nrmse, 
                                          color = Imputation_Method), size = 2) +
  labs(title = "Weighted NRMSE Boxplots and Mean Trends for Imputation Methods (MCAR)",
       x = "MCAR Proportion (%)", y = "Weighted NRMSE", fill = "Imputation Method", color = "Imputation Method") +
  theme_minimal() +
  theme(legend.position = "right")

# 1. Boxplots with facet wrap (one panel per method)
boxplot_facet <- ggplot(mcar_nrmse_melted, aes(x = as.factor(MCAR_Prop), y = NRMSE, fill = Imputation_Method)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  facet_wrap(~ Imputation_Method, scales = "free") +  # Facet by method
  ylim(0, 0.5) +  # Set y-axis limits for all facets
  labs(title = "Weighted NRMSE Boxplots for Each Imputation Method (MCAR)",
       x = "MCAR Proportion (%)", y = "Weighted NRMSE") +
  theme_minimal() +
  theme(legend.position = "none")  # No legend for the boxplot facet

# 2. Line plot for mean Weighted NRMSE across all methods
line_plot_means <- ggplot(mcar_nrmse_means, aes(x = as.factor(MCAR_Prop), y = mean_nrmse, color = Imputation_Method, group = Imputation_Method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Mean Weighted NRMSE for All Imputation Methods",
       x = "MCAR Proportion (%)", y = "Mean Weighted NRMSE", color = "Imputation Method") +
  theme_minimal() +
  theme(legend.position = "right")

# Combine the plots using patchwork
final_plot <- boxplot_facet / line_plot_means

# Display the combined plot
print(final_plot)


##### EVALUATION of MNAR imputation methods (weighted NRMSE)------
# List of imputation functions (names only)
MNAR_imputation_methods <- list("Min_0.5","RF","kNN","QRILC")
V5_unique_species <- readRDS('V5_unique_species.rds')
MNAR_datasets <- readRDS('MNAR_datasets.rds')
MNAR_imputed_dataframes <- readRDS('MNAR_imputed_dataframes_v2.rds')

# Empty list to store NRMSE results for each dataset and method
MNAR_nrmse_results <- list()

# Loop over each dataset in MNAR_datasets
for (dataset_name in names(MNAR_datasets)) {
  # Extract the dataset with missing values
  missing_df <- MNAR_datasets[[dataset_name]]
  
  # The full dataset is always V5_unique_species
  full_df <- V5_unique_species
  
  # Create a list to store the NRMSE results for this dataset
  dataset_nrmse <- list()
  print(paste("Processing dataset", dataset_name))
  
  # Loop over each imputation method (using pre-imputed data)
  for (method_name in MNAR_imputation_methods) {
    # Retrieve the pre-imputed dataframe
    imputed_df <- MNAR_imputed_dataframes[[dataset_name]][[method_name]]
    
    # Calculate the weighted NRMSE for each column in this imputation method
    weighted_nrmse_values <- weighted_nrmse_by_column(imputed_df, full_df, missing_df)
    
    # Store weighted NRMSE for each column under the method and dataset
    dataset_nrmse[[method_name]] <- weighted_nrmse_values
  }
  
  # Store the weighted NRMSE results for this dataset
  MNAR_nrmse_results[[dataset_name]] <- dataset_nrmse
}

# Convert the Weighted NRMSE results into a long format for plotting
mnar_nrmse_df <- do.call(rbind, lapply(names(MNAR_nrmse_results), function(dataset_name) {
  # Create a data frame for each dataset, including the Weighted NRMSE values per column for each imputation method
  data.frame(
    Dataset = dataset_name,
    Column = rep(1:ncol(V5_unique_species), length(MNAR_imputation_methods)),
    kNN = MNAR_nrmse_results[[dataset_name]]$kNN,
    Min_0.5 = MNAR_nrmse_results[[dataset_name]]$Min_0.5,
    QRILC = MNAR_nrmse_results[[dataset_name]]$QRILC,
    RF = MNAR_nrmse_results[[dataset_name]]$RF
  )
}))

# Extract the MNAR proportion from the dataset names (e.g., "MNAR_0.1_rep_1" -> 0.1)
mnar_nrmse_df$MNAR_Prop <- as.numeric(gsub("MNAR_var_prop_([0-9.]+).*", "\\1", mnar_nrmse_df$Dataset))

# Reshape the data for plotting (long format), keeping the column information
mnar_nrmse_melted <- reshape2::melt(mnar_nrmse_df, id.vars = c("Dataset", "MNAR_Prop", "Column"), 
                                    variable.name = "Imputation_Method", value.name = "Weighted_NRMSE")

# Plot boxplots showing the variance across columns for each method and MNAR proportion
ggplot(mnar_nrmse_melted, aes(x = as.factor(MNAR_Prop), y = Weighted_NRMSE, fill = Imputation_Method)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +  # Boxplot with transparency to show variance across columns
  labs(title = "Weighted NRMSE Boxplots Across Columns for Imputation Methods (MNAR)",
       x = "MNAR Proportion (%)", y = "Weighted NRMSE", fill = "Imputation Method") +
  theme_minimal() +
  theme(legend.position = "right")

# Calculate mean Weighted NRMSE for each method and MNAR proportion (to plot the trend lines)
mnar_nrmse_means <- mnar_nrmse_melted %>%
  dplyr::group_by(MNAR_Prop, Imputation_Method) %>%
  dplyr::summarise(mean_nrmse = mean(Weighted_NRMSE, na.rm = TRUE))

# Create the combined plot (boxplots + mean trend lines)
ggplot(mnar_nrmse_melted, aes(x = as.factor(MNAR_Prop), y = Weighted_NRMSE, fill = Imputation_Method)) +
  # Create the boxplot to show variance across columns
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +  
  # Overlay with a line plot showing the mean trends
  geom_line(data = mnar_nrmse_means, aes(x = as.factor(MNAR_Prop), y = mean_nrmse, 
                                         group = Imputation_Method, color = Imputation_Method), size = 1) +
  geom_point(data = mnar_nrmse_means, aes(x = as.factor(MNAR_Prop), y = mean_nrmse, 
                                          color = Imputation_Method), size = 2) +
  labs(title = "Weighted NRMSE Boxplots and Mean Trends for Imputation Methods (MNAR)",
       x = "MNAR Proportion (%)", y = "Weighted NRMSE", fill = "Imputation Method", color = "Imputation Method") +
  theme_minimal() +
  ylim(0, 0.6) +
  theme(legend.position = "right")

# 1. Boxplots with facet wrap (one panel per method)
boxplot_facet_mnar <- ggplot(mnar_nrmse_melted, aes(x = as.factor(MNAR_Prop), y = Weighted_NRMSE, fill = Imputation_Method)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  facet_wrap(~ Imputation_Method, scales = "free") +  # Facet by method
  ylim(0, 0.5) +  # Set y-axis limits for all facets (adjust as needed for MNAR)
  labs(title = "Weighted NRMSE Boxplots for Each Imputation Method (MNAR)",
       x = "MNAR Proportion (%)", y = "Weighted NRMSE") +
  theme_minimal() +
  theme(legend.position = "none")  # No legend for the boxplot facet

# 2. Line plot for mean Weighted NRMSE across all methods
line_plot_means_mnar <- ggplot(mnar_nrmse_means, aes(x = as.factor(MNAR_Prop), y = mean_nrmse, color = Imputation_Method, group = Imputation_Method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Mean Weighted NRMSE for All Imputation Methods (MNAR)",
       x = "MNAR Proportion (%)", y = "Mean Weighted NRMSE", color = "Imputation Method") +
  theme_minimal() +
  ylim(0, 0.5) +
  theme(legend.position = "right")

# Combine the plots using patchwork
final_plot_mnar <- boxplot_facet_mnar / line_plot_means_mnar

# Display the combined MNAR plot
print(final_plot_mnar)

##### PROCRUSTES FOR FULL DS----------------------
# Initialize a data frame to store Procrustes SS and missing proportions
procrustes_results <- data.frame()

# Loop through each dataset in MNAR_datasets
for (i in seq_along(MNAR_datasets)) {
  # Extract the current dataset with missing values and the full dataset
  missing_df <- MNAR_datasets[[i]]
  full_df <- V5_unique_species  # Assuming V5_unique_species is the full dataset
  
  # Get the proportion of missing data in this dataset
  miss_prop <- mean(is.na(missing_df))
  
  # Loop through each imputation method
  for (method_name in names(MNAR_imputed_dataframes[[i]])) {
    # Get the imputed dataset for this method
    imputed_df <- MNAR_imputed_dataframes[[i]][[method_name]]
    
    # Perform PCA on both the full dataset and the imputed dataset
    full_pca <- prcomp(full_df, center = TRUE, scale. = TRUE)
    imputed_pca <- prcomp(imputed_df, center = TRUE, scale. = TRUE)
    
    # Use the first two principal components for Procrustes analysis
    full_pca_reduced <- full_pca$x[, 1:2]      # First 2 principal components from full data
    imputed_pca_reduced <- imputed_pca$x[, 1:2]  # First 2 principal components from imputed data
    
    # Perform symmetric Procrustes analysis (allow both datasets to transform)
    procrustes_result <- procrustes(full_pca_reduced, imputed_pca_reduced, symmetric = TRUE)
    
    # Calculate the Procrustes sum of squares (SS)
    procrustes_ss <- procrustes_result$ss
    
    # Append the results (method, missing proportion, and SS) to the results data frame
    procrustes_results <- rbind(procrustes_results, 
                                data.frame(Method = method_name, 
                                           Miss_Prop = miss_prop, 
                                           Pro_SS = procrustes_ss))
  }
}

# Now that we have collected the Procrustes SS values, let's reshape and plot the results

# Plot the Procrustes SS for each method against the missing proportion
ggplot(procrustes_results, aes(x = Miss_Prop, y = Pro_SS, color = Method, group = Method)) +
  geom_point() +
  geom_line() +
  labs(title = "Procrustes SS vs Missing Proportion for Each Imputation Method",
       x = "Missing Proportion",
       y = "Procrustes Sum of Squares (SS)",
       color = "Imputation Method") +
  theme_minimal()

##### PROCRUSTES TEST FOR MCAR ------------------
# Initialize a data frame to store Procrustes SS and missing proportions for MCAR data
procrustes_ss_results_MCAR <- data.frame()
MCAR_datasets <- mcar_dataframes
# Loop through each dataset in MCAR_datasets
for (i in seq_along(MCAR_datasets)) {
  # Extract the current dataset with missing values and the full dataset
  missing_df <- MCAR_datasets[[i]]
  full_df <- V5_unique_species  # Assuming V5_unique_species is the full dataset
  
  # Get the proportion of missing data in this dataset
  miss_prop <- mean(is.na(missing_df))
  
  # Loop through each imputation method
  for (method_name in names(MCAR_imputed_dataframes[[i]])) {
    # Get the imputed dataset for this method
    imputed_df <- MCAR_imputed_dataframes[[i]][[method_name]]
    
    # Perform PCA on both the full dataset and the imputed dataset
    full_pca <- prcomp(full_df, center = TRUE, scale. = TRUE)
    imputed_pca <- prcomp(imputed_df, center = TRUE, scale. = TRUE)
    
    # Use the first two principal components for Procrustes analysis
    full_pca_reduced <- full_pca$x[, 1:2]      # First 2 principal components from full data
    imputed_pca_reduced <- imputed_pca$x[, 1:2]  # First 2 principal components from imputed data
    
    # Perform symmetric Procrustes analysis (allow both datasets to transform)
    procrustes_result <- procrustes(full_pca_reduced, imputed_pca_reduced, symmetric = TRUE)
    
    # Calculate the Procrustes sum of squares (SS)
    procrustes_ss <- procrustes_result$ss
    
    # Append the results (method, missing proportion, and SS) to the results data frame
    procrustes_ss_results_MCAR <- rbind(procrustes_ss_results_MCAR, 
                                        data.frame(Method = method_name, 
                                                   Miss_Prop = miss_prop, 
                                                   Pro_SS = procrustes_ss))
  }
}

# Now that we have collected the Procrustes SS values, let's reshape and plot the results
# Summarize data: Calculate mean and standard error (SE) for each method and missing proportion
procrustes_summary <- procrustes_ss_results_MCAR %>%
  group_by(Miss_Prop, Method) %>%
  summarise(mean_SS = mean(Pro_SS),
            se_SS = sd(Pro_SS) / sqrt(n()))

# Plot the mean Procrustes SS with error bars (SE)
ggplot(procrustes_summary, aes(x = Miss_Prop, y = mean_SS, color = Method, group = Method)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_SS - se_SS, ymax = mean_SS + se_SS), width = 0.02) +
  labs(title = "Procrustes SS vs Missing Proportion (Mean ± SE)",
       x = "Missing Proportion",
       y = "Mean Procrustes Sum of Squares (SS)",
       color = "Imputation Method") +
  theme_minimal()



##### PROCRUSTES FOR FULL MANR ----------------------
# Initialize a data frame to store Procrustes SS and missing proportions for MNAR data
procrustes_ss_results_MNAR <- data.frame()

# Loop through each dataset in MNAR_datasets
for (i in seq_along(MNAR_datasets)) {
  # Extract the current dataset with missing values and the full dataset
  missing_df <- MNAR_datasets[[i]]
  full_df <- V5_unique_species  # Assuming V5_unique_species is the full dataset
  
  # Get the proportion of missing data in this dataset
  miss_prop <- mean(is.na(missing_df))
  
  # Loop through each imputation method
  for (method_name in names(MNAR_imputed_dataframes[[i]])) {
    # Get the imputed dataset for this method
    imputed_df <- MNAR_imputed_dataframes[[i]][[method_name]]
    
    # Perform PCA on both the full dataset and the imputed dataset
    full_pca <- prcomp(full_df, center = TRUE, scale. = TRUE)
    imputed_pca <- prcomp(imputed_df, center = TRUE, scale. = TRUE)
    
    # Use the first two principal components for Procrustes analysis
    full_pca_reduced <- full_pca$x[, 1:2]      # First 2 principal components from full data
    imputed_pca_reduced <- imputed_pca$x[, 1:2]  # First 2 principal components from imputed data
    
    # Perform symmetric Procrustes analysis (allow both datasets to transform)
    procrustes_result <- procrustes(full_pca_reduced, imputed_pca_reduced, symmetric = TRUE)
    
    # Calculate the Procrustes sum of squares (SS)
    procrustes_ss <- procrustes_result$ss
    
    # Append the results (method, missing proportion, and SS) to the results data frame
    procrustes_ss_results_MNAR <- rbind(procrustes_ss_results_MNAR, 
                                        data.frame(Method = method_name, 
                                                   Miss_Prop = miss_prop, 
                                                   Pro_SS = procrustes_ss))
  }
}

# Now that we have collected the Procrustes SS values, we reshape and plot the results

# Plot the Procrustes SS for each method against the missing proportion
ggplot(procrustes_ss_results_MNAR, aes(x = Miss_Prop, y = Pro_SS, color = Method, group = Method)) +
  geom_point() +
  geom_line() +
  labs(title = "Procrustes SS vs Missing Proportion for Each Imputation Method (MNAR)",
       x = "Missing Proportion",
       y = "Procrustes Sum of Squares (SS)",
       color = "Imputation Method") +
  theme_minimal()



##### PROCRUSTES FOR 1 METHOD+PROPORTION--------
# Load the vegan package for Procrustes analysis
perform_procrustes_analysis <- function(missing_df, full_df, imputed_df, method_name = "Imputation Method", plot_result = TRUE) {
  # Load the necessary library
  library(vegan)
  
  # Perform PCA on both the full dataset and the imputed dataset
  full_pca <- prcomp(full_df, center = TRUE, scale. = TRUE)
  imputed_pca <- prcomp(imputed_df, center = TRUE, scale. = TRUE)
  
  # Use the first two principal components for Procrustes analysis
  full_pca_reduced <- full_pca$x[, 1:2]      # First 2 principal components from full data
  imputed_pca_reduced <- imputed_pca$x[, 1:2]  # First 2 principal components from imputed data
  
  # Perform symmetric Procrustes analysis (allow both datasets to transform)
  procrustes_result <- procrustes(full_pca_reduced, imputed_pca_reduced, symmetric = TRUE)
  
  # Calculate the Procrustes sum of squares (SS)
  procrustes_ss <- procrustes_result$ss
  
  # Output the Procrustes sum of squares (SS)
  cat("Procrustes SS for", method_name, ":", procrustes_ss, "\n")
  
  # Plot the Procrustes result to visualize the alignment
  if (plot_result) {
    plot(procrustes_result, main = paste("Procrustes errors for", method_name))
  }
  
  # Return the Procrustes SS
  return(procrustes_ss)
}

# MNAR HM METHODS (QRILC and HM) at 0.1:
missing_df <- MNAR_datasets[[2]]  # Specify which dataset
full_df <- V5_unique_species      # Full dataset
perform_procrustes_analysis(missing_df, full_df, 
                            imputed_df = MNAR_imputed_dataframes[[2]][["Min_0.5"]],
                            method_name = "Min_0.5", 
                            plot_result = TRUE)

perform_procrustes_analysis(missing_df, full_df, 
                            imputed_df = MNAR_imputed_dataframes[[2]][["QRILC"]],
                            method_name = "QRILC", 
                            plot_result = TRUE)

# MNAR HM METHODS (QRILC and HM) at 0.4:
missing_df <- MNAR_datasets[[5]]  # Specify which dataset
full_df <- V5_unique_species      # Full dataset
perform_procrustes_analysis(missing_df, full_df, 
                            imputed_df = MNAR_imputed_dataframes[[2]][["Min_0.5"]],
                            method_name = "Min_0.5", 
                            plot_result = TRUE)

perform_procrustes_analysis(missing_df, full_df, 
                            imputed_df = MNAR_imputed_dataframes[[2]][["QRILC"]],
                            method_name = "QRILC", 
                            plot_result = TRUE)

# MCAR HM METHODS (QRILC and HM) at 0.15:
missing_df <- MCAR_datasets[[13]]  # Specify which dataset
full_df <- V5_unique_species      # Full dataset
perform_procrustes_analysis(missing_df, full_df, 
                            imputed_df = MCAR_imputed_dataframes[[13]][["Mean"]],
                            method_name = "Mean", 
                            plot_result = TRUE)
perform_procrustes_analysis(missing_df, full_df, 
                            imputed_df = MCAR_imputed_dataframes[[13]][["RF"]],
                            method_name = "RF", 
                            plot_result = TRUE)

# MCAR HM METHODS (QRILC and HM) at 0.45:
missing_df <- MCAR_datasets[[43]]  # Specify which dataset
full_df <- V5_unique_species      # Full dataset
perform_procrustes_analysis(missing_df, full_df, 
                            imputed_df = MCAR_imputed_dataframes[[43]][["Mean"]],
                            method_name = "Mean", 
                            plot_result = TRUE)
perform_procrustes_analysis(missing_df, full_df, 
                            imputed_df = MCAR_imputed_dataframes[[43]][["RF"]],
                            method_name = "RF", 
                            plot_result = TRUE)

