# Load required packages
library(car)  # For Levene’s test
library(ggplot2)

# ANOVA for MNAR------------------------------------------------
# Step 1: Fit a two-way ANOVA model (with interaction)
anova_model <- aov(sqrt(Weighted_NRMSE) ~ Imputation_Method * MNAR_Prop, data = mnar_nrmse_melted)

# Step 2: Check residuals for normality
# 2.1: Histogram of residuals
ggplot(data.frame(residuals = residuals(anova_model)), aes(x = residuals)) +
  geom_histogram(binwidth = 0.05, fill = "blue", alpha = 0.7) +
  labs(title = "Histogram of Residuals", x = "Residuals", y = "Frequency") +
  theme_minimal()

par(mfrow=c(1,1))
# 2.2: Q-Q plot of residuals
qqnorm(residuals(anova_model))
qqline(residuals(anova_model), col = "red")

plot(fitted(anova_model), resid(anova_model), main="Tukey-anscombe")
qqPlot(resid(anova_model), dist="norm", mean=mean(resid(anova_model)), sd=sd(resid(anova_model)), 
       main="QQ Plot residues")

# Step 4: Display the ANOVA summary
summary(anova_model)

# Perform Kruskal-Wallis test
kruskal_test_result <- kruskal.test(Weighted_NRMSE ~ Imputation_Method, data = mnar_nrmse_melted)

# Print the results
print(kruskal_test_result)

# Pairwise test for Weighted NRMSE, MNAR ------------

# Subset the data to include only the columns you need (Imputation_Method, NRMSE, MNAR_Prop)

# Perform the pairwise Wilcoxon Rank-Sum Test
pairwise_wilcox_test <- pairwise.wilcox.test(mnar_nrmse_melted$Weighted_NRMSE, mnar_nrmse_melted$Imputation_Method, 
                                             p.adjust.method = "bonferroni")
# Print the results
print(pairwise_wilcox_test)

# Loop through each MNAR proportion and perform the pairwise Wilcoxon test
for (prop in unique(mnar_nrmse_melted$MNAR_Prop)) {
  # Subset the data for the current MNAR proportion
  subset_data <- mnar_nrmse_melted[mnar_nrmse_melted$MNAR_Prop == prop, ]
  
  # Perform the pairwise Wilcoxon Rank-Sum Test for this MNAR proportion
  wilcox_test <- pairwise.wilcox.test(subset_data$Weighted_NRMSE, subset_data$Imputation_Method, 
                                      p.adjust.method = "bonferroni")
  
  # Print results for each MNAR proportion
  cat("Wilcoxon Rank-Sum Test for MNAR Proportion:", prop, "\n")
  print(wilcox_test)
}

# Step 1: Fit a two-way ANOVA model (with interaction)
anova_model <- aov( log(Pro_SS) ~ Method * Miss_Prop, data = procrustes_ss_results_MNAR)

# Step 2: Check residuals for normality
# 2.1: Histogram of residuals
ggplot(data.frame(residuals = residuals(anova_model)), aes(x = residuals)) +
  geom_histogram(binwidth = 0.05, fill = "blue", alpha = 0.7) +
  labs(title = "Histogram of Residuals", x = "Residuals", y = "Frequency") +
  theme_minimal()

par(mfrow=c(1,1))
# 2.2: Q-Q plot of residuals
qqnorm(residuals(anova_model))
qqline(residuals(anova_model), col = "red")

plot(fitted(anova_model), resid(anova_model), main="Tukey-anscombe")
qqPlot(resid(anova_model), dist="norm", mean=mean(resid(anova_model)), sd=sd(resid(anova_model)), 
       main="QQ Plot residues")

# Step 4: Display the ANOVA summary
summary(anova_model)



# ANOVA for MCAR ------------------------------------------------
# Step 1: Fit a two-way ANOVA model (with interaction)
anova_model_mcar <- aov(log(NRMSE) ~ Imputation_Method * MCAR_Prop, data = mcar_nrmse_melted)

# Step 2: Check residuals for normality
# 2.1: Histogram of residuals
ggplot(data.frame(residuals = residuals(anova_model_mcar)), aes(x = residuals)) +
  geom_histogram(binwidth = 0.05, fill = "blue", alpha = 0.7) +
  labs(title = "Histogram of Residuals", x = "Residuals", y = "Frequency") +
  theme_minimal()

par(mfrow=c(1,1))
# 2.2: Q-Q plot of residuals
qqnorm(residuals(anova_model_mcar))
qqline(residuals(anova_model_mcar), col = "red")

# Tukey-Anscombe plot for checking residuals vs fitted values
plot(fitted(anova_model_mcar), resid(anova_model_mcar), main="Tukey-Anscombe Plot")

# Optional: Q-Q plot for residuals using `car` package
library(car)
qqPlot(resid(anova_model_mcar), dist="norm", mean=mean(resid(anova_model_mcar)), sd=sd(resid(anova_model_mcar)), 
       main="QQ Plot of Residuals")

# Step 3: Display the ANOVA summary
summary(anova_model_mcar)

# Perform Kruskal-Wallis test for MCAR ------------------------------------------------
kruskal_test_mcar <- kruskal.test(NRMSE ~ Imputation_Method, data = mcar_nrmse_melted)

# Print the Kruskal-Wallis results
print(kruskal_test_mcar)

# Pairwise test for Weighted NRMSE, MCAR ------------

# Perform the pairwise Wilcoxon Rank-Sum Test for MCAR
pairwise_wilcox_test_mcar <- pairwise.wilcox.test(mcar_nrmse_melted$NRMSE, mcar_nrmse_melted$Imputation_Method, 
                                                  p.adjust.method = "bonferroni")
# Print the pairwise Wilcoxon results
print(pairwise_wilcox_test_mcar)

# Loop through each MCAR proportion and perform the pairwise Wilcoxon test
for (prop in unique(mcar_nrmse_melted$MCAR_Prop)) {
  # Subset the data for the current MCAR proportion
  subset_data_mcar <- mcar_nrmse_melted[mcar_nrmse_melted$MCAR_Prop == prop, ]
  
  # Perform the pairwise Wilcoxon Rank-Sum Test for this MCAR proportion
  wilcox_test_mcar <- pairwise.wilcox.test(subset_data_mcar$NRMSE, subset_data_mcar$Imputation_Method, 
                                           p.adjust.method = "bonferroni")
  
  # Print results for each MCAR proportion
  cat("Wilcoxon Rank-Sum Test for MCAR Proportion:", prop, "\n")
  print(wilcox_test_mcar)
}

# ANOVA for Procrustes SS (MCAR) ------------------------------------------------

# Step 1: Fit a two-way ANOVA model (with interaction)
anova_model_procrustes_mcar <- aov(log(Pro_SS) ~ Method * Miss_Prop, data = procrustes_ss_results_MCAR)

# Step 2: Check residuals for normality
# 2.1: Histogram of residuals
ggplot(data.frame(residuals = residuals(anova_model_procrustes_mcar)), aes(x = residuals)) +
  geom_histogram(binwidth = 0.05, fill = "blue", alpha = 0.7) +
  labs(title = "Histogram of Residuals", x = "Residuals", y = "Frequency") +
  theme_minimal()

par(mfrow=c(1,1))
# 2.2: Q-Q plot of residuals
qqnorm(residuals(anova_model_procrustes_mcar))
qqline(residuals(anova_model_procrustes_mcar), col = "red")

# Tukey-Anscombe plot for residuals vs fitted values
plot(fitted(anova_model_procrustes_mcar), resid(anova_model_procrustes_mcar), main="Tukey-Anscombe Plot")

# Q-Q plot for residuals using `car` package
qqPlot(resid(anova_model_procrustes_mcar), dist="norm", mean=mean(resid(anova_model_procrustes_mcar)), sd=sd(resid(anova_model_procrustes_mcar)), 
       main="QQ Plot of Residuals")

# Step 4: Display the ANOVA summary
summary(anova_model_procrustes_mcar)

# Wilcoxon Rank-Sum Test for MCAR PRO_SS
# Loop through each MCAR proportion and perform the pairwise Wilcoxon test
for (prop in unique(procrustes_ss_results_MCAR$Miss_Prop)) {
  # Subset the data for the current MCAR proportion
  subset_data_mcar <- procrustes_ss_results_MCAR[procrustes_ss_results_MCAR$Miss_Prop == prop, ]
  
  # Perform the pairwise Wilcoxon Rank-Sum Test for this MCAR proportion
  wilcox_test_mcar <- pairwise.wilcox.test(subset_data_mcar$Pro_SS, subset_data_mcar$Method, 
                                           p.adjust.method = "bonferroni")
  
  # Print results for each MCAR proportion
  cat("Wilcoxon Rank-Sum Test for MCAR Proportion:", prop, "\n")
  print(wilcox_test_mcar)
}

