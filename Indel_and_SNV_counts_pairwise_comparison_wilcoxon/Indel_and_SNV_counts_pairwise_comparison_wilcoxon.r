# MMRd Colorectal Neoplasia Statistical Analysis Script
# This script performs pairwise Wilcoxon tests for Indels.Exons and SNPs_exons.
# The script has been adjusted for public availability and general use.
# Author: Hagay Ladany
# Date: 4.Dec.2024
# ================================
# Load Required Libraries
# ================================
library(readxl)    # For reading Excel files
library(dplyr)      # For data manipulation
library(tidyr)      # For handling missing data
library(purrr)      # For functional programming
library(stringr)    # For string manipulation
# Function to perform pairwise Wilcoxon tests and adjust p-values
perform_pairwise_wilcox <- function(data, group_var, test_var) {
  groups <- unique(data[[group_var]])
  comparisons <- combn(groups, 2, simplify = FALSE)
  # Initialize an empty dataframe to store results
  results <- data.frame(group1 = character(),
                        group2 = character(),
                        p.value = numeric(),
                        stringsAsFactors = FALSE)
  # Perform Wilcoxon test for each pair
  for (comp in comparisons) {
    group1 <- comp[1]
    group2 <- comp[2]
    x1 <- data %>% filter((!!sym(group_var)) == group1) %>% pull(!!sym(test_var))
    x2 <- data %>% filter((!!sym(group_var)) == group2) %>% pull(!!sym(test_var))
    # Ensure both groups have data
    if(length(x1) > 0 & length(x2) > 0){
      wt <- wilcox.test(x1, x2)
      results <- rbind(results, data.frame(
        group1 = group1,
        group2 = group2,
        p.value = wt$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }
  # Adjust p-values using Benjamini-Hochberg (BH) method
  results <- results %>%
    mutate(p.adj = p.adjust(p.value, method = "BH"))
  return(results)
}
# ================================
# Load and Prepare Data
# ================================
# Define the path to your Excel file
metadata_file <- "Supplementary Table 1 - Metadata.xlsx"
# Read the Excel file
MetaData <- read_excel(metadata_file,na = "NA")
# Subset the data to include only rows where SNPs are not NA
no_COVTHRESH <- MetaData %>% filter(!is.na(SNPs))
# Replace NA in Lynch_germLine_mut with "Unknown" to ensure group names are complete
no_COVTHRESH <- no_COVTHRESH %>%
  mutate(Lynch_germLine_mut = replace_na(Lynch_germLine_mut, "Unknown"))
# Create 'group' column by combining 'class.', 'Histology', and 'Lynch_germLine_mut'
no_COVTHRESH <- no_COVTHRESH %>%
  mutate(group = paste(class., Histology, Lynch_germLine_mut, sep = "\n"))
# Clean up group names by trimming whitespace and removing any trailing newlines
no_COVTHRESH <- no_COVTHRESH %>%
  mutate(group = str_trim(group))
# Verify that the 'group' column has been created correctly
# table(no_COVTHRESH$group)
# ================================
# Perform Pairwise Wilcoxon Tests
# ================================
# Ensure 'group' is a factor
no_COVTHRESH$group <- as.factor(no_COVTHRESH$group)
# Perform pairwise Wilcoxon tests for Indels.Exons
results_indels <- perform_pairwise_wilcox(data = no_COVTHRESH, 
                                          group_var = "group", 
                                          test_var = "Indels.Exons")
# Perform pairwise Wilcoxon tests for SNPs_exons
results_snps <- perform_pairwise_wilcox(data = no_COVTHRESH, 
                                        group_var = "group", 
                                        test_var = "SNPs_exons")
# Save Indels.Exons pairwise comparison results
rez <- cbind(results_indels,results_snps[,c("p.value","p.adj")])
# Save SNPs.Exons pairwise comparison results
colnames(rez)[(ncol(rez)-3):ncol(rez)] <- c(c("indels_p.value","indels_p.adj","snv_p.value","snv_p.adj"))
cat("Pairwise Wilcoxon tests completed.\n")
# ================================
# End of Script
# ================================
