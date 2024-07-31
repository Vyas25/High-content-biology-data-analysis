
#### Phospho related sites among all phases for PI3K 

# Load necessary libraries
library(reshape2)
library(dplyr)

# Load the data
data <- read.csv("C:/Users/avani/Downloads/Updated_Data_Mutated-June24.csv")

# Define PI3K signaling genes of interest
pi3k_signaling_genes <- c("Pik3ca", "Pik3cb", "Pik3cd", "Akt1", "Akt2", "Mtor")

# Define all phases of interest
phases_of_interest <- c("phase_earlyG1_A", "phase_earlyG1_B", "phase_earlyG1_C", 
                        "phase_lateG1_A", "phase_lateG1_B", "phase_lateG1_C",
                        "phase_S_A", "phase_S_B", "phase_S_C",
                        "phase_G2_A", "phase_G2_B", "phase_G2_C",
                        "phase_M_A", "phase_M_B", "phase_M_C")

# Melt the data to long format
phospho_data_long <- melt(data, 
                          id.vars = c("Gene", "New_Gene_Residue", "PhosR_Rownames"), 
                          measure.vars = phases_of_interest,
                          variable.name = "Phase", value.name = "Phosphorylation")

# Create a new column to categorize phases
phospho_data_long$PhaseCategory <- ifelse(grepl("earlyG1", phospho_data_long$Phase), "eG1",
                                          ifelse(grepl("lateG1", phospho_data_long$Phase), "LG1",
                                                 ifelse(grepl("S", phospho_data_long$Phase), "S",
                                                        ifelse(grepl("G2", phospho_data_long$Phase), "G2",
                                                               ifelse(grepl("M", phospho_data_long$Phase), "M", NA)))))

# Handle missing values by imputing with column means
phospho_data_long$Phosphorylation[is.na(phospho_data_long$Phosphorylation)] <- 
  ave(phospho_data_long$Phosphorylation, phospho_data_long$Phase, FUN = function(x) mean(x, na.rm = TRUE))

# Check the PhaseCategory column and filter out NA values
phospho_data_long <- phospho_data_long[!is.na(phospho_data_long$PhaseCategory), ]

# Function to print phosphorylation sites and their levels for each gene
print_gene_phospho_sites <- function(gene) {
  gene_data <- phospho_data_long[phospho_data_long$Gene == gene, ]
  gene_sites <- unique(gene_data$New_Gene_Residue)
  
  cat("Number of unique phosphorylation sites for", gene, ":", length(gene_sites), "\n")
  cat("Phosphorylation sites for", gene, ":\n")
  print(gene_sites)
  
  cat("Phosphorylation levels by phase for", gene, ":\n")
  print(gene_data)
  cat("\n\n")
}

# Print phosphorylation sites and their levels for each gene
for (gene in pi3k_signaling_genes) {
  print_gene_phospho_sites(gene)
}

#### phospho related sites among all phases for ERK

## Load necessary libraries
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# Load the data
data <- read.csv("C:/Users/avani/Downloads/Updated_Data_Mutated-June24.csv")

# Define ERK signaling genes
erk_signaling_genes <- c("Mapk1", "Mapk3", "Raf1", "Braf", "Map2k1", "Map2k2")

# Extract phosphorylation data for ERK signaling genes
phospho_data_erk <- data[data$Gene %in% erk_signaling_genes, ]

# Check if data is extracted properly
print(phospho_data_erk)

# If the data frame is empty, stop here and print a message
if (nrow(phospho_data_erk) == 0) {
  stop("No data found for the specified ERK signaling genes.")
}

# Define all phases of interest
phases_of_interest <- c("phase_earlyG1_A", "phase_earlyG1_B", "phase_earlyG1_C", 
                        "phase_lateG1_A", "phase_lateG1_B", "phase_lateG1_C",
                        "phase_S_A", "phase_S_B", "phase_S_C",
                        "phase_G2_A", "phase_G2_B", "phase_G2_C",
                        "phase_M_A", "phase_M_B", "phase_M_C")

# Prepare the data for heatmap
phospho_data_wide_erk <- phospho_data_erk[, c("Gene", "New_Gene_Residue", phases_of_interest)]
phospho_matrix_erk <- as.matrix(phospho_data_wide_erk[,-(1:2)])
rownames(phospho_matrix_erk) <- paste(phospho_data_wide_erk$Gene, phospho_data_wide_erk$New_Gene_Residue, sep = "_")

# Create a heatmap for ERK signaling genes
Heatmap(
  phospho_matrix_erk, 
  name = "Phosphorylation",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "Cell Cycle Phase",
  row_title = "Gene and Residue",
  color = colorRamp2(c(min(phospho_matrix_erk, na.rm = TRUE), 0, max(phospho_matrix_erk, na.rm = TRUE)), c("blue", "white", "red"))
)

# Prepare the data for plotting
phospho_data_long_erk <- melt(phospho_data_erk, 
                              id.vars = c("Gene", "New_Gene_Residue"), 
                              measure.vars = phases_of_interest,
                              variable.name = "Phase", value.name = "Phospho_Values")

# Create a new column to categorize phases
phospho_data_long_erk$PhaseCategory <- ifelse(grepl("earlyG1", phospho_data_long_erk$Phase), "earlyG1",
                                              ifelse(grepl("lateG1", phospho_data_long_erk$Phase), "lateG1",
                                                     ifelse(grepl("S", phospho_data_long_erk$Phase), "S",
                                                            ifelse(grepl("G2", phospho_data_long_erk$Phase), "G2",
                                                                   ifelse(grepl("M", phospho_data_long_erk$Phase), "M", NA)))))

# Handle missing values by imputing with column means
phospho_data_long_erk$Phospho_Values[is.na(phospho_data_long_erk$Phospho_Values)] <- 
  ave(phospho_data_long_erk$Phospho_Values, phospho_data_long_erk$Phase, FUN = function(x) mean(x, na.rm = TRUE))

# Check the final data for plotting
print(phospho_data_long_erk)

# Plot phosphorylation changes for ERK signaling
ggplot(phospho_data_long_erk, aes(y = Phase, x = Phospho_Values, color = PhaseCategory)) +
  geom_point(size = 3) +
  facet_wrap(~ Gene, scales = "free_x") +
  labs(title = "Phosphorylation Changes in ERK Signaling Pathway (All Phases)",
       y = "Cell Cycle Phase",
       x = "Phosphorylation Level") +
  scale_color_manual(values = c("earlyG1" = "blue", "lateG1" = "green", "S" = "red", "G2" = "purple", "M" = "orange")) +
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 45, hjust = 1),
        legend.position = "right")

# Function to print phosphorylation sites and their levels for each gene
print_gene_phospho_sites <- function(gene) {
  gene_data <- phospho_data_long_erk[phospho_data_long_erk$Gene == gene, ]
  gene_sites <- unique(gene_data$New_Gene_Residue)
  
  cat("Number of unique phosphorylation sites for", gene, ":", length(gene_sites), "\n")
  cat("Phosphorylation sites for", gene, ":\n")
  print(gene_sites)
  
  cat("Phosphorylation levels by phase for", gene, ":\n")
  print(gene_data)
  cat("\n\n")
}

# Print phosphorylation sites and their levels for each ERK signaling gene
for (gene in erk_signaling_genes) {
  print_gene_phospho_sites(gene)
}


#### phospho related sites among all phases for WNT

# Load necessary libraries
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# Load the data
data <- read.csv("C:/Users/avani/Downloads/Updated_Data_Mutated-June24.csv")

# Define WNT signaling genes
wnt_signaling_genes <- c("Ctnnb1", "Apc", "Gsk3b", "Axin1", "Lef1", "Tcf7l2")

# Extract phosphorylation data for WNT signaling genes
phospho_data_wnt <- data[data$Gene %in% wnt_signaling_genes, ]

# Check if data is extracted properly
print(phospho_data_wnt)

# If the data frame is empty, stop here and print a message
if (nrow(phospho_data_wnt) == 0) {
  stop("No data found for the specified WNT signaling genes.")
}

# Define all phases of interest
phases_of_interest <- c("phase_earlyG1_A", "phase_earlyG1_B", "phase_earlyG1_C", 
                        "phase_lateG1_A", "phase_lateG1_B", "phase_lateG1_C",
                        "phase_S_A", "phase_S_B", "phase_S_C",
                        "phase_G2_A", "phase_G2_B", "phase_G2_C",
                        "phase_M_A", "phase_M_B", "phase_M_C")

# Prepare the data for heatmap
phospho_data_wide_wnt <- phospho_data_wnt[, c("Gene", "New_Gene_Residue", phases_of_interest)]
phospho_matrix_wnt <- as.matrix(phospho_data_wide_wnt[,-(1:2)])
rownames(phospho_matrix_wnt) <- paste(phospho_data_wide_wnt$Gene, phospho_data_wide_wnt$New_Gene_Residue, sep = "_")

# Create a heatmap for WNT signaling genes
Heatmap(
  phospho_matrix_wnt, 
  name = "Phosphorylation",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "Cell Cycle Phase",
  row_title = "Gene and Residue",
  color = colorRamp2(c(min(phospho_matrix_wnt, na.rm = TRUE), 0, max(phospho_matrix_wnt, na.rm = TRUE)), c("blue", "white", "red"))
)

# Prepare the data for plotting
phospho_data_long_wnt <- melt(phospho_data_wnt, 
                              id.vars = c("Gene", "New_Gene_Residue"), 
                              measure.vars = phases_of_interest,
                              variable.name = "Phase", value.name = "Phospho_Values")

# Create a new column to categorize phases
phospho_data_long_wnt$PhaseCategory <- ifelse(grepl("earlyG1", phospho_data_long_wnt$Phase), "earlyG1",
                                              ifelse(grepl("lateG1", phospho_data_long_wnt$Phase), "lateG1",
                                                     ifelse(grepl("S", phospho_data_long_wnt$Phase), "S",
                                                            ifelse(grepl("G2", phospho_data_long_wnt$Phase), "G2",
                                                                   ifelse(grepl("M", phospho_data_long_wnt$Phase), "M", NA)))))

# Handle missing values by imputing with column means
phospho_data_long_wnt$Phospho_Values[is.na(phospho_data_long_wnt$Phospho_Values)] <- 
  ave(phospho_data_long_wnt$Phospho_Values, phospho_data_long_wnt$Phase, FUN = function(x) mean(x, na.rm = TRUE))

# Check the final data for plotting
print(phospho_data_long_wnt)

# Plot phosphorylation changes for WNT signaling
ggplot(phospho_data_long_wnt, aes(y = Phase, x = Phospho_Values, color = PhaseCategory)) +
  geom_point(size = 3) +
  facet_wrap(~ Gene, scales = "free_x") +
  labs(title = "Phosphorylation Changes in WNT Signaling Pathway (All Phases)",
       y = "Cell Cycle Phase",
       x = "Phosphorylation Level") +
  scale_color_manual(values = c("earlyG1" = "blue", "lateG1" = "green", "S" = "red", "G2" = "purple", "M" = "orange")) +
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 45, hjust = 1),
        legend.position = "right")

# Function to print phosphorylation sites and their levels for each gene
print_gene_phospho_sites <- function(gene) {
  gene_data <- phospho_data_long_wnt[phospho_data_long_wnt$Gene == gene, ]
  gene_sites <- unique(gene_data$New_Gene_Residue)
  
  cat("Number of unique phosphorylation sites for", gene, ":", length(gene_sites), "\n")
  cat("Phosphorylation sites for", gene, ":\n")
  print(gene_sites)
  
  cat("Phosphorylation levels by phase for", gene, ":\n")
  print(gene_data)
  cat("\n\n")
}

# Print phosphorylation sites and their levels for each WNT signaling gene
for (gene in wnt_signaling_genes) {
  print_gene_phospho_sites(gene)
}



###################################################
###############################################

## Fold change



# Load necessary libraries
library(reshape2) # For melt function
library(ggplot2) # For ggplot2
library(dplyr) # For data manipulation
library(tidyr) # For pivot_wider function

# Load the data
data <- read.csv("C:/Users/avani/Downloads/Updated_Data_Mutated-June24.csv")

# Define ERK signaling genes of interest
erk_signaling_genes <- c("Mapk1", "Mapk3", "Raf1", "Braf", "Map2k1", "Map2k2")

# Extract phosphorylation data for ERK signaling genes
phospho_data_erk <- data[data$Gene %in% erk_signaling_genes, ]

# Define the phases of interest (LG1 to S)
phases_of_interest <- c("phase_lateG1_A", "phase_lateG1_B", "phase_lateG1_C", 
                        "phase_S_A", "phase_S_B", "phase_S_C")

# Prepare the data for plotting
phospho_data_long_erk <- melt(phospho_data_erk, 
                              id.vars = c("Gene", "New_Gene_Residue"), 
                              measure.vars = phases_of_interest,
                              variable.name = "Phase", value.name = "Phosphorylation")

# Create a new column to categorize phases into LG1 and S
phospho_data_long_erk$PhaseCategory <- ifelse(grepl("lateG1", phospho_data_long_erk$Phase), "LG1", 
                                              ifelse(grepl("S", phospho_data_long_erk$Phase), "S", NA))

# Handle missing values by imputing with column means
phospho_data_long_erk$Phosphorylation[is.na(phospho_data_long_erk$Phosphorylation)] <- 
  ave(phospho_data_long_erk$Phosphorylation, phospho_data_long_erk$Phase, FUN = function(x) mean(x, na.rm = TRUE))

# Check the PhaseCategory column and filter out NA values
phospho_data_long_erk <- phospho_data_long_erk[!is.na(phospho_data_long_erk$PhaseCategory), ]

# Calculate mean phosphorylation levels for each Gene and PhaseCategory
mean_phospho <- phospho_data_long_erk %>%
  group_by(Gene, New_Gene_Residue, PhaseCategory) %>%
  summarise(MeanPhosphorylation = mean(Phosphorylation, na.rm = TRUE), .groups = 'drop')

# Pivot data to wide format for comparison
mean_phospho_wide <- pivot_wider(mean_phospho, 
                                 names_from = PhaseCategory, 
                                 values_from = MeanPhosphorylation,
                                 names_prefix = "MeanPhosphorylation_")

# Check the column names
colnames(mean_phospho_wide)

# Calculate log2 fold change between LG1 and S phases
if ("MeanPhosphorylation_LG1" %in% colnames(mean_phospho_wide) && 
    "MeanPhosphorylation_S" %in% colnames(mean_phospho_wide)) {
  
  mean_phospho_wide$Log2FoldChange <- log2((mean_phospho_wide$MeanPhosphorylation_LG1 + 1) / 
                                             (mean_phospho_wide$MeanPhosphorylation_S + 1))
  
} else {
  stop("Required columns are missing in the reshaped data.")
}

# Prepare data for plotting
plot_data <- mean_phospho_wide %>%
  select(Gene, New_Gene_Residue, Log2FoldChange) %>%
  melt(id.vars = c("Gene", "New_Gene_Residue"), variable.name = "Metric", value.name = "Value")

# Plot log2 fold change
ggplot(plot_data, aes(x = New_Gene_Residue, y = Value, color = Gene)) +
  geom_point(size = 3) +
  facet_wrap(~ Gene, scales = "free_y") +
  labs(title = "Log2 Fold Change in Phosphorylation (LG1 to S Phase)",
       x = "Phosphorylation Site",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")







# Load necessary libraries
library(reshape2) # For melt function
library(ggplot2) # For ggplot2
library(dplyr) # For data manipulation
library(tidyr) # For pivot_wider function

# Load the data
data <- read.csv("C:/Users/avani/Downloads/Updated_Data_Mutated-June24.csv")

# Define ERK signaling genes of interest
erk_signaling_genes <- c("Mapk1", "Mapk3", "Raf1", "Braf", "Map2k1", "Map2k2")

# Extract phosphorylation data for ERK signaling genes
phospho_data_erk <- data[data$Gene %in% erk_signaling_genes, ]

# Define the phases of interest
phases_of_interest <- c("phase_M_A", "phase_M_B", "phase_M_C", 
                        "phase_lateG1_A", "phase_lateG1_B", "phase_lateG1_C",
                        "phase_S_A", "phase_S_B", "phase_S_C")

# Prepare the data for plotting
phospho_data_long_erk <- melt(phospho_data_erk, 
                              id.vars = c("Gene", "New_Gene_Residue"), 
                              measure.vars = phases_of_interest,
                              variable.name = "Phase", value.name = "Phosphorylation")

# Create a new column to categorize phases
phospho_data_long_erk$PhaseCategory <- case_when(
  grepl("M", phospho_data_long_erk$Phase) ~ "M",
  grepl("lateG1", phospho_data_long_erk$Phase) ~ "LG1",
  grepl("S", phospho_data_long_erk$Phase) ~ "S",
  TRUE ~ NA_character_
)

# Handle missing values by imputing with column means
phospho_data_long_erk$Phosphorylation[is.na(phospho_data_long_erk$Phosphorylation)] <- 
  ave(phospho_data_long_erk$Phosphorylation, phospho_data_long_erk$Phase, FUN = function(x) mean(x, na.rm = TRUE))

# Check the PhaseCategory column and filter out NA values
phospho_data_long_erk <- phospho_data_long_erk[!is.na(phospho_data_long_erk$PhaseCategory), ]

# Calculate mean phosphorylation levels for each Gene and PhaseCategory
mean_phospho <- phospho_data_long_erk %>%
  group_by(Gene, New_Gene_Residue, PhaseCategory) %>%
  summarise(MeanPhosphorylation = mean(Phosphorylation, na.rm = TRUE), .groups = 'drop')

# Pivot data to wide format for comparison
mean_phospho_wide <- pivot_wider(mean_phospho, 
                                 names_from = PhaseCategory, 
                                 values_from = MeanPhosphorylation,
                                 names_prefix = "MeanPhosphorylation_")

# Check the column names to ensure correct naming
print("Column names in mean_phospho_wide:")
print(colnames(mean_phospho_wide))

# Calculate log10 fold change between each pair of phases
phases <- c("M", "LG1", "S")

# Initialize an empty list to store fold change data
log10_fold_changes <- list()

for (i in 1:(length(phases) - 1)) {
  for (j in (i + 1):length(phases)) {
    phase1 <- phases[i]
    phase2 <- phases[j]
    
    col1 <- paste0("MeanPhosphorylation_", phase1)
    col2 <- paste0("MeanPhosphorylation_", phase2)
    
    # Print column names being used for debugging
    print(paste("Processing:", col1, "and", col2))
    
    if (col1 %in% colnames(mean_phospho_wide) && col2 %in% colnames(mean_phospho_wide)) {
      log10_fold_changes[[paste(phase1, "vs", phase2)]] <- mean_phospho_wide %>%
        mutate(Log10FoldChange = log10((get(col1) + 1) / (get(col2) + 1))) %>%
        select(Gene, New_Gene_Residue, Log10FoldChange) %>%
        mutate(PhaseComparison = paste(phase1, "vs", phase2))
    } else {
      warning(paste("Columns", col1, "or", col2, "do not exist."))
    }
  }
}

# Combine all fold changes into one data frame
log10_fold_changes_df <- bind_rows(log10_fold_changes)

# Plot log10 fold changes
ggplot(log10_fold_changes_df, aes(x = New_Gene_Residue, y = Log10FoldChange, color = Gene)) +
  geom_point(size = 3) +
  facet_wrap(~ PhaseComparison, scales = "free_y") +
  labs(title = "Log10 Fold Change in Phosphorylation Across Phases",
       x = "Phosphorylation Site",
       y = "Log10 Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")


## Each dot shows the log10 fold change in phosphorylation for a particular
##phosphorylation site. A positive value indicates increased phosphorylation
##in the first phase compared to the second phase, while a negative value indicates 
##decreased phosphorylation.



