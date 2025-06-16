############################################################################################
# This script processes an RNA-seq counts matrix and applies an unsupervised machine learning
# technique—Principal Component Analysis (PCA)—for dimensionality reduction. The goal is to 
# visualize the natural clustering of *E. coli* samples treated with 17 different antibiotics 
# and reduce the complexity of the gene expression data.
#
# This highlights how machine learning can uncover pharmacological patterns in organisms,
# aiding drug discovery efforts.
#
# The analysis continues in Python with supervised learning (e.g., Random Forest)
# to predict the mechanism of action (MoA), assessing the model's ability to distinguish
# and classify antibiotics based on their gene expression profiles.
############################################################################################

setwd("C:/Users/cndobela/Documents/bioinfo_Masters_Project/Project/understanding MOA/ML/analysis/salmon workflow")

library(tximport)
library(GenomicFeatures)
library(readr)
library(txdbmaker)
library(tidyverse)
library(ggplot2)



#Generating counts from Salmon quant output
samples <- read.table("filereport_read_run_PRJNA432937_tsv.txt", header = TRUE)
files <- file.path("quant", samples$run_accession, "quant.sf")
names(files) <- paste0(samples$run_accession)
txi.salmon <- tximport(files, type = "salmon", txIn = FALSE, geneIdCol = "Name")
head(txi.salmon$counts)
counts <- txi.salmon$counts

#checking for missing values
# 1. Check for NA values
sum(is.na(counts))
# 2. Check for "-" as missing values 
sum(counts == "-", na.rm = TRUE)

#save counts as csv
write.csv(counts, "Counts_salmon.csv",row.names = T)

#Shape of data
ncol(counts)
nrow(counts)

#transposing data
transposed_data <- t(counts)
write.csv(transposed_data, "Counts_salmon_T.csv",row.names = T)

Counts_salmon_T <- read.csv("Counts_salmon_T.csv", header = T)
counts_T <- t(Counts_salmon_T, header = T)
tail(counts_T)


# Set the first row as header
colnames(counts_T) <- as.character(unlist(counts_T[1, ]))

# Remove the first row
counts_T<- counts_T[-1, ]

write.csv(counts_T, "counts_salmon2.csv",row.names = T)

#######################normalise with TPM ################
CountsR <- read.csv("CountsR.csv", header = T)

# Sum counts for duplicate genes and keep first length
CountsR <- CountsR %>%
  group_by(Symbol) %>%
  summarise(
    across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),  # Sum counts
    length = first(length)  # Keep the first gene length
  ) %>%
  ungroup()

# Separate counts and gene lengths
gene_lengths <- CountsR$length
count_matrix <- CountsR[, !(colnames(CountsR) %in% c('length', 'Symbol'))]  # remove gene_len from counts

# check for duplicate genes
dup <- CountsR %>%
  group_by(Symbol) %>%
  filter(n() > 1) %>%
  arrange(Symbol)


# Before summing (list duplicates)
 duplicates <- CountsR$Symbol[duplicated(CountsR$Symbol)]
 print(unique(duplicates))


# calculate tpm
 
calculate_tpm_matrix <- function(count_matrix, gene_lengths) {
  # Convert gene lengths to kilobases
  gene_lengths_kb <- gene_lengths / 1000
  
  # Calculate RPK for each sample
  rpk_matrix <- sweep(count_matrix, 1, gene_lengths_kb, "/")
  
  # Calculate scaling factor for each sample
  scaling_factors <- colSums(rpk_matrix) / 1e6
  
  # Calculate TPM for each sample
  tpm_matrix <- sweep(rpk_matrix, 2, scaling_factors, "/")
  
  return(tpm_matrix)
}
tpm_matrix <- calculate_tpm_matrix(count_matrix, gene_lengths)
print(tpm_matrix)

#sum of each column to varyfy tpm, they should all be eqaul!
colSums(tpm_matrix)

#save tpm_matrix
write.csv(tpm_matrix, "tpm_matrix.csv",row.names = T)

# Extract the "Symbol" column from CountsR
#symbols <- CountsR$Symbol  

# Add it as the first column in tpm_matrix
#tpm_matrix <- cbind(Symbol = symbols, tpm_matrix)



##############z-scoring standardization#########
# Log2 transform
log_tpm <- log2(tpm_matrix + 1)

# check for zeros the program rounded off 
rows_all_zeros <- apply(log_tpm, 1, function(row) all(row == 0))
which(rows_all_zeros)  # Indices of fully zero rows

# Extract the "Symbol" column from CountsR
symbols <- CountsR$Symbol  

# Add it as the first column in tpm_matrix
log_tpm <- cbind(Symbol = symbols, log_tpm)

#remove rows with zero
log_tpm_filtered <- log_tpm[!rows_all_zeros, ]

# Check if any fully zero rows remain
any(apply(log_tpm_filtered, 1, function(row) all(row == 0)))  # Should return FALSE

# Dimensions before vs. after
dim(log_tpm)          # Original
dim(log_tpm_filtered)  # Filtered

# remove column symbol from log_tpm_filtered
Gene <- log_tpm_filtered$Symbol
log_tpm <- log_tpm_filtered[, !colnames(log_tpm_filtered) %in% "Symbol"]

# Z-score per gene
z_scores <- t(apply(log_tpm, 1, function(x) (x - mean(x)) / sd(x)))

# Save result
write.csv(z_scores, "zscore_matrix.csv")
###################PCA#######################
pca_data <- prcomp(t(z_scores), scale. = FALSE)

# View summary of PCA
summary(pca_data)

pca_scores = as.data.frame(pca_data$x)

# Get variance explained by each PC
var_explained <- pca_data$sdev^2 / sum(pca_data$sdev^2)


# Create a data frame from var_explained
data <- data.frame(Principal_Component = 1:length(var_explained),
                   Variance_Explained = var_explained)

#Extract summary statistics
eigenvalues <- pca_data$sdev^2
proportion <- eigenvalues / sum(eigenvalues)
cumulative <- cumsum(proportion)

# Create summary table
summary_table <- data.frame(
  PC = 1:length(eigenvalues),
  Eigenvalue = eigenvalues,
  Proportion = proportion,
  Cumulative = cumulative
)

##########plots############################
# Scree plot
Scree_plot <- ggplot(data, aes(x = Principal_Component, y = Variance_Explained)) +
  geom_point(color = "steelblue", size = 0.8) +  # Add points
  geom_line(color = "steelblue") +             # Add line
  labs(x = "Principal Component", y = "Eigenvalues") + 
  theme_minimal()

ggsave("Scree_plot.png", 
       plot = Scree_plot,
       width = 4,
       height = 2,
       dpi = 300)
###########################renaming colunms of sample to MoA###########################
samples <- c("PDD_P2_01", "PDD_P2_02", "PDD_P2_03", "PDD_P2_04", "PDD_P2_05", "PDD_P2_06", 
             "PDD_P2_07", "PDD_P2_09", "PDD_P2_10", "PDD_P2_11", "PDD_P2_12", "PDD_P2_14", 
             "PDD_P2_15", "PDD_P2_18", "PDD_P2_20", "PDD_P2_21", "PDD_P2_23", "PDD_P2_24", 
             "PDD_P2_26", "PDD_P2_27", "PDD_P2_28", "PDD_P2_29", "PDD_P2_30", "PDD_P2_31", 
             "PDD_P2_32", "PDD_P2_33", "PDD_P2_35", "PDD_P2_36", "PDD_P2_37", "PDD_P2_38", 
             "PDD_P2_40", "PDD_P2_41", "PDD_P2_44", "PDD_P2_46", "PDD_P2_47", "PDD_P2_49", 
             "PDD_P2_50", "PDD_P2_52", "PDD_P2_53", "PDD_P2_54", "PDD_P2_55", "PDD_P2_56", 
             "PDD_P2_57", "PDD_P2_58", "PDD_P2_59", "PDD_P2_61", "PDD_P2_62", "PDD_P2_63", 
             "PDD_P2_64", "PDD_P2_66", "PDD_P2_67", "PDD_P2_70", "PDD_P2_72", "PDD_P2_73", 
             "PDD_P2_75", "PDD_P2_76", "PDD_P2_78")
moa <- c(
  "control",
  "control",
  "DNA damage",
  "cell wall synthesis inhibitor",
  "Protein synthesis",
  "DNA replication inhibitor",
  "Protein synthesis",
  "Protein synthesis",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "new hit",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "DNA replication inhibitor",
  "Folic Acid synthesis inhibitor",
  "control",
  "control",
  "DNA damage",
  "cell wall synthesis inhibitor",
  "Protein synthesis",
  "DNA replication inhibitor",
  "Protein synthesis",
  "Protein synthesis",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "new hit",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "DNA replication inhibitor",
  "Folic Acid synthesis inhibitor",
  "control",
  "control",
  "DNA damage",
  "cell wall synthesis inhibitor",
  "Protein synthesis",
  "DNA replication inhibitor",
  "Protein synthesis",
  "Protein synthesis",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "new hit",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "DNA replication inhibitor",
  "Folic Acid synthesis inhibitor"
)


# Create a named vector mapping samples to MOA
names(moa) <- samples

# Assuming 'tpm_data' is your data frame, and samples are in columns 3 to 59, we will rename those columns based on the MOA
colnames(z_scores)[1:57] <- moa[match(colnames(z_scores)[1:57], samples)]

# View the updated column names
head(colnames(z_scores)[1:57])

#check how many colomuns i have in each position
data.frame(Position = 1:ncol(z_scores), Column_Name = colnames(z_scores))

z_scores = as.data.frame(z_scores)



# Define color palette for each MOA (ensure the color vector is the same length as moa)
colors <- c(
  "control" = "blue",
  "DNA damage" = "green",
  "cell wall synthesis inhibitor" = "purple",
  "Protein synthesis" = "orange",
  "DNA replication inhibitor" = "red",
  "Membrane perturbation" = "yellow",
  "new hit" = "pink",
  "Folic Acid synthesis inhibitor" = "brown"
)

# Create color vector corresponding to each entry in moa
moa_colors <- colors[moa]

# Plotting (assuming pca_data$x is already available)
library(ggplot2)
#######################pc plot for MOA ##################
# Convert PCA data to a data frame for ggplot
pca_df <- data.frame(
  PC1 = pca_data$x[,1],
  PC2 = pca_data$x[,2],
  MOA = moa  # Add MOA as a column for coloring
)

# Create ggplot
pc_2D_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = MOA)) +
  geom_point(size = 3) +  # Adjust point size
  scale_color_manual(
    name = "Mechanism of Action",
    values = colors
  ) +  # <-- Fixed: Removed extra closing parenthesis
  labs(
    x = "PC1 (29.58%)",
    y = "PC2 (22.73%)"
  ) +
  theme_minimal() +
  theme(legend.position = "right")
# Position legend on the right

# Save the plot
ggsave("pc_2D_plot.png", 
       plot = pc_2D_plot,
       width = 8,
       height = 5,
       dpi = 300)

# Create 3D scatter plot with color-coded markers
plot_ly(
  x = pca_data$x[,1],
  y = pca_data$x[,2],
  z = pca_data$x[,3],
  type = "scatter3d",
  mode = "markers",
  marker = list(
    size = 5,
    opacity = 0.8,
    color = moa_colors
  ),
  text = moa,  # optional: hover text to show MOA
  hoverinfo = "text"
) %>%
  layout(
    scene = list(
      xaxis = list(title = "PC1 (29.58%)"),
      yaxis = list(title = "PC2 (22.73%)"),
      zaxis = list(title = "PC3 (11.36%)")
    )
  )

##########################cummulative_plot#############

cum_plot <- ggplot(summary_table, aes(x = PC, y = Cumulative)) +
  # Plot elements (keep your original sizes)
  geom_line(color = "steelblue", linewidth = 0.5) +
  geom_point(color = "steelblue", size = 1.5) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "red", linewidth = 0.6) +
  
  # Y-axis (unchanged)
  scale_y_continuous(
    name = "Cumulative Variance",
    limits = c(0, 1),
    sec.axis = sec_axis(~ . * 100, name = "Cumulative Variance (%)", 
                        labels = scales::percent_format(scale = 1))
  ) +
  
  # X-axis modified to show PCs in steps of 2
  scale_x_continuous(
    name = "Principal Component (PC)",
    breaks = seq(2, max(summary_table$PC), by = 2),  # Key change: intervals of 2
    expand = expansion(mult = 0.05)
  ) +
  
  # Your original theme settings
  theme_minimal() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 8, face = "bold"),
    axis.title.y.right = element_text(color = "black", margin = margin(l = 10)),
    panel.grid.minor = element_blank()
  )

# Save (unchanged)
ggsave("cumulative_plot.png", 
       plot = cum_plot,
       width = 8,
       height = 5,
       dpi = 300)


##################Z-score_curve#############
library(ggplot2)

# Freedman-Diaconis bin calculation
bin_width <- 2 * IQR(z_scores_vector) / (length(z_scores_vector))^(1/3)
n_bins_fd <- round((max(z_scores_vector) - min(z_scores_vector)) / bin_width)

# Print for debugging
cat("Freedman-Diaconis bins:", n_bins_fd, "\n")

# Plot
Z_score_curve <- ggplot(data.frame(z = z_scores_vector), aes(x = z)) +
  geom_histogram(
    aes(y = ..density..), 
    bins = n_bins_fd,  # according to theoery on no.of samples
    alpha = 0.5, 
    fill = "gray",
    color = "steelblue"    # Border for clarity
  ) +
  stat_function(
    fun = dnorm,
    args = list(mean = mean(z_scores_vector), sd = sd(z_scores_vector)),
    color = "blue", 
    linewidth = 1
  ) +
  labs(
    x = "Z-Score",
    y = "Density"
  ) +
  theme_minimal()

ggsave("Z_score_curve.png", 
       plot = Z_score_curve,
       width = 8,
       height = 5,
       dpi = 300)
###############################for compounds#######################################
samples <- c("PDD_P2_01", "PDD_P2_02", "PDD_P2_03", "PDD_P2_04", "PDD_P2_05", "PDD_P2_06", 
             "PDD_P2_07", "PDD_P2_09", "PDD_P2_10", "PDD_P2_11", "PDD_P2_12", "PDD_P2_14", 
             "PDD_P2_15", "PDD_P2_18", "PDD_P2_20", "PDD_P2_21", "PDD_P2_23", "PDD_P2_24", 
             "PDD_P2_26", "PDD_P2_27", "PDD_P2_28", "PDD_P2_29", "PDD_P2_30", "PDD_P2_31", 
             "PDD_P2_32", "PDD_P2_33", "PDD_P2_35", "PDD_P2_36", "PDD_P2_37", "PDD_P2_38", 
             "PDD_P2_40", "PDD_P2_41", "PDD_P2_44", "PDD_P2_46", "PDD_P2_47", "PDD_P2_49", 
             "PDD_P2_50", "PDD_P2_52", "PDD_P2_53", "PDD_P2_54", "PDD_P2_55", "PDD_P2_56", 
             "PDD_P2_57", "PDD_P2_58", "PDD_P2_59", "PDD_P2_61", "PDD_P2_62", "PDD_P2_63", 
             "PDD_P2_64", "PDD_P2_66", "PDD_P2_67", "PDD_P2_70", "PDD_P2_72", "PDD_P2_73", 
             "PDD_P2_75", "PDD_P2_76", "PDD_P2_78")

compounds <- c(
  "DMSO", "DMSO", "Nitrofurantoin", "Meropenem", "Chloramphenicol", "Nitroxolin",
  "Clarithromycin", "Doxycycline", "Levofloxacin", "Ceftriaxone", "PolymyxineB",
  "Mecillinam", "Colistin", "AZ-LolCDE", "B-01", "Ciprofloxacin", "Globomycin",
  "Norfloxacin", "Trimethoprim", "DMSO", "DMSO", "Nitrofurantoin", "Meropenem",
  "Chloramphenicol", "Nitroxolin", "Clarithromycin", "Doxycycline", "Levofloxacin",
  "Ceftriaxone", "PolymyxineB", "Mecillinam", "Colistin", "AZ-LolCDE", "B-01",
  "Ciprofloxacin", "Globomycin", "Norfloxacin", "Trimethoprim", "DMSO", "DMSO",
  "Nitrofurantoin", "Meropenem", "Chloramphenicol", "Nitroxolin", "Clarithromycin",
  "Doxycycline", "Levofloxacin", "Ceftriaxone", "PolymyxineB", "Mecillinam",
  "Colistin", "AZ-LolCDE", "B-01", "Ciprofloxacin", "Globomycin", "Norfloxacin",
  "Trimethoprim"
)

moa <- c(
  "control",
  "control",
  "DNA damage",
  "cell wall synthesis inhibitor",
  "Protein synthesis",
  "DNA replication inhibitor",
  "Protein synthesis",
  "Protein synthesis",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "new hit",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "DNA replication inhibitor",
  "Folic Acid synthesis inhibitor",
  "control",
  "control",
  "DNA damage",
  "cell wall synthesis inhibitor",
  "Protein synthesis",
  "DNA replication inhibitor",
  "Protein synthesis",
  "Protein synthesis",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "new hit",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "DNA replication inhibitor",
  "Folic Acid synthesis inhibitor",
  "control",
  "control",
  "DNA damage",
  "cell wall synthesis inhibitor",
  "Protein synthesis",
  "DNA replication inhibitor",
  "Protein synthesis",
  "Protein synthesis",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "new hit",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "DNA replication inhibitor",
  "Folic Acid synthesis inhibitor"
)

# Create a named vector mapping samples to MOA
names(compounds) <- samples

# Assuming 'tpm_data' is your data frame, and samples are in columns 3 to 59, we will rename those columns based on the MOA
colnames(z_scores)[1:57] <- compounds[match(colnames(z_scores)[1:57], samples)]

# View the updated column names
head(colnames(z_scores)[1:57])

#check how many colomuns i have in each position
data.frame(Position = 1:ncol(z_scores), Column_Name = colnames(z_scores))

z_scores = as.data.frame(z_scores)



# Define compound color
compounds <- c(
  "DMSO", "DMSO", "Nitrofurantoin", "Meropenem", "Chloramphenicol", "Nitroxolin",
  "Clarithromycin", "Doxycycline", "Levofloxacin", "Ceftriaxone", "PolymyxineB",
  "Mecillinam", "Colistin", "AZ-LolCDE", "B-01", "Ciprofloxacin", "Globomycin",
  "Norfloxacin", "Trimethoprim", "DMSO", "DMSO", "Nitrofurantoin", "Meropenem",
  "Chloramphenicol", "Nitroxolin", "Clarithromycin", "Doxycycline", "Levofloxacin",
  "Ceftriaxone", "PolymyxineB", "Mecillinam", "Colistin", "AZ-LolCDE", "B-01",
  "Ciprofloxacin", "Globomycin", "Norfloxacin", "Trimethoprim", "DMSO", "DMSO",
  "Nitrofurantoin", "Meropenem", "Chloramphenicol", "Nitroxolin", "Clarithromycin",
  "Doxycycline", "Levofloxacin", "Ceftriaxone", "PolymyxineB", "Mecillinam",
  "Colistin", "AZ-LolCDE", "B-01", "Ciprofloxacin", "Globomycin", "Norfloxacin",
  "Trimethoprim"
)

# Define compound color palette
colors <- c(
  "DMSO" = "blue",          # Light gray for control
  "Nitrofurantoin" = "#1F78B4",  # Dark blue
  "Meropenem" = "#33A02C",    # Dark green
  "Chloramphenicol" = "#E31A1C", # Red
  "Nitroxolin" = "#FF7F00",    # Orange
  "Clarithromycin" = "#6A3D9A", # Purple
  "Doxycycline" = "#B15928",   # Brown
  "Levofloxacin" = "#A6CEE3",  # Light blue
  "Ceftriaxone" = "#B2DF8A",   # Light green
  "PolymyxineB" = "#FB9A99",   # Pink
  "Mecillinam" = "#FDBF6F",    # Light orange
  "Colistin" = "#CAB2D6",      # Light purple
  "AZ-LolCDE" = "#FFFF99",     # Light yellow
  "B-01" = "purple",          # Medium gray
  "Ciprofloxacin" = "#1B9E77", # Teal
  "Globomycin" = "#D95F02",    # Dark orange
  "Norfloxacin" = "#7570B3",   # Medium purple
  "Trimethoprim" = "#E7298A",  # Magenta
  "Rifampicin" = "#66A61E",    # Olive green
  "Vancomycin" = "#E6AB02"     # Gold
)

# Create color vector corresponding to each entry in moa
moa_colors <- colors[compounds]

# Plotting (assuming pca_data$x is already available)
library(ggplot2)
library(plotly)
#######################pc plot##################
# Convert PCA data to a data frame for ggplot
pca_df <- data.frame(
  PC1 = pca_data$x[,1],
  PC2 = pca_data$x[,2],
  Compounds = compounds  # Add MOA as a column for coloring
)

# Create ggplot
pc_2D_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = compounds)) +
  geom_point(size = 3) +  # Adjust point size
  scale_color_manual(
    name = "Compounds",
    values = colors
  ) +
  labs(
    x = "PC1 (29.58%)",
    y = "PC2 (22.73%)"
  ) +
  theme_minimal() +
  theme(legend.position = "right")  # Position legend on the right

# Save the plot
ggsave("pc_2D_plot_1.png", 
       plot = pc_2D_plot,
       width = 8,
       height = 5,
       dpi = 300)

# Create 3D scatter plot with color-coded markers
plot_ly(
  x = pca_data$x[,1],
  y = pca_data$x[,2],
  z = pca_data$x[,3],
  type = "scatter3d",
  mode = "markers",
  marker = list(
    size = 5,
    opacity = 0.8,
    color = moa_colors
  ),
  text = compounds,  # optional: hover text to show MOA
  hoverinfo = "text"
) %>%
  layout(
    scene = list(
      xaxis = list(title = "PC1 (29.58%)"),
      yaxis = list(title = "PC2 (22.73%)"),
      zaxis = list(title = "PC3 (11.36%)")
    )
  )


############################selecting pc to export#############
samples <- c("PDD_P2_01", "PDD_P2_02", "PDD_P2_03", "PDD_P2_04", "PDD_P2_05", "PDD_P2_06", 
             "PDD_P2_07", "PDD_P2_09", "PDD_P2_10", "PDD_P2_11", "PDD_P2_12", "PDD_P2_14", 
             "PDD_P2_15", "PDD_P2_18", "PDD_P2_20", "PDD_P2_21", "PDD_P2_23", "PDD_P2_24", 
             "PDD_P2_26", "PDD_P2_27", "PDD_P2_28", "PDD_P2_29", "PDD_P2_30", "PDD_P2_31", 
             "PDD_P2_32", "PDD_P2_33", "PDD_P2_35", "PDD_P2_36", "PDD_P2_37", "PDD_P2_38", 
             "PDD_P2_40", "PDD_P2_41", "PDD_P2_44", "PDD_P2_46", "PDD_P2_47", "PDD_P2_49", 
             "PDD_P2_50", "PDD_P2_52", "PDD_P2_53", "PDD_P2_54", "PDD_P2_55", "PDD_P2_56", 
             "PDD_P2_57", "PDD_P2_58", "PDD_P2_59", "PDD_P2_61", "PDD_P2_62", "PDD_P2_63", 
             "PDD_P2_64", "PDD_P2_66", "PDD_P2_67", "PDD_P2_70", "PDD_P2_72", "PDD_P2_73", 
             "PDD_P2_75", "PDD_P2_76", "PDD_P2_78")


moa <- c(
  "control",
  "control",
  "DNA damage",
  "cell wall synthesis inhibitor",
  "Protein synthesis",
  "DNA replication inhibitor",
  "Protein synthesis",
  "Protein synthesis",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "new hit",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "DNA replication inhibitor",
  "Folic Acid synthesis inhibitor",
  "control",
  "control",
  "DNA damage",
  "cell wall synthesis inhibitor",
  "Protein synthesis",
  "DNA replication inhibitor",
  "Protein synthesis",
  "Protein synthesis",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "new hit",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "DNA replication inhibitor",
  "Folic Acid synthesis inhibitor",
  "control",
  "control",
  "DNA damage",
  "cell wall synthesis inhibitor",
  "Protein synthesis",
  "DNA replication inhibitor",
  "Protein synthesis",
  "Protein synthesis",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "new hit",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "DNA replication inhibitor",
  "Folic Acid synthesis inhibitor"
)

# Create a named vector mapping samples to MOA
names(moa) <- samples

t_pca_scores<- t(pca_scores)

colnames(t_pca_scores)
dim(t_pca_scores)
ncol(t_pca_scores)
nrow(t_pca_scores)

# Assuming 'tpm_data' is your data frame, and samples are in columns 3 to 59, we will rename those columns based on the MOA
colnames(t_pca_scores)[1:57] <- moa[match(colnames(t_pca_scores)[1:57], samples)]

# View the updated column names
head(colnames(t_pca_scores)[1:57])

#check how many colomuns i have in each position
data.frame(Position = 1:ncol(t_pca_scores), Column_Name = colnames(t_pca_scores))

t_pca_scores = as.data.frame(t_pca_scores)
pca_scores <- t(t_pca_scores)%>%
  as.data.frame()

write.csv(pca_scores, "pca_scores_ml.csv")

#########################log_filtered data comparison###########

samples <- c("PDD_P2_01", "PDD_P2_02", "PDD_P2_03", "PDD_P2_04", "PDD_P2_05", "PDD_P2_06", 
             "PDD_P2_07", "PDD_P2_09", "PDD_P2_10", "PDD_P2_11", "PDD_P2_12", "PDD_P2_14", 
             "PDD_P2_15", "PDD_P2_18", "PDD_P2_20", "PDD_P2_21", "PDD_P2_23", "PDD_P2_24", 
             "PDD_P2_26", "PDD_P2_27", "PDD_P2_28", "PDD_P2_29", "PDD_P2_30", "PDD_P2_31", 
             "PDD_P2_32", "PDD_P2_33", "PDD_P2_35", "PDD_P2_36", "PDD_P2_37", "PDD_P2_38", 
             "PDD_P2_40", "PDD_P2_41", "PDD_P2_44", "PDD_P2_46", "PDD_P2_47", "PDD_P2_49", 
             "PDD_P2_50", "PDD_P2_52", "PDD_P2_53", "PDD_P2_54", "PDD_P2_55", "PDD_P2_56", 
             "PDD_P2_57", "PDD_P2_58", "PDD_P2_59", "PDD_P2_61", "PDD_P2_62", "PDD_P2_63", 
             "PDD_P2_64", "PDD_P2_66", "PDD_P2_67", "PDD_P2_70", "PDD_P2_72", "PDD_P2_73", 
             "PDD_P2_75", "PDD_P2_76", "PDD_P2_78")


moa <- c(
  "control",
  "control",
  "DNA damage",
  "cell wall synthesis inhibitor",
  "Protein synthesis",
  "DNA replication inhibitor",
  "Protein synthesis",
  "Protein synthesis",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "new hit",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "DNA replication inhibitor",
  "Folic Acid synthesis inhibitor",
  "control",
  "control",
  "DNA damage",
  "cell wall synthesis inhibitor",
  "Protein synthesis",
  "DNA replication inhibitor",
  "Protein synthesis",
  "Protein synthesis",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "new hit",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "DNA replication inhibitor",
  "Folic Acid synthesis inhibitor",
  "control",
  "control",
  "DNA damage",
  "cell wall synthesis inhibitor",
  "Protein synthesis",
  "DNA replication inhibitor",
  "Protein synthesis",
  "Protein synthesis",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "Membrane perturbation",
  "cell wall synthesis inhibitor",
  "new hit",
  "DNA replication inhibitor",
  "cell wall synthesis inhibitor",
  "DNA replication inhibitor",
  "Folic Acid synthesis inhibitor"
)

#remove rows with zero
log_tpm_filtered <- log_tpm[!rows_all_zeros, ]

# Create a named vector mapping samples to MOA
names(moa) <- samples

dim(log_tpm_filtered)


# Assuming 'tpm_data' is your data frame, and samples are in columns 3 to 59, we will rename those columns based on the MOA
colnames(log_tpm_filtered)[1:58] <- moa[match(colnames(log_tpm_filtered)[1:58], samples)]

# View the updated column names
head(colnames(log_tpm_filtered)[1:58])

#check how many colomuns i have in each position
data.frame(Position = 1:ncol(log_tpm_filtered), Column_Name = colnames(log_tpm_filtered))

log_tpm_filtered = as.data.frame(log_tpm_filtered)
log_tpm_filtered <- t(log_tpm_filtered)%>%
  as.data.frame()


write.csv(log_tpm_filtered, "log_tpm.csv")
