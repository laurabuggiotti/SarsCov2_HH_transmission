# ==============================================================================
# SARS-CoV-2 Transmission Direction Inference - Supplementary Material
# Supplementary Figure 6: TransPhylo-based transmission probability matrices
# 
# This script performs transmission direction inference using TransPhylo and 
# creates heatmap visualizations showing transmission probabilities between
# individuals within households, organized by viral lineage.
#
# The analysis includes:
# 1. Phylogenetic tree dating using BactDating
# 2. Transmission tree inference using TransPhylo
# 3. Who-infected-whom (WIW) probability matrix calculation
# 4. Lineage-specific transmission probability heatmaps
#
# Input: 
#   - vw_noProbSiteMutPos_full.12.raxml.bestTree - Phylogenetic tree
#   - metadata_tree.txt - Sample collection dates
#   - lineage_VL.csv - Sample metadata with household and lineage information
# Output: 
#   - Transmission probability heatmaps by viral lineage
#   - TransPhylo analysis plots and traces
# ==============================================================================

# Load required libraries
library(ape)          # for phylogenetic tree analysis
library(BactDating)   # for bacterial/viral dating
library(TransPhylo)   # for transmission tree inference
library(lubridate)    # for date handling
library(dplyr)        # for data manipulation
library(ggplot2)      # for plotting
library(reshape2)     # for data reshaping
library(lattice)      # for levelplot
library(coda)         # for MCMC diagnostics

# Set file paths (adjust as needed for your directory structure)
tree_path <- "data/vw_noProbSiteMutPos_full.12.raxml.bestTree"
metadata_path <- "data/metadata_tree.txt"
lineage_path <- "data/lineage_VL.csv"
output_dir <- "figures"
analysis_dir <- "analysis_output"

# Create output directories if they don't exist
for (dir in c(output_dir, analysis_dir)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

# ==============================================================================
# PART 1: Load and prepare phylogenetic tree
# ==============================================================================

cat("Loading phylogenetic tree from:", tree_path, "\n")
t_full <- read.tree(tree_path)

# Process tree for TransPhylo analysis
t <- t_full

cat("Original tree branch length sum:", sum(t$edge.length), "\n")

# Scale branch lengths to genome length (substitutions per genome)
# SARS-CoV-2 genome length: ~29,903 nucleotides
t$edge.length <- t$edge.length * 29903

cat("Scaled tree branch length sum:", sum(t$edge.length), "\n")

# Clean tip labels (remove suffix after dash)
t$tip.label <- sapply(strsplit(t$tip.label, "\\-"), '[', 1)

cat("Tree contains", Ntip(t), "tips after processing\n")

# ==============================================================================
# PART 2: Load and process collection date metadata
# ==============================================================================

cat("Loading collection date metadata from:", metadata_path, "\n")
metad <- read.table(metadata_path, header = TRUE)

cat("Metadata contains", nrow(metad), "records\n")

# Convert collection dates to decimal years
metad$CollectionDate <- decimal_date(as.Date(metad$CollectionDate, '%d/%m/%Y'))

# Create named vector linking tip labels to dates
dates <- setNames(metad$CollectionDate, metad$tipLabel)

# Reorder to match tree tip labels
dates <- dates[t$tip.label]

cat("Date range:", min(dates, na.rm = TRUE), "to", max(dates, na.rm = TRUE), "\n")
cat("Number of dated samples:", sum(!is.na(dates)), "\n")

# ==============================================================================
# PART 3: Phylogenetic dating with BactDating
# ==============================================================================

cat("Performing phylogenetic dating analysis...\n")

# Find optimal root for the tree
cat("Finding optimal root position...\n")
rooted_full <- initRoot(t, dates)

# Assess temporal signal with root-to-tip regression
cat("Assessing root-to-tip temporal signal...\n")
pdf(file.path(analysis_dir, "root_to_tip_regression.pdf"))
r <- roottotip(rooted_full, dates)
dev.off()

cat("Root-to-tip correlation: R² =", round(r$R2, 4), "\n")

# Perform Bayesian dating analysis
cat("Running BactDating MCMC analysis (this may take several minutes)...\n")
res_rooted_full <- bactdate(rooted_full, dates, nbIts = 1e5)

# Save results
saveRDS(res_rooted_full, file = file.path(analysis_dir, "bactdating_results.rds"))

# Generate diagnostic plots
cat("Creating BactDating diagnostic plots...\n")

# MCMC trace plots
pdf(file.path(analysis_dir, "bactdating_trace.pdf"))
plot(res_rooted_full, 'trace')
dev.off()

# Tree with confidence intervals
pdf(file.path(analysis_dir, "dated_tree_CI.pdf"))
plot(res_rooted_full, 'treeCI', cex = 0.2)
dev.off()

# Inferred root
pdf(file.path(analysis_dir, "inferred_root.pdf"))
plot(res_rooted_full, 'treeRoot', cex = 0.3)
dev.off()

# Extract dated tree
dated_tree <- res_rooted_full$tree

# ==============================================================================
# PART 4: Load sample metadata and prepare for TransPhylo
# ==============================================================================

cat("Loading sample metadata from:", lineage_path, "\n")
vw_metadata <- read.csv(lineage_path)

cat("Sample metadata contains", nrow(vw_metadata), "records\n")

# Filter for household transmission pairs (hh2)
metad_hh2 <- vw_metadata %>% filter(hh_type == 'hh2')
cat("Household transmission pairs (hh2):", nrow(metad_hh2), "samples\n")

# Update tree tip labels with unique household identifiers
dated_tree$tip.label[match(vw_metadata$SAMPLE, dated_tree$tip.label)] <- vw_metadata$hh_type_unique_b

# ==============================================================================
# PART 5: TransPhylo transmission tree inference
# ==============================================================================

cat("Preparing data for TransPhylo analysis...\n")

# Convert to ptree format for TransPhylo
obs_end <- lubridate::decimal_date(as.Date('2023/7/1'))
pt <- ptreeFromPhylo(dated_tree, obs_end)

cat("Created ptree object for TransPhylo\n")

# Quick initial inference
cat("Running initial TransPhylo inference...\n")
res_initial <- inferTTree(pt, w.shape = 0.01, w.scale = 1)

# Main MCMC inference with specific parameters
cat("Running full TransPhylo MCMC inference (this will take considerable time)...\n")
cat("Using parameters optimized for within-household transmission\n")

# Set pi to high value as we're interested in within-household directionality
mcmc_Tree2 <- inferTTree(pt, 
                        w.mean = 0.01,
                        w.std = 0.009,
                        startPi = 0.99, 
                        updatePi = FALSE,
                        mcmcIterations = 1000000)

# Save TransPhylo results
saveRDS(mcmc_Tree2, file = file.path(analysis_dir, "transphylo_results.rds"))

# Plot transmission tree
pdf(file.path(analysis_dir, "transmission_tree.pdf"))
plot(mcmc_Tree2)
dev.off()

# ==============================================================================
# PART 6: MCMC diagnostics and parameter extraction
# ==============================================================================

cat("Performing MCMC diagnostics...\n")

# Extract MCMC trace
mcmc_trace <- extractParmTrace(mcmc_Tree2)

# Calculate effective sample sizes
eff_sizes <- effectiveSize(mcmc_trace)
cat("MCMC effective sample sizes:\n")
print(eff_sizes)

# Consensus transmission tree
cons <- consTTree(mcmc_Tree2)
ttree <- extractTTree(cons)

# ==============================================================================
# PART 7: Who-infected-whom probability matrix calculation
# ==============================================================================

cat("Computing who-infected-whom probability matrix...\n")

# Calculate transmission probability matrix
mat <- computeMatWIW(mcmc_Tree2)

cat("WIW matrix dimensions:", nrow(mat), "x", ncol(mat), "\n")

# Process matrix data for analysis
matp_m <- reshape2:::melt.matrix(mat, na.rm = TRUE)
matp_m_f <- matp_m[matp_m$Var1 != matp_m$Var2, ]  # Remove diagonal

# Merge with household metadata
metaD <- vw_metadata %>% select(SAMPLE, hh, hh_type) %>% distinct()

# Add infector metadata
names(matp_m_f)[names(matp_m_f) == 'Var1'] <- 'Infector'
names(matp_m_f)[names(matp_m_f) == 'Var2'] <- 'Infectee'

df <- merge(matp_m_f, metaD, by.x = 'Infector', by.y = 'SAMPLE', all = FALSE)
names(df)[names(df) == 'hh'] <- 'hh_Infector'
names(df)[names(df) == 'hh_type'] <- 'hh_type_Infector'

# Add infectee metadata
df1 <- merge(df, metaD, by.x = 'Infectee', by.y = 'SAMPLE', all = FALSE)
names(df1)[names(df1) == 'hh'] <- 'hh_Infectee'
names(df1)[names(df1) == 'hh_type'] <- 'hh_type_Infectee'

# Filter for within-household transmissions
df2 <- df1 %>% 
  filter(hh_Infector == hh_Infectee, 
         hh_Infector != 'hh1', 
         hh_Infectee != 'hh1', 
         hh_type_Infector != 'hh2_noTransmission')

# Determine transmission direction
df3 <- df2 %>%
  group_by(hh_Infectee) %>% 
  mutate(
    TransPhylo_d_r = case_when(
      value == max(value) & value != min(value) ~ "Donor",
      value == min(value) & value != max(value) ~ "Recipient",
      TRUE ~ "Unclear"
    )
  )

# Save transmission direction results
write.csv(df3, file.path(analysis_dir, "transmission_directions.csv"), row.names = FALSE)

# ==============================================================================
# PART 8: Create lineage-specific transmission probability heatmaps
# ==============================================================================

cat("Creating lineage-specific transmission probability heatmaps...\n")

# Define household samples for analysis
list_hh2 <- unique(metad_hh2$SAMPLE)

# Subset matrix to household samples
mat_hh2 <- mat[list_hh2, list_hh2, drop = FALSE]

# Group samples by lineage
grouped_samples <- split(metad_hh2$SAMPLE, metad_hh2$LINEAGE)

# Filter groups to include only existing samples in matrix
grouped_samples <- lapply(grouped_samples, function(samples) {
  valid_samples <- samples[samples %in% rownames(mat_hh2) & samples %in% colnames(mat_hh2)]
  return(valid_samples)
})

# Remove empty groups
grouped_samples <- grouped_samples[lengths(grouped_samples) > 0]

# Create matrices for each lineage group
group_matrices <- lapply(grouped_samples, function(samples) {
  mat_hh2[samples, samples, drop = FALSE]
})

cat("Creating heatmaps for", length(group_matrices), "lineage groups\n")

# Define consistent color scale
mid_value <- 0.5
fill_scale <- scale_fill_gradient2(
  low = 'orchid1', 
  mid = "white", 
  high = "turquoise2",
  midpoint = mid_value, 
  limits = c(0, 1),
  breaks = seq(0, 1, by = 0.2), 
  labels = as.character(seq(0, 1, by = 0.2))
)

# Generate heatmap for each lineage group
for (group in names(group_matrices)) {
  mat_lin <- group_matrices[[group]]
  
  # Skip if matrix is too small
  if (nrow(mat_lin) < 2 || ncol(mat_lin) < 2) {
    cat("Skipping", group, "- insufficient samples\n")
    next
  }
  
  cat("Creating heatmap for lineage:", group, "(", nrow(mat_lin), "x", ncol(mat_lin), ")\n")
  
  # Convert matrix to long format
  mat_long <- melt(mat_lin)
  colnames(mat_long) <- c("Row", "Column", "Value")
  
  # Create heatmap
  plot <- ggplot(mat_long, aes(x = Column, y = Row, fill = Value)) +
    geom_tile() +
    fill_scale +
    theme_minimal() +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 8, angle = 90),
      axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 8),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    ) +
    # Gray diagonal (self-infection impossible)
    geom_tile(data = mat_long[mat_long$Row == mat_long$Column, ], 
              aes(x = Column, y = Row), fill = "gray41") +
    labs(title = paste("Transmission Probabilities:", group)) +
    scale_x_discrete(labels = colnames(mat_lin)) +
    scale_y_discrete(labels = rownames(mat_lin))
  
  # Save individual lineage heatmap
  filename <- file.path(output_dir, paste0("SupplementaryFigure6_", group, "_transmission_matrix.pdf"))
  ggsave(plot = plot, height = 8, width = 10, filename = filename)
  
  cat("Saved heatmap:", filename, "\n")
}

# ==============================================================================
# PART 9: Summary statistics and interpretation
# ==============================================================================

cat("\n=== SUPPLEMENTARY FIGURE 6 ANALYSIS SUMMARY ===\n")

# Dating analysis summary
cat("BactDating analysis:\n")
cat("- Root-to-tip R²:", round(r$R2, 4), "\n")
cat("- Temporal signal strength:", ifelse(r$R2 > 0.3, "Strong", ifelse(r$R2 > 0.1, "Moderate", "Weak")), "\n")

# TransPhylo analysis summary
cat("\nTransPhylo analysis:\n")
cat("- MCMC iterations:", 1000000, "\n")
cat("- Within-household transmission probability (π):", 0.99, "\n")
cat("- Generation time mean:", 0.01, "\n")
cat("- Generation time std:", 0.009, "\n")

# Matrix analysis summary
cat("\nTransmission probability matrix:\n")
cat("- Matrix size:", nrow(mat), "x", ncol(mat), "\n")
cat("- Household transmission pairs analyzed:", nrow(df2), "\n")
cat("- Clear transmission directions identified:", sum(df3$TransPhylo_d_r %in% c("Donor", "Recipient")), "\n")

# Lineage analysis summary
cat("\nLineage-specific analysis:\n")
cat("- Total lineage groups:", length(grouped_samples), "\n")
cat("- Groups with sufficient data for heatmaps:", sum(lengths(grouped_samples) >= 2), "\n")

lineage_counts <- sapply(grouped_samples, length)
cat("- Samples per lineage:\n")
for (i in 1:min(5, length(lineage_counts))) {
  cat("  ", names(lineage_counts)[i], ":", lineage_counts[i], "samples\n")
}

cat("\nMethodological validation:\n")
if (all(eff_sizes > 200)) {
  cat("- MCMC convergence: Good (all ESS > 200)\n")
} else {
  cat("- MCMC convergence: Check required (some ESS < 200)\n")
}

cat("\nOutput files generated:\n")
cat("- Individual lineage transmission probability heatmaps\n")
cat("- BactDating diagnostic plots\n")
cat("- TransPhylo analysis results\n")
cat("- Transmission direction assignments\n")

# ==============================================================================
# End of script
# ==============================================================================
