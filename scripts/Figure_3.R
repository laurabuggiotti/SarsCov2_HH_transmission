# ==============================================================================
# SARS-CoV-2 Phylogenetic Tree Visualization
# Figure 3: Maximum likelihood phylogenetic tree with metadata annotations
# 
# This script creates publication-ready phylogenetic tree visualizations including:
# 1. Fan layout tree with household labels and bootstrap support values
# 2. Circular plot with household type heatmap
# 3. Enhanced circular plot with mutation data and viral load information
#
# Input: 
#   - VW_full_ivar_1k.raxml.support - Maximum likelihood tree with 1000 bootstrap
#   - lineage_VL.csv - Sample metadata including lineage, household, and viral load data
#   - mut - Mutation data (consensus vs minority variants)
# Output: 
#   - Fan layout phylogenetic tree (PDF)
#   - Circular plot with household annotations (PDF)
#   - Enhanced circular plot with mutations and viral load (PDF)
# ==============================================================================

# Load required libraries
library(ape)          # for phylogenetic tree analysis
library(ggtree)       # for tree visualization
library(ggtreeExtra)  # for additional tree annotations
library(dplyr)        # for data manipulation
library(tibble)       # for data frames
library(ggplot2)      # for plotting
library(ggnewscale)   # for multiple color scales

# Set file paths (adjust as needed for your directory structure)
tree_path <- "data/VW_full_ivar_1k.raxml.support"
metadata_path <- "data/lineage_VL.csv"
output_dir <- "figures"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==============================================================================
# PART 1: Load and process phylogenetic tree
# ==============================================================================

cat("Reading maximum likelihood tree with bootstrap support from:", tree_path, "\n")
vw_tree <- read.tree(tree_path)

cat("Tree contains", Ntip(vw_tree), "tips and", Nnode(vw_tree), "internal nodes\n")

# Process bootstrap support values
cat("Processing bootstrap support values...\n")
bs_tibble <- tibble(
  # Internal nodes start after tip nodes
  node = 1:Nnode(vw_tree) + Ntip(vw_tree), 
  # Convert bootstrap values to numeric
  bootstrap = vw_tree$node.label
)

# Clean bootstrap values
bs_tibble$bootstrap <- as.numeric(bs_tibble$bootstrap)
bs_tibble$bootstrap[is.na(bs_tibble$bootstrap)] <- 0

# Convert to dataframe for compatibility
df_bs_tibble <- as.data.frame(bs_tibble)
row.names(df_bs_tibble) <- df_bs_tibble$node
bs <- data.frame("bootstrap" = df_bs_tibble[, c("bootstrap")])

# Clean tree node labels
vw_tree$node.label <- as.numeric(vw_tree$node.label)
vw_tree$node.label[is.na(vw_tree$node.label)] <- 0
tr <- vw_tree

# ==============================================================================
# PART 2: Load and process metadata
# ==============================================================================

cat("Reading metadata from:", metadata_path, "\n")
vw_metadata <- read.csv(metadata_path)

cat("Processing metadata for", nrow(vw_metadata), "samples\n")

# Convert household labels to character
vw_metadata$hh_label <- as.character(vw_metadata$hh_label)

# Set lineage factor levels for consistent ordering
vw_metadata$LINEAGE <- factor(vw_metadata$LINEAGE, 
                             levels = c("BA.2", "BA.2.1", "BA.2.10", "BA.2.37", "BA.2.73", "BM.1.1.1",
                                       "BN.1.7", "BQ.1 ", "BQ.1.1", "BQ.1.1.2", "BQ.1.1.24", "BQ.1.1.5",
                                       "BQ.1.18", "BQ.1.8.2", "CH.1.1", "CH.1.1.1", "CH.1.1.2", "XBB.1", "XBB.1.4",
                                       "XBB.1.5", "XBB.2", "XBC.1", "XBF"))

# Set row names to SampleID for merging with tree
row.names(vw_metadata) <- vw_metadata$SampleID

# Create individual data frames for different annotations
vw_sampID <- data.frame("vw_sampID" = vw_metadata[, c("SampleID")])
lineage <- data.frame("lineage" = vw_metadata[, c("LINEAGE")])
hh <- data.frame('hh' = vw_metadata[, c("hh")])
hh_type <- data.frame('hh_type' = vw_metadata[, c("hh_type")])
voc <- data.frame('VOC_VUI' = vw_metadata[, c("VOC_VUI")])
hh_label <- data.frame('hh_label' = vw_metadata[, c("hh_label")])
VL_cp.ml <- data.frame('VL_cp.ml' = vw_metadata[, c("VL_cp.ml")])

# Set consistent row names
rownames(lineage) <- vw_metadata$SampleID
rownames(vw_sampID) <- vw_metadata$SampleID
rownames(hh) <- vw_metadata$SampleID
rownames(voc) <- vw_metadata$SampleID
rownames(hh_type) <- vw_metadata$SampleID
rownames(hh_label) <- vw_metadata$SampleID
rownames(VL_cp.ml) <- vw_metadata$SampleID

metadata <- vw_metadata

cat("Metadata processing complete\n")
cat("Household types:", paste(unique(vw_metadata$hh_type), collapse = ", "), "\n")

# ==============================================================================
# PART 3: Create base phylogenetic tree plot
# ==============================================================================

cat("Creating base phylogenetic tree visualization...\n")

# Create base tree plot with fan layout and bootstrap support
p <- ggtree(tr, size = 1, layout = "fan", open.angle = 40) + 
  geom_treescale() +
  geom_nodepoint(aes(color = as.numeric(label))) +
  scale_colour_continuous(low = '#8B8000', high = 'darkgreen', 
                         name = "Bootstrap\nSupport") +
  new_scale_colour()

# Add metadata to the plot
p <- p %<+% metadata 

# ==============================================================================
# PART 4: Create detailed fan layout tree
# ==============================================================================

cat("Creating detailed fan layout tree with household annotations...\n")

p1 <- p +
  geom_tiplab(aes(label = hh_label), parse = T, size = 3, hjust = -2, 
              color = 'darkred', fontface = 'bold') +
  geom_tippoint(mapping = aes(color = hh_type), size = 2.5, stroke = 0) +
  scale_color_manual(values = c('darkgray', "darkorange", 'darkred'),
                     name = "Household\nType") +
  new_scale_colour() +
  theme(legend.position = "right",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"))

# Display and save the fan layout plot
print(p1)
# ==============================================================================
# PART 5: Create circular plot with heatmap of household types
# ==============================================================================

cat("Creating circular plot with household type heatmap...\n")

# Create base circular plot
p_circular <- ggtree(tr, size = 1, layout = "fan", open.angle = 40) + 
  geom_treescale() +
  geom_nodepoint(aes(color = as.numeric(label))) +
  scale_colour_continuous(low = '#8B8000', high = 'darkgreen',
                         name = "Bootstrap\nSupport") +
  new_scale_colour()

# Add metadata
p_circular <- p_circular %<+% metadata 

# Add heatmap for household types
pc <- gheatmap(p_circular, hh_type, width = 0.1, font.size = 0) +
  scale_fill_manual(values = c('darkgray', "darkorange", 'darkred'),
                    name = "Household\nType") +
  theme(legend.position = "right")

# Display and save circular plot
print(pc)

# ==============================================================================
# PART 6: Enhanced circular plot with mutations 
# ==============================================================================

cat("Creating enhanced circular plot with mutation data...\n")

# Note: This section requires mutation data ('mut') to be loaded separately
# Uncomment and modify the following code when mutation data is available:

# # Ensure mutation factor levels are set
 if(exists("mut")) {
   mut$mutation <- factor(mut$mutation, levels = c('consVar', 'minorityVar'))
   
   pc1 <- pc + new_scale_fill() +
     geom_fruit(data = mut, geom = geom_bar,
                mapping = aes(y = SampleID, x = n, fill = mutation),
                pwidth = 0.2, offset = .15,
                orientation = "y", 
                stat = "identity") +
     scale_fill_manual(values = c("#B6D0E2", "#6082B6"),
                       name = "Mutation\nType") +
     new_scale_fill() +
     geom_fruit(geom = geom_bar,
                mapping = aes(y = SampleID, x = VL_cp.ml, fill = VL_cp.ml),
                pwidth = 0.2, offset = .005,
                orientation = "y", 
                stat = "identity") + 
     scale_fill_continuous(low = 'blue', high = 'red',
                          name = "Viral Load\n(log10)",
                          guide = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 4))
   
   print(pc1)
   ggsave(plot = pc1, height = 14, width = 16, 
          filename = file.path(output_dir, "Figure3_enhanced_circular_plot.pdf"))
   
   cat("Enhanced circular plot saved\n")
 } else {
   cat("Warning: Mutation data ('mut') not found. Skipping enhanced circular plot.\n")
   cat("To create Figure 3C, load mutation data and uncomment the relevant code section.\n")
 }

# ==============================================================================
# PART 7: Summary statistics and diagnostics
# ==============================================================================

cat("\n=== PHYLOGENETIC TREE ANALYSIS SUMMARY ===\n")
cat("Tree statistics:\n")
cat("- Number of tips (samples):", Ntip(tr), "\n")
cat("- Number of internal nodes:", Nnode(tr), "\n")
cat("- Tree is rooted:", is.rooted(tr), "\n")
cat("- Tree is binary:", is.binary(tr), "\n")

cat("\nBootstrap support summary:\n")
bootstrap_values <- as.numeric(tr$node.label[!is.na(tr$node.label)])
if(length(bootstrap_values) > 0) {
  cat("- Mean bootstrap support:", round(mean(bootstrap_values), 2), "\n")
  cat("- Median bootstrap support:", round(median(bootstrap_values), 2), "\n")
  cat("- Nodes with >70% support:", sum(bootstrap_values > 70), "/", length(bootstrap_values), "\n")
  cat("- Nodes with >95% support:", sum(bootstrap_values > 95), "/", length(bootstrap_values), "\n")
}

cat("\nSample metadata summary:\n")
cat("- Total samples with metadata:", nrow(vw_metadata), "\n")
cat("- Household types:", paste(names(table(vw_metadata$hh_type)), collapse = ", "), "\n")
cat("- Household type counts:\n")
print(table(vw_metadata$hh_type))

if("VL_cp.ml" %in% colnames(vw_metadata)) {
  vl_summary <- summary(vw_metadata$VL_cp.ml[!is.na(vw_metadata$VL_cp.ml)])
  cat("\nViral load summary (log10 copies/ml):\n")
  print(vl_summary)
}

cat("\nOutput files generated:\n")
cat("- Figure3A_fan_layout_tree.pdf\n")
cat("- Figure3B_circular_plot.pdf\n")
cat("- Figure3C_enhanced_circular_plot.pdf (if mutation data available)\n")

# ==============================================================================
# End of script
# ==============================================================================
