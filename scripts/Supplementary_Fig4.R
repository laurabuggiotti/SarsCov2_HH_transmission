# ==============================================================================
# SARS-CoV-2 Phylogenetic Tree Comparison - Supplementary Material
# Supplementary Figure 4: Tanglegram comparing phylogenetic trees
# 
# This script creates a tanglegram visualization comparing two phylogenetic trees:
# 1. Tree A: Full tree with bootstrap support (maximum likelihood)
# 2. Tree B: Consensus tree with different parameters
# 
# The tanglegram shows how the same samples cluster differently between trees,
# providing insight into phylogenetic reconstruction consistency and the
# impact of different analytical parameters on tree topology.
#
# Input: 
#   - VW_full_ivar_1k.raxml.support - Full tree with bootstrap support
#   - vw_cons_ivar01.13.raxml.bestTree - Consensus tree
# Output: 
#   - Tanglegram visualization showing tree comparison with connecting lines
# ==============================================================================

# Load required libraries
library(ape)          # for phylogenetic tree analysis
library(phytools)     # for tanglegram and tree comparison functions

# Set file paths (adjust as needed for your directory structure)
tree_full_path <- "data/VW_full_ivar_1k.raxml.support"
tree_consensus_path <- "data/vw_cons_ivar01.13.raxml.bestTree"
output_dir <- "figures"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==============================================================================
# PART 1: Load phylogenetic trees
# ==============================================================================

cat("Loading phylogenetic trees...\n")

# Load Tree A: Full tree with bootstrap support
cat("Reading full tree from:", tree_full_path, "\n")
treeA <- read.tree(tree_full_path)

# Load Tree B: Consensus tree
cat("Reading consensus tree from:", tree_consensus_path, "\n")
treeB <- read.tree(tree_consensus_path)

cat("Tree A (full) contains", Ntip(treeA), "tips and", Nnode(treeA), "internal nodes\n")
cat("Tree B (consensus) contains", Ntip(treeB), "tips and", Nnode(treeB), "internal nodes\n")

# ==============================================================================
# PART 2: Process and prepare trees for comparison
# ==============================================================================

cat("Processing trees for comparison...\n")

# Convert trees to cladograms (equal branch lengths for cleaner visualization)
cat("Converting to cladograms with equal branch lengths...\n")
treeA <- compute.brlen(treeA)
treeB <- compute.brlen(treeB)

# Check if trees are identical
trees_identical <- identical(treeA, treeB)
cat("Trees are identical:", trees_identical, "\n")

# Validate tip labels match between trees
treeA_tips <- sort(treeA$tip.label)
treeB_tips <- sort(treeB$tip.label)
matching_tips <- identical(treeA_tips, treeB_tips)
cat("Tip labels match between trees:", matching_tips, "\n")

if (!matching_tips) {
  # Identify differences in tip labels
  only_in_A <- setdiff(treeA_tips, treeB_tips)
  only_in_B <- setdiff(treeB_tips, treeA_tips)
  
  if (length(only_in_A) > 0) {
    cat("Tips only in Tree A:", length(only_in_A), "samples\n")
    cat("First few:", head(only_in_A, 5), "\n")
  }
  
  if (length(only_in_B) > 0) {
    cat("Tips only in Tree B:", length(only_in_B), "samples\n")
    cat("First few:", head(only_in_B, 5), "\n")
  }
  
  # Use common tips only
  common_tips <- intersect(treeA_tips, treeB_tips)
  cat("Using", length(common_tips), "common tips for comparison\n")
  
  # Prune trees to common tips
  treeA <- keep.tip(treeA, common_tips)
  treeB <- keep.tip(treeB, common_tips)
}

# ==============================================================================
# PART 3: Create association matrix for tanglegram
# ==============================================================================

cat("Creating association matrix for tanglegram...\n")

# Create association matrix linking corresponding tips between trees
association <- cbind(treeB$tip.label, treeA$tip.label)
colnames(association) <- c("Tree_B", "Tree_A")

cat("Association matrix dimensions:", nrow(association), "x", ncol(association), "\n")
cat("Sample association matrix (first 5 rows):\n")
print(head(association, 5))

# ==============================================================================
# PART 4: Generate cophylogenetic object and tanglegram
# ==============================================================================

cat("Generating cophylogenetic analysis object...\n")

# Create cophylogenetic object for tanglegram
obj <- cophylo(treeA, treeB, assoc = association, print = TRUE)

cat("Cophylogenetic object created successfully\n")
cat("Object summary:\n")
print(obj)

# ==============================================================================
# PART 5: Create and save tanglegram visualization
# ==============================================================================

cat("Creating tanglegram visualization...\n")

# Set up PDF output
pdf_path <- file.path(output_dir, "SupplementaryFigure4_tanglegram.pdf")
pdf(pdf_path, width = 16, height = 12)

# Create tanglegram plot
plot(obj, 
     link.type = 'curved',           # Curved connecting lines
     link.lwd = 1,                   # Line width
     link.lty = 'solid',             # Solid lines
     link.col = make.transparent('blue', 0.25),  # Semi-transparent blue lines
     fsize = 0.2,                    # Font size for tip labels
     pts = FALSE,                    # No points at tips
     ftype = "i")                    # Italic font for tip labels

# Add title and labels
title(main = "Tanglegram: Phylogenetic Tree Comparison", 
      sub = "Left: Full tree with bootstrap | Right: Consensus tree",
      cex.main = 1.2, cex.sub = 1.0)

# Close PDF device
dev.off()

cat("Tanglegram saved to:", pdf_path, "\n")

# Also create PNG version for presentations
png_path <- file.path(output_dir, "SupplementaryFigure4_tanglegram.png")
png(png_path, width = 1600, height = 1200, res = 150)

plot(obj, 
     link.type = 'curved',
     link.lwd = 1,
     link.lty = 'solid',
     link.col = make.transparent('blue', 0.25),
     fsize = 0.2,
     pts = FALSE,
     ftype = "i")

title(main = "Tanglegram: Phylogenetic Tree Comparison", 
      sub = "Left: Full tree with bootstrap | Right: Consensus tree",
      cex.main = 1.2, cex.sub = 1.0)

dev.off()

cat("PNG version saved to:", png_path, "\n")

# ==============================================================================
# PART 6: Tree comparison statistics and analysis
# ==============================================================================

cat("\n=== SUPPLEMENTARY FIGURE 4 ANALYSIS SUMMARY ===\n")

# Basic tree statistics
cat("Tree comparison statistics:\n")
cat("- Tree A (full) tips:", Ntip(treeA), "\n")
cat("- Tree A (full) nodes:", Nnode(treeA), "\n")
cat("- Tree B (consensus) tips:", Ntip(treeB), "\n")
cat("- Tree B (consensus) nodes:", Nnode(treeB), "\n")
cat("- Trees identical:", trees_identical, "\n")

# Robinson-Foulds distance
if (require(phangorn, quietly = TRUE)) {
  rf_distance <- RF.dist(treeA, treeB)
  cat("- Robinson-Foulds distance:", rf_distance, "\n")
  
  # Normalized RF distance
  max_rf <- 2 * (Ntip(treeA) - 3)  # Maximum possible RF distance
  normalized_rf <- rf_distance / max_rf
  cat("- Normalized RF distance:", round(normalized_rf, 3), "\n")
} else {
  cat("- Install 'phangorn' package for Robinson-Foulds distance calculation\n")
}

# Tree comparison metrics using ape
if (Ntip(treeA) == Ntip(treeB)) {
  # Cophenetic correlation
  dist_A <- cophenetic.phylo(treeA)
  dist_B <- cophenetic.phylo(treeB)
  
  # Calculate correlation between distance matrices
  cophenetic_cor <- cor(as.vector(dist_A), as.vector(dist_B), use = "complete.obs")
  cat("- Cophenetic correlation:", round(cophenetic_cor, 3), "\n")
  
  if (cophenetic_cor > 0.8) {
    cat("  -> High correlation: trees are topologically similar\n")
  } else if (cophenetic_cor > 0.6) {
    cat("  -> Moderate correlation: some topological differences\n")
  } else {
    cat("  -> Low correlation: substantial topological differences\n")
  }
}

# Tip label analysis
common_tips <- intersect(treeA$tip.label, treeB$tip.label)
cat("\nTip label analysis:\n")
cat("- Tips in both trees:", length(common_tips), "\n")
cat("- Tips only in Tree A:", length(setdiff(treeA$tip.label, treeB$tip.label)), "\n")
cat("- Tips only in Tree B:", length(setdiff(treeB$tip.label, treeA$tip.label)), "\n")

# ==============================================================================
# End of script
# ==============================================================================
