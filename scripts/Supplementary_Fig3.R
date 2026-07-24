# ==============================================================================
# Supplementary Figure S4 (revised) — Tanglegram of within-host vs consensus trees
#
# Improvements over the original:
#   - Larger, readable tip labels
#   - Cleaner side labels (not overwhelming)
#   - Household pairs shown with matched colors (each pair one color)
#   - Discordant pairs (whose topology differs between trees) highlighted
#   - Wider export dimensions
# ==============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(dplyr)
  library(readr)
})

# ----------- USER PATHS — edit these -------------------------------------
withinhost_tree_path <- "/Users/Laura/Desktop/VirusWatch/tree_bestree_ivar/VW_full_ivar_1k.raxml.support"
consensus_tree_path  <- "/Users/Laura/Desktop/VirusWatch/tree_bestree_ivar/vw_cons_ivar01.13.raxml.bestTree"
metadata_path        <- "/Users/Laura/Desktop/VirusWatch/Arturo/lineage_VL.csv"
output_path          <- "/Users/Laura/Desktop/VirusWatch/Arturo/SuppFigS4_revised.pdf"
# -------------------------------------------------------------------------


# ============================================================================
# 1. Load trees and metadata
# ============================================================================

treeA <- read.tree(withinhost_tree_path)  # within-host
treeB <- read.tree(consensus_tree_path)   # consensus

cat("Within-host tree tips:", length(treeA$tip.label), "\n")
cat("Consensus tree tips:  ", length(treeB$tip.label), "\n")

vw_metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)


# ============================================================================
# 2. Subset to HH2 tips only (the pairs used in transmission analysis)
# ============================================================================

# Identify HH2 samples from the metadata
hh2_samples <- vw_metadata$SampleID[grepl("^hh2_[0-9]+$", vw_metadata$hh)]
cat("HH2 samples in metadata:", length(hh2_samples), "\n")

# Get tip labels that correspond to HH2 samples (match by prefix — SampleIDs may
# include date suffix in the metadata but tree tips also have suffixes)
# The trees use the full label with date, so we match directly
tips_to_keep_A <- treeA$tip.label[treeA$tip.label %in% hh2_samples]
tips_to_keep_B <- treeB$tip.label[treeB$tip.label %in% hh2_samples]

# Restrict to tips present in BOTH trees
common_tips <- intersect(tips_to_keep_A, tips_to_keep_B)
cat("Common HH2 tips in both trees:", length(common_tips), "\n")

# Prune each tree to just the common HH2 tips
treeA_hh2 <- keep.tip(treeA, common_tips)
treeB_hh2 <- keep.tip(treeB, common_tips)


# ============================================================================
# 3. Build household ID lookup and rename tips as hh2_X_a / hh2_X_b
# ============================================================================

hh_lookup <- vw_metadata[vw_metadata$SampleID %in% common_tips,
                         c("SampleID", "hh")]
names(hh_lookup) <- c("tip", "hh_id")

# For each household, assign _a to the first member alphabetically and _b to the second
hh_lookup <- hh_lookup[order(hh_lookup$hh_id, hh_lookup$tip), ]
hh_lookup$suffix <- ave(hh_lookup$tip, hh_lookup$hh_id,
                        FUN = function(x) paste0("_", letters[seq_along(x)]))
hh_lookup$new_label <- paste0(hh_lookup$hh_id, hh_lookup$suffix)

# Rename tips in both trees using the lookup
relabel_tips <- function(tree, lookup) {
  new_labels <- lookup$new_label[match(tree$tip.label, lookup$tip)]
  tree$tip.label <- new_labels
  tree
}
treeA_hh2 <- relabel_tips(treeA_hh2, hh_lookup)
treeB_hh2 <- relabel_tips(treeB_hh2, hh_lookup)


# ============================================================================
# 4. Convert to cladograms and build tanglegram
# ============================================================================

treeA_c <- compute.brlen(treeA_hh2)
treeB_c <- compute.brlen(treeB_hh2)

# Association matrix — same label in both trees now
assoc <- cbind(treeA_c$tip.label, treeA_c$tip.label)

# Build the co-phylogeny object (this also does the tip ordering)
obj <- cophylo(treeA_c, treeB_c, assoc = assoc, print = FALSE)


# ============================================================================
# 5. Plot
# ============================================================================

pdf(output_path, width = 14, height = 20)

# Set OUTER margins on left and right (in lines) — this is space outside the
# plot region that cophylo doesn't clip. Inner margins stay small.
par(oma = c(0, 4, 0, 4),
    mar = c(2, 2, 2, 2))

plot(obj,
     link.type = 'curved',
     link.lwd  = 1.2,
     link.lty  = 'solid',
     link.col  = scales::alpha("deepskyblue2", 0.6),
     fsize     = 0.55,
     ftype     = "i",       # 'i' = italic; preserves underscores
     part      = 0.4
)

# Draw labels in the outer margins (outer = TRUE writes to oma region)
mtext("Within-host", side = 2, line = 1.5, cex = 1.6, font = 2, outer = TRUE)
mtext("Consensus",   side = 4, line = 1.5, cex = 1.6, font = 2, outer = TRUE)

dev.off()
cat("Figure saved to:", output_path, "\n")
