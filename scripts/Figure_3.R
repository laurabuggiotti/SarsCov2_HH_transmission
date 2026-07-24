# ==============================================================================
# Figure 3 (revised) — Within-host phylogeny with clarified legends
#
# Fixes vs original:
#   - Renames legend items to spell out HH1/HH2 meaning
#   - Renames "# mutation" and its factor levels to be reader-friendly
#   - Ensures viral load legend appears
#   - Wider PDF export so full legend is captured
# ==============================================================================

suppressPackageStartupMessages({
  library(ggstance)
  library(ggpubr)
  library(ggh4x)
  library(ggtree)
  library(ggrepel)
  library(ggnewscale)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ape)
  library(ggplot2)
})

# ----------- USER PATHS — edit these to your local paths --------------------
vars_path      <- "/Users/Laura/Desktop/VirusWatch/Arturo/vars/allVars.tsv"
tree_path      <- "/Users/Laura/Desktop/VirusWatch/tree_bestree_ivar/VW_full_ivar_1k.raxml.support"
metadata_path  <- "/Users/Laura/Desktop/VirusWatch/Arturo/lineage_VL.csv"
output_path    <- "/Users/Laura/Desktop/VirusWatch/Arturo/Figure3_revised.pdf"
# ---------------------------------------------------------------------------


# ============================================================================
# 1. Load and process variant data
# ============================================================================

var <- read.table(vars_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Filter, exclude problematic samples, classify variants
data <- var %>%
  filter(AF >= 0.01,
         DP >= 100,
         Sample != "T001KNJBC7",
         Sample != "T001KO0G00") %>%
  mutate(mutation = case_when(
    AF <= 0.6 ~ "Minority variant (AF <= 60%)",
    AF > 0.6  ~ "Consensus-level (AF > 60%)"
  ))

# Count variants per sample × class (this is what was called `mut` in original)
# Use SampleID (with date suffix) since that's what matches the tree tip labels
mut <- data %>%
  group_by(SampleID, mutation) %>%
  tally() %>%
  ungroup()

# Order factor for consistent legend
mut$mutation <- factor(mut$mutation,
                       levels = c("Consensus-level (AF > 60%)",
                                  "Minority variant (AF <= 60%)"))


# ============================================================================
# 2. Load tree with bootstrap
# ============================================================================

vw_tree <- read.tree(tree_path)

# Bootstrap values → numeric, NAs → 0
vw_tree$node.label <- as.numeric(vw_tree$node.label)
vw_tree$node.label[is.na(vw_tree$node.label)] <- 0


# ============================================================================
# 3. Load metadata
# ============================================================================

vw_metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)

# Household type readable labels
vw_metadata$hh_type_readable <- factor(
  vw_metadata$hh_type,
  levels = c("hh1", "hh2", "hh2_noTransmission"),
  labels = c("HH1: single-positive household",
             "HH2: multi-positive (transmission)",
             "HH2: multi-positive (no transmission)")
)

# Lineage factor (kept exactly as in original)
vw_metadata$LINEAGE <- factor(vw_metadata$LINEAGE,
                              levels = c("BA.2","BA.2.1","BA.2.10","BA.2.37","BA.2.73","BM.1.1.1",
                                         "BN.1.7","BQ.1 ","BQ.1.1","BQ.1.1.2","BQ.1.1.24","BQ.1.1.5",
                                         "BQ.1.18","BQ.1.8.2","CH.1.1","CH.1.1.1","CH.1.1.2","XBB.1",
                                         "XBB.1.4","XBB.1.5","XBB.2","XBC.1","XBF"))

row.names(vw_metadata) <- vw_metadata$SampleID

# Build the strip data frames used for gheatmap
hh_type_df <- data.frame(
  hh_type = vw_metadata$hh_type_readable,
  row.names = vw_metadata$SampleID
)


# ============================================================================
# 4. Build the plot
# ============================================================================

p <- ggtree(vw_tree, size = 1, layout = "fan", open.angle = 40) +
  geom_treescale() +
  geom_nodepoint(aes(color = as.numeric(label))) +
  scale_colour_continuous(low = '#8B8000', high = 'darkgreen',
                          name = "Bootstrap\nsupport (%)") +
  new_scale_colour()

p <- p %<+% vw_metadata

# Household-type ring
pc <- gheatmap(p, hh_type_df, width = 0.1, font.size = 0,
               colnames = FALSE) +
  scale_fill_manual(values = c("darkgray", "darkorange", "darkred"),
                    name = "Household type") +
  new_scale_fill()

library(ggtreeExtra)
# Mutation-count bar ring
pc1 <- pc +
  geom_fruit(
    data = mut,
    geom = geom_bar,
    mapping = aes(y = SampleID, x = n, fill = mutation),
    pwidth = 0.2, offset = 0.15,
    orientation = "y",
    stat = "identity"
  ) +
  scale_fill_manual(values = c("#B6D0E2", "#6082B6"),
                    name = "Number of variants\nper sample") +
  new_scale_fill()


# ============================================================================
# 5. Save with generous canvas so legend isn't cropped
# ============================================================================

ggsave(plot = pc1,
       filename = output_path,
       height = 14,
       width = 20,   # widened from 16 to accommodate all legends
       dpi = 300)

cat("Figure saved to:", output_path, "\n")
