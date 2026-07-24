# ==============================================================================
# Supplementary Figure 5 (revised) — Global phylogeny with Virus Watch samples
# in context of Nextstrain UK background sequences.
# ==============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
})

# ----------- USER PATHS — edit these ------------------------------------
tree_path     <- "~/UCL Dropbox/Laura Buggiotti/Mac/Desktop/VirusWatch/Arturo/genbank/tree/Nextstrain_HighCov_VW_consREF_tree.treefile"
metadata_path <- "~/UCL Dropbox/Laura Buggiotti/Mac/Desktop/VirusWatch/Arturo/genbank/NextstrainHighCov_VW_metadata.csv"
output_path   <- "/Users/Laura/Desktop/VirusWatch/Arturo/SuppFig5_revised.pdf"
# ------------------------------------------------------------------------


# ============================================================================
# 1. Load tree
# ============================================================================

tree <- read.tree(tree_path)
cat("Tree tips:", length(tree$tip.label), "\n")
cat("First 5 tip labels:", paste(head(tree$tip.label, 5), collapse = ", "), "\n\n")


# ============================================================================
# 2. Load metadata (auto-detect tab vs comma)
# ============================================================================

# Force CSV reader (your metadata is comma-separated)
vw_metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)

# If it somehow came in with a single mashed column, try tab separator
if (ncol(vw_metadata) < 3) {
  vw_metadata <- read.table(metadata_path, sep = "\t", header = TRUE,
                            stringsAsFactors = FALSE, quote = "", fill = TRUE)
}

cat("Metadata columns:", paste(names(vw_metadata), collapse = ", "), "\n")
cat("Metadata rows:", nrow(vw_metadata), "\n\n")


# ============================================================================
# 3. Build tip-annotation data frame
# ============================================================================

tip_df <- data.frame(
  label = tree$tip.label,
  category = "Nextstrain background",
  stringsAsFactors = FALSE
)

# For each metadata row, find the corresponding tip and categorise it
for (i in seq_len(nrow(vw_metadata))) {
  sid <- vw_metadata$SampleID[i]
  hh_type <- vw_metadata$hh_type[i]
  if (is.na(sid) || sid == "") next
  
  # Try exact match
  match_idx <- which(tip_df$label == sid)
  # Fall back to matching prefix (before any date suffix)
  if (length(match_idx) == 0) {
    stripped_tips <- sub("-[0-9]{8}$", "", tip_df$label)
    stripped_sid <- sub("-[0-9]{8}$", "", sid)
    match_idx <- which(stripped_tips == stripped_sid)
  }
  
  if (length(match_idx) >= 1) {
    match_idx <- match_idx[1]
    if (hh_type == "hh1") {
      tip_df$category[match_idx] <- "HH1: single-positive household"
    } else if (hh_type == "hh2") {
      tip_df$category[match_idx] <- "HH2: multi-positive (transmission)"
    } else if (hh_type == "hh2_noTransmission") {
      tip_df$category[match_idx] <- "HH2: multi-positive (no transmission)"
    }
  }
}

cat("Tips assigned to each category:\n")
print(table(tip_df$category))
cat("\n")

tip_df$category <- factor(tip_df$category,
                          levels = c("Nextstrain background",
                                     "HH1: single-positive household",
                                     "HH2: multi-positive (transmission)",
                                     "HH2: multi-positive (no transmission)"))


# ============================================================================
# 4. Build the plot
# ============================================================================

category_colors <- c(
  "Nextstrain background"                   = "grey80",
  "HH1: single-positive household"          = "#1E90FF",
  "HH2: multi-positive (transmission)"      = "#8A2BE2",
  "HH2: multi-positive (no transmission)"   = "#FF8C00"
)

tip_df$point_size <- ifelse(tip_df$category == "Nextstrain background", 0.4, 2.5)
tip_df$point_alpha <- ifelse(tip_df$category == "Nextstrain background", 0.4, 0.95)

p <- ggtree(tree, layout = "fan", open.angle = 15, size = 0.12) %<+% tip_df +
  geom_tippoint(aes(color = category, size = point_size, alpha = point_alpha)) +
  scale_color_manual(values = category_colors,
                     name = "Sample category",
                     guide = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  scale_size_identity() +
  scale_alpha_identity() +
  theme(legend.position = c(0.90, 0.55),   # inside the plot, near top-right
        legend.justification = c(0, 0.5),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12),
        legend.key = element_blank(),
        legend.background = element_rect(fill = "white", colour = NA),
        plot.margin = ggplot2::margin(5, 5, 5, 5))


# ============================================================================
# 5. Save at large canvas
# ============================================================================

ggsave(plot = p,
       filename = output_path,
       width = 18,
       height = 18,
       dpi = 300)

cat("\nFigure saved to:", output_path, "\n")
