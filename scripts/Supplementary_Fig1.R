# ==============================================================================
# SARS-CoV-2 Transmission Bottleneck Analysis - Supplementary Material
# Supplementary Figure 1: Density distribution of bottleneck size estimates
# 
# This script creates density plots showing the distribution of transmission 
# bottleneck size estimates across different variant frequency thresholds.
# The analysis examines how different minority variant detection thresholds
# affect bottleneck size estimates, providing methodological validation.
#
# Input: 
#   - BB_approx_hh2_final.csv - Bottleneck size estimates with different thresholds
#     (or bb_1 object from main Figure 4 analysis)
# Output: 
#   - Density plot showing bottleneck size distributions by threshold
# ==============================================================================

# Load required libraries
library(ggplot2)      # for plotting
library(ggpubr)       # for publication-ready themes
library(dplyr)        # for data manipulation

# Set file paths (adjust as needed for your directory structure)
bottleneck_path <- "data/BB_approx_hh2_final.csv"
output_dir <- "figures"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==============================================================================
# PART 1: Load and process bottleneck data (if not already loaded)
# ==============================================================================

# Check if bb_1 object exists from previous analysis (Figure 4)
# If not, load and process the data
if (!exists("bb_1")) {
  cat("Loading bottleneck data from:", bottleneck_path, "\n")
  bb <- read.csv(bottleneck_path)
  
  cat("Processing bottleneck data for", nrow(bb), "household transmission pairs\n")
  
  # Log-transform bottleneck sizes for visualization
  bb_1 <- bb %>%
    mutate(
      bb_log10 = log10(BB_size), 
      CI_low_log10 = log10(CI_low), 
      CI_high_log10 = log10(CI_high)
    )
  
  cat("Data processing complete\n")
} else {
  cat("Using existing bb_1 object from previous analysis\n")
}

# Verify required columns exist
required_cols <- c("bb_log10", "bb_varThr")
missing_cols <- setdiff(required_cols, colnames(bb_1))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

cat("Data structure:\n")
str(bb_1)

# ==============================================================================
# PART 2: Explore variant threshold categories
# ==============================================================================

cat("\nVariant threshold analysis:\n")
cat("Unique threshold categories:", paste(unique(bb_1$bb_varThr), collapse = ", "), "\n")
cat("Sample counts by threshold:\n")
print(table(bb_1$bb_varThr))

# Summary statistics by threshold
threshold_summary <- bb_1 %>%
  group_by(bb_varThr) %>%
  summarise(
    n_samples = n(),
    mean_log10 = mean(bb_log10, na.rm = TRUE),
    median_log10 = median(bb_log10, na.rm = TRUE),
    sd_log10 = sd(bb_log10, na.rm = TRUE),
    .groups = 'drop'
  )

cat("\nSummary statistics by variant threshold:\n")
print(threshold_summary)

# ==============================================================================
# PART 3: Create density plot
# ==============================================================================

cat("\nCreating density distribution plot...\n")

plot_density <- ggplot(data = bb_1, aes(x = bb_log10, group = bb_varThr, fill = bb_varThr)) +
  geom_density(adjust = 1.5, alpha = .4) +
  scale_fill_manual(
    values = c("lightblue", 'pink', 'olivedrab2'),
    name = "Variant\nThreshold"
  ) +
  scale_colour_manual(
    values = c("darkblue", 'hotpink', 'olivedrab'),
    name = "Variant\nThreshold"
  ) +
  labs(
    title = "", 
    x = 'Bottleneck approx (Log10)', 
    y = "Density"
  ) +
  theme_pubr() +
  theme(
    legend.position = 'right',
    axis.text.y.left = element_text(size = 10),
    axis.text.x.bottom = element_text(size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.ticks.x = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10, face = "bold"),
    strip.background = element_rect(colour = "black", fill = "black")
  )

# Display the plot
print(plot_density)

# ==============================================================================
# PART 4: Save the plot
# ==============================================================================

output_path <- file.path(output_dir, "SupplementaryFigure1_bottleneck_density.pdf")
ggsave(plot = plot_density, height = 8, width = 14, filename = output_path)

cat("Supplementary Figure 1 saved to:", output_path, "\n")

# Also save as PNG for presentations/web use
png_path <- file.path(output_dir, "SupplementaryFigure1_bottleneck_density.png")
ggsave(plot = plot_density, height = 8, width = 14, filename = png_path, dpi = 300)

cat("PNG version saved to:", png_path, "\n")

# ==============================================================================
# PART 5: Statistical analysis and interpretation
# ==============================================================================

cat("\n=== SUPPLEMENTARY FIGURE 1 ANALYSIS SUMMARY ===\n")

# Overall distribution statistics
cat("Overall bottleneck size distribution (log10 scale):\n")
cat("- Mean:", round(mean(bb_1$bb_log10, na.rm = TRUE), 3), "\n")
cat("- Median:", round(median(bb_1$bb_log10, na.rm = TRUE), 3), "\n")
cat("- Standard deviation:", round(sd(bb_1$bb_log10, na.rm = TRUE), 3), "\n")
cat("- Range:", round(min(bb_1$bb_log10, na.rm = TRUE), 3), "to", 
    round(max(bb_1$bb_log10, na.rm = TRUE), 3), "\n")

# Convert back to original scale for interpretation
original_scale_summary <- bb_1 %>%
  summarise(
    mean_original = mean(10^bb_log10, na.rm = TRUE),
    median_original = median(10^bb_log10, na.rm = TRUE),
    min_original = min(10^bb_log10, na.rm = TRUE),
    max_original = max(10^bb_log10, na.rm = TRUE)
  )

cat("\nOriginal scale interpretation:\n")
cat("- Mean bottleneck size:", round(original_scale_summary$mean_original, 1), "variants\n")
cat("- Median bottleneck size:", round(original_scale_summary$median_original, 1), "variants\n")
cat("- Range:", round(original_scale_summary$min_original, 1), "to", 
    round(original_scale_summary$max_original, 1), "variants\n")

# Statistical tests between thresholds (if multiple thresholds exist)
if (length(unique(bb_1$bb_varThr)) > 1) {
  cat("\nStatistical comparison between variant thresholds:\n")
  
  # Perform Kruskal-Wallis test for differences between groups
  kw_test <- kruskal.test(bb_log10 ~ bb_varThr, data = bb_1)
  cat("Kruskal-Wallis test p-value:", format(kw_test$p.value, scientific = TRUE), "\n")
  
  if (kw_test$p.value < 0.05) {
    cat("Significant differences detected between threshold groups\n")
  } else {
    cat("No significant differences between threshold groups\n")
  }
  
  # Pairwise comparisons if more than 2 groups
  if (length(unique(bb_1$bb_varThr)) > 2) {
    cat("\nPairwise comparisons needed for detailed analysis\n")
  }
}

# Distribution shape analysis
cat("\nDistribution characteristics:\n")
library(moments)  # for skewness and kurtosis
if (require(moments, quietly = TRUE)) {
  skew_value <- skewness(bb_1$bb_log10, na.rm = TRUE)
  kurt_value <- kurtosis(bb_1$bb_log10, na.rm = TRUE)
  cat("- Skewness:", round(skew_value, 3), 
      ifelse(skew_value > 0, "(right-skewed)", "(left-skewed)"), "\n")
  cat("- Kurtosis:", round(kurt_value, 3), 
      ifelse(kurt_value > 3, "(heavy-tailed)", "(light-tailed)"), "\n")
} else {
  cat("- Install 'moments' package for skewness/kurtosis analysis\n")
}

cat("\nMethodological interpretation:\n")
cat("- This supplementary figure validates bottleneck estimates across different variant frequency thresholds\n")
cat("- Consistent distributions suggest robust methodology\n")
cat("- Different colored densities represent different analytical parameters\n")
cat("- Log10 transformation normalizes the wide range of bottleneck size estimates\n")

if (length(unique(bb_1$bb_varThr)) == 3) {
  cat("- Three threshold categories suggest sensitivity analysis across detection limits\n")
}

cat("\nOutput files generated:\n")
cat("- SupplementaryFigure1_bottleneck_density.pdf\n")
cat("- SupplementaryFigure1_bottleneck_density.png\n")

# ==============================================================================
# End of script
# ==============================================================================
