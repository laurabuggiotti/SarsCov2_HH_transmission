# ==============================================================================
# SARS-CoV-2 Transmission Bottleneck Analysis
# Figure 4: Transmission bottleneck size estimates and posterior probabilities
# 
# This script creates a two-panel figure showing:
# Panel A: Estimated transmission bottleneck sizes (Nb) with confidence intervals
# Panel B: Posterior probabilities and proportion of shared variants between households
#
# The analysis examines transmission bottlenecks in household pairs where 
# transmission likely occurred, providing insights into the genetic diversity
# that passes through transmission events.
#
# Input: 
#   - BB_approx_hh2_final.csv - Bottleneck size estimates with confidence intervals
#   - PP_propShared_value.txt - Posterior probabilities and shared variant proportions
# Output: 
#   - Combined two-panel figure showing bottleneck analysis results
# ==============================================================================

# Load required libraries
library(tidyverse)    # for data manipulation and visualization
library(forcats)      # for factor handling
library(cowplot)      # for combining plots
library(ggplot2)      # for plotting

# Set file paths (adjust as needed for your directory structure)
bottleneck_path <- "data/BB_approx_hh2_final.csv"
posterior_path <- "data/PP_propShared_value.txt"
output_dir <- "figures"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==============================================================================
# PART 1: Load and process bottleneck size data
# ==============================================================================

cat("Reading transmission bottleneck data from:", bottleneck_path, "\n")
bb <- read.csv(bottleneck_path)

cat("Bottleneck data contains", nrow(bb), "household transmission pairs\n")
cat("Data structure:\n")
str(bb)

# Log-transform bottleneck sizes and confidence intervals for visualization
cat("Log-transforming bottleneck size estimates...\n")
bb_1 <- bb %>%
  mutate(
    bb_log10 = log10(BB_size), 
    CI_low_log10 = log10(CI_low), 
    CI_high_log10 = log10(CI_high)
  )

# Reorder households by bottleneck size for better visualization
bb_1$hh <- fct_reorder(bb_1$hh, bb_1$bb_log10)

cat("Bottleneck size summary (original scale):\n")
print(summary(bb$BB_size))
cat("Bottleneck size range:", min(bb$BB_size, na.rm = TRUE), "to", max(bb$BB_size, na.rm = TRUE), "\n")

# ==============================================================================
# PART 2: Create bottleneck size plot (Panel A)
# ==============================================================================

cat("Creating transmission bottleneck size plot (Panel A)...\n")

plot_bottleneck <- ggplot(bb_1, aes(x = hh, y = bb_log10, ymin = CI_low_log10, ymax = CI_high_log10)) + 
  geom_point(position = position_dodge(width = 0.5), size = 2.5, alpha = 0.7, colour = 'olivedrab') + 
  geom_linerange(position = position_dodge(width = 0.5), size = 0.5, colour = 'olivedrab') +
  labs(title = "", x = 'HH pairs', y = "Transmission bottleneck size (Nb)") +
  scale_y_continuous(
    breaks = c(1, 2, 3),
    labels = c("10", "100", "1000")
  ) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title.x = element_blank(),  # Remove x-axis title for cleaner look
    axis.title.y = element_text(face = "bold", size = 18),
    axis.text.y.left = element_text(size = 14),
    axis.ticks.x = element_blank(),
    axis.text.x.bottom = element_text(hjust = 0.5, vjust = 0.5, size = 10, angle = 90),
    strip.text.x = element_text(face = "bold", colour = "white", size = 10),
    strip.text.y = element_text(face = "bold", colour = "white", size = 10, angle = 360),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    strip.background = element_rect(colour = "black", fill = "black")
  )

# Display Panel A
print(plot_bottleneck)

# ==============================================================================
# PART 3: Load and process posterior probability data
# ==============================================================================

cat("Reading posterior probability and shared variants data from:", posterior_path, "\n")
df <- read.table(posterior_path, sep = '\t', stringsAsFactors = FALSE, header = TRUE)

cat("Posterior probability data contains", nrow(df), "household pairs\n")
cat("Data structure:\n")
str(df)

# Prepare shared proportions in long format for plotting
df_shared <- df %>%
  pivot_longer(cols = starts_with("prop_shared"),
               names_to = "type", values_to = "value")

# Ensure consistent factor ordering between datasets
hh_levels <- levels(bb_1$hh)
df$hh2_type <- factor(df$hh2_type, levels = hh_levels)
df_shared$hh2_type <- factor(df_shared$hh2_type, levels = hh_levels)

# Filter out missing data
df_shared <- df_shared %>% filter(hh2_type != 'NA')
df <- df %>% filter(hh2_type != 'NA')

cat("Final datasets after filtering:\n")
cat("- Shared variants data:", nrow(df_shared), "observations\n")
cat("- Posterior probability data:", nrow(df), "household pairs\n")

# ==============================================================================
# PART 4: Create posterior probability and shared variants plot (Panel B)
# ==============================================================================

cat("Creating posterior probability and shared variants plot (Panel B)...\n")

plot_posterior <- ggplot() +
  # Lines connecting prop_shared_a and prop_shared_b for each household
  geom_line(data = df_shared,
            aes(x = hh2_type, y = value, group = hh2_type),
            color = "steelblue", size = 0.5) +
  
  # Points for proportion of shared variants (both directions)
  geom_point(data = df_shared,
             aes(x = hh2_type, y = value),
             color = "steelblue", size = 3) +
  
  # Triangular points for posterior probability (separate layer)
  geom_point(data = df,
             aes(x = hh2_type, y = posteriorProb),
             color = "darkorange", size = 3.5, shape = 17) +
  
  scale_shape_manual(values = c("prop_shared_a" = 16, "prop_shared_b" = 1)) +
  ylim(0, 1) +
  labs(y = "Probability / Proportion") +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y.left = element_text(size = 14),
    axis.ticks.x = element_blank(),
    axis.text.x.bottom = element_text(hjust = 0.5, vjust = 0.5, size = 10, angle = 90),
    strip.text.x = element_text(face = "bold", colour = "white", size = 10),
    strip.text.y = element_text(face = "bold", colour = "white", size = 10, angle = 360),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    strip.background = element_rect(colour = "black", fill = "black")
  )

# Display Panel B
print(plot_posterior)

# ==============================================================================
# PART 5: Combine plots into final figure
# ==============================================================================

cat("Combining plots into final two-panel figure...\n")

plot_final <- cowplot::plot_grid(
  plot_bottleneck, plot_posterior, 
  ncol = 1, 
  rel_heights = c(1, 0.5),
  align = 'v', 
  axis = 'rlbt',
  labels = c("a", "b"),
  label_size = 20
) 

# Display final combined plot
print(plot_final)

# Save the final combined figure
output_path <- file.path(output_dir, "Figure4_transmission_bottleneck_analysis.pdf")
ggsave(plot = plot_final, height = 10, width = 16, filename = output_path)

cat("Final figure saved to:", output_path, "\n")

# Also save individual panels for flexibility
ggsave(plot = plot_bottleneck, height = 8, width = 14, 
       filename = file.path(output_dir, "Figure4A_bottleneck_sizes.pdf"))

ggsave(plot = plot_posterior, height = 4, width = 14, 
       filename = file.path(output_dir, "Figure4B_posterior_probabilities.pdf"))

# ==============================================================================
# PART 6: Summary statistics and interpretation
# ==============================================================================

cat("\n=== TRANSMISSION BOTTLENECK ANALYSIS SUMMARY ===\n")

# Bottleneck size statistics
cat("Bottleneck size analysis:\n")
cat("- Number of household transmission pairs analyzed:", nrow(bb_1), "\n")
cat("- Median bottleneck size:", round(median(bb$BB_size, na.rm = TRUE), 1), "variants\n")
cat("- Mean bottleneck size:", round(mean(bb$BB_size, na.rm = TRUE), 1), "variants\n")
cat("- Range:", round(min(bb$BB_size, na.rm = TRUE), 1), "to", round(max(bb$BB_size, na.rm = TRUE), 1), "variants\n")

# Categorize bottleneck sizes
small_bottleneck <- sum(bb$BB_size < 10, na.rm = TRUE)
medium_bottleneck <- sum(bb$BB_size >= 10 & bb$BB_size < 100, na.rm = TRUE)
large_bottleneck <- sum(bb$BB_size >= 100, na.rm = TRUE)

cat("Bottleneck size distribution:\n")
cat("- Small (<10 variants):", small_bottleneck, "pairs\n")
cat("- Medium (10-99 variants):", medium_bottleneck, "pairs\n")
cat("- Large (≥100 variants):", large_bottleneck, "pairs\n")

# Posterior probability statistics
if(nrow(df) > 0 && "posteriorProb" %in% colnames(df)) {
  cat("\nPosterior probability analysis:\n")
  cat("- Mean posterior probability:", round(mean(df$posteriorProb, na.rm = TRUE), 3), "\n")
  cat("- Median posterior probability:", round(median(df$posteriorProb, na.rm = TRUE), 3), "\n")
  high_confidence <- sum(df$posteriorProb > 0.8, na.rm = TRUE)
  cat("- High confidence pairs (>0.8):", high_confidence, "/", nrow(df), "\n")
}

# Shared variants statistics
if(nrow(df_shared) > 0) {
  shared_summary <- df_shared %>%
    group_by(type) %>%
    summarise(
      mean_shared = mean(value, na.rm = TRUE),
      median_shared = median(value, na.rm = TRUE),
      .groups = 'drop'
    )
  
  cat("\nShared variants analysis:\n")
  print(shared_summary)
}

cat("\nOutput files generated:\n")
cat("- Figure4_transmission_bottleneck_analysis.pdf (combined figure)\n")
cat("- Figure4A_bottleneck_sizes.pdf (Panel A only)\n")
cat("- Figure4B_posterior_probabilities.pdf (Panel B only)\n")

cat("\nInterpretation notes:\n")
cat("- Panel A shows estimated transmission bottleneck sizes with 95% confidence intervals\n")
cat("- Panel B shows posterior probabilities (orange triangles) and shared variant proportions (blue circles)\n")
cat("- Lines in Panel B connect bidirectional sharing proportions between household pairs\n")
cat("- Smaller bottlenecks suggest more constrained transmission events\n")

# ==============================================================================
# End of script
# ==============================================================================
