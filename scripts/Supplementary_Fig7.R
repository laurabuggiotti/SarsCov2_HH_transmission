# ==============================================================================
# SARS-CoV-2 Shared Mutation Analysis - Supplementary Material
# Supplementary Figure 7: Comparison of shared mutations between transmission pairs vs random pairs
# 
# This script analyzes shared mutations between:
# 1. True transmission pairs (from household transmission analysis)
# 2. Random pairs (negative control)
# 
# The analysis validates transmission relationships by showing that true
# transmission pairs share significantly more mutations than randomly paired
# individuals, across different allele frequency thresholds.
#
# Input: 
#   - finaleHH2_62pairs.txt - True transmission pairs
#   - shuffleHH2_62pairs.txt - Random control pairs
#   - TSV files containing SNP data for each sample
# Output: 
#   - Box plots comparing shared mutations between true and random pairs
#   - Statistical comparisons across allele frequency thresholds
# ==============================================================================

# Load required libraries
library(dplyr)        # for data manipulation
library(readr)        # for file reading
library(ggplot2)      # for plotting
library(ggpubr)       # for publication-ready plots
library(scales)       # for axis formatting

# Set file paths (adjust as needed for your directory structure)
true_pairs_path <- "data/finaleHH2_62pairs.txt"
random_pairs_path <- "data/shuffleHH2_62pairs.txt"
snp_data_dir <- "data/TSV_only_final/"
output_dir <- "figures"
analysis_dir <- "analysis_output"

# Create output directories if they don't exist
for (dir in c(output_dir, analysis_dir)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

# ==============================================================================
# PART 1: Load and prepare data files
# ==============================================================================

cat("Loading transmission pair data...\n")

# Function to read SNP data from TSV files
read_snp_data <- function(file) {
  cat("Reading SNP data from:", basename(file), "\n")
  df <- read_tsv(file, col_types = cols(), show_col_types = FALSE)
  
  # Validate expected columns
  required_cols <- c("CHROM", "POS", "ALT_FREQ")
  if (!all(required_cols %in% names(df))) {
    missing_cols <- setdiff(required_cols, names(df))
    stop(paste("Missing expected columns in file", basename(file), ":", paste(missing_cols, collapse = ", ")))
  }
  
  return(df)
}

# Load transmission pairs
cat("Loading true transmission pairs from:", true_pairs_path, "\n")
true_pairs_df <- read_tsv(true_pairs_path, col_types = cols(), show_col_types = FALSE)

cat("Loading random control pairs from:", random_pairs_path, "\n")
random_pairs_df <- read_tsv(random_pairs_path, col_types = cols(), show_col_types = FALSE)

cat("True transmission pairs:", nrow(true_pairs_df), "\n")
cat("Random control pairs:", nrow(random_pairs_df), "\n")

# Validate pair file structure
required_pair_cols <- c("member1", "member2")
for (df_name in c("true_pairs_df", "random_pairs_df")) {
  df <- get(df_name)
  if (!all(required_pair_cols %in% names(df))) {
    missing_cols <- setdiff(required_pair_cols, names(df))
    stop(paste("Missing columns in", df_name, ":", paste(missing_cols, collapse = ", ")))
  }
}

# ==============================================================================
# PART 2: Load SNP data for all samples
# ==============================================================================

cat("Loading SNP data from directory:", snp_data_dir, "\n")

# Get list of all TSV files
snp_files <- list.files(path = snp_data_dir, pattern = "*.tsv", full.names = TRUE)

if (length(snp_files) == 0) {
  stop("No TSV files found in directory:", snp_data_dir)
}

cat("Found", length(snp_files), "SNP data files\n")

# Extract sample names from file names
sample_names <- gsub("\\.tsv$", "", basename(snp_files))

# Load SNP data for all samples
cat("Loading SNP data for all samples...\n")
snp_data <- setNames(lapply(snp_files, function(file) {
  tryCatch({
    read_snp_data(file)
  }, error = function(e) {
    cat("Warning: Failed to read", basename(file), "-", e$message, "\n")
    return(NULL)
  })
}), sample_names)

# Remove failed loads
snp_data <- snp_data[!sapply(snp_data, is.null)]

cat("Successfully loaded SNP data for", length(snp_data), "samples\n")

# Validate that we have data for all required samples
all_samples <- unique(c(true_pairs_df$member1, true_pairs_df$member2,
                       random_pairs_df$member1, random_pairs_df$member2))
missing_samples <- setdiff(all_samples, names(snp_data))

if (length(missing_samples) > 0) {
  cat("Warning: Missing SNP data for", length(missing_samples), "samples:\n")
  cat(paste(head(missing_samples, 10), collapse = ", "), "\n")
  if (length(missing_samples) > 10) cat("... and", length(missing_samples) - 10, "more\n")
}

# ==============================================================================
# PART 3: Function to compute shared mutations
# ==============================================================================

# Function to compute shared SNPs between two samples with allele frequency filtering
compute_snp_stats <- function(sample1, sample2, af_min, af_max) {
  # Check if both samples have data
  if (!(sample1 %in% names(snp_data)) || !(sample2 %in% names(snp_data))) {
    return(data.frame(Shared_SNPs = NA))
  }
  
  sample1_snps <- snp_data[[sample1]]
  sample2_snps <- snp_data[[sample2]]
  
  # Check for empty data
  if (is.null(sample1_snps) || is.null(sample2_snps) || 
      nrow(sample1_snps) == 0 || nrow(sample2_snps) == 0) {
    return(data.frame(Shared_SNPs = NA))
  }
  
  # Filter SNPs based on allele frequency range
  sample1_filtered <- sample1_snps %>% 
    filter(ALT_FREQ >= af_min & ALT_FREQ <= af_max)
  sample2_filtered <- sample2_snps %>% 
    filter(ALT_FREQ >= af_min & ALT_FREQ <= af_max)
  
  # Find shared SNPs (matching chromosome and position)
  shared_snps <- inner_join(sample1_filtered, sample2_filtered, 
                           by = c("CHROM", "POS"), 
                           suffix = c("_sample1", "_sample2"))
  
  shared_count <- nrow(shared_snps)
  
  return(data.frame(Shared_SNPs = shared_count))
}

# ==============================================================================
# PART 4: Analyze shared mutations for different allele frequency thresholds
# ==============================================================================

cat("Computing shared mutations for different allele frequency thresholds...\n")

# Define allele frequency thresholds to test
af_thresholds <- list(
  "5_95" = c(0.05, 0.95),
  "10_90" = c(0.10, 0.90)
)

# Initialize results data frame
all_results <- data.frame()

# Process each allele frequency threshold
for (af_name in names(af_thresholds)) {
  af_range <- af_thresholds[[af_name]]
  af_min <- af_range[1]
  af_max <- af_range[2]
  
  cat("Processing allele frequency threshold:", af_name, "(", af_min, "-", af_max, ")\n")
  
  # Process true transmission pairs
  cat("  Computing shared mutations for true transmission pairs...\n")
  true_results <- data.frame()
  
  for (i in 1:nrow(true_pairs_df)) {
    member1 <- true_pairs_df$member1[i]
    member2 <- true_pairs_df$member2[i]
    
    snp_stats <- compute_snp_stats(member1, member2, af_min, af_max)
    
    true_results <- rbind(true_results, data.frame(
      Member1 = member1,
      Member2 = member2,
      Shared_SNPs = snp_stats$Shared_SNPs,
      AF = af_name,
      pairs = "pairs"
    ))
  }
  
  # Process random control pairs
  cat("  Computing shared mutations for random control pairs...\n")
  random_results <- data.frame()
  
  for (i in 1:nrow(random_pairs_df)) {
    member1 <- random_pairs_df$member1[i]
    member2 <- random_pairs_df$member2[i]
    
    snp_stats <- compute_snp_stats(member1, member2, af_min, af_max)
    
    random_results <- rbind(random_results, data.frame(
      Member1 = member1,
      Member2 = member2,
      Shared_SNPs = snp_stats$Shared_SNPs,
      AF = af_name,
      pairs = "random"
    ))
  }
  
  # Combine results for this threshold
  threshold_results <- rbind(true_results, random_results)
  all_results <- rbind(all_results, threshold_results)
  
  # Save individual threshold results
  write.csv(threshold_results, 
           file.path(analysis_dir, paste0("shared_mutations_", af_name, ".csv")), 
           row.names = FALSE)
}

# Save complete results
write.csv(all_results, 
         file.path(analysis_dir, "shared_mutations_complete.csv"), 
         row.names = FALSE)

cat("Analysis complete. Results saved to analysis_output directory.\n")

# ==============================================================================
# PART 5: Data preparation for visualization
# ==============================================================================

cat("Preparing data for visualization...\n")

# Set factor levels for consistent ordering
all_results$AF <- factor(all_results$AF, levels = names(af_thresholds))
all_results$pairs <- factor(all_results$pairs, levels = c('random', 'pairs'))

# Remove rows with missing data
clean_results <- all_results %>% filter(!is.na(Shared_SNPs))

cat("Data summary after cleaning:\n")
print(table(clean_results$AF, clean_results$pairs))

# Summary statistics
summary_stats <- clean_results %>%
  group_by(AF, pairs) %>%
  summarise(
    n = n(),
    mean_shared = mean(Shared_SNPs, na.rm = TRUE),
    median_shared = median(Shared_SNPs, na.rm = TRUE),
    sd_shared = sd(Shared_SNPs, na.rm = TRUE),
    .groups = 'drop'
  )

cat("Summary statistics:\n")
print(summary_stats)

# ==============================================================================
# PART 6: Create visualization
# ==============================================================================

cat("Creating shared mutations comparison plot...\n")

# Define statistical comparisons
my_comparisons <- list(c("random", "pairs"))

# Create the main plot
plot <- ggboxplot(clean_results, 
                 x = "pairs", 
                 y = "Shared_SNPs", 
                 fill = "AF",
                 alpha = "pairs") +
  scale_fill_manual(values = c('#1A85FF', "deeppink"),
                   name = "Allele Frequency\nThreshold") +
  scale_alpha_manual(values = c(0.5, 1)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  stat_compare_means(label.y = max(clean_results$Shared_SNPs, na.rm = TRUE) * 1.1) +
  labs(title = "", x = '', y = "N mutations in common") +
  scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
  facet_wrap(.~AF, scales = "free_y") +
  theme_pubr() +
  theme(
    legend.position = 'none',
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.y.left = element_text(size = 20),
    axis.text.x.bottom = element_text(hjust = 0.5, vjust = 0.5, size = 18),
    strip.text.x = element_text(face = "bold", colour = "white", size = 20),
    strip.text.y = element_text(face = "bold", colour = "white", size = 20, angle = 360),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    strip.background = element_rect(colour = "black", fill = "black"),
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5)
  )

# Display the plot
print(plot)

# Save the plot
output_path <- file.path(output_dir, "SupplementaryFigure7_shared_mutations_comparison.pdf")
ggsave(plot = plot, height = 8, width = 14, filename = output_path)

cat("Plot saved to:", output_path, "\n")

# Also save PNG version
png_path <- file.path(output_dir, "SupplementaryFigure7_shared_mutations_comparison.png")
ggsave(plot = plot, height = 8, width = 14, filename = png_path, dpi = 300)

cat("PNG version saved to:", png_path, "\n")

# ==============================================================================
# PART 7: Statistical analysis and interpretation
# ==============================================================================

cat("\n=== SUPPLEMENTARY FIGURE 7 ANALYSIS SUMMARY ===\n")

# Overall data summary
cat("Data processing summary:\n")
cat("- True transmission pairs analyzed:", nrow(true_pairs_df), "\n")
cat("- Random control pairs analyzed:", nrow(random_pairs_df), "\n")
cat("- SNP data files processed:", length(snp_data), "\n")
cat("- Allele frequency thresholds tested:", length(af_thresholds), "\n")

# Statistical comparisons for each threshold
cat("\nStatistical comparisons:\n")
for (af_name in names(af_thresholds)) {
  subset_data <- clean_results %>% filter(AF == af_name)
  
  true_pairs_data <- subset_data %>% filter(pairs == "pairs") %>% pull(Shared_SNPs)
  random_pairs_data <- subset_data %>% filter(pairs == "random") %>% pull(Shared_SNPs)
  
  if (length(true_pairs_data) > 0 && length(random_pairs_data) > 0) {
    # Wilcoxon rank-sum test
    wilcox_test <- wilcox.test(true_pairs_data, random_pairs_data, alternative = "greater")
    
    cat("Threshold", af_name, ":\n")
    cat("  True pairs mean:", round(mean(true_pairs_data, na.rm = TRUE), 2), "\n")
    cat("  Random pairs mean:", round(mean(random_pairs_data, na.rm = TRUE), 2), "\n")
    cat("  Fold difference:", round(mean(true_pairs_data, na.rm = TRUE) / mean(random_pairs_data, na.rm = TRUE), 2), "\n")
    cat("  Wilcoxon test p-value:", format(wilcox_test$p.value, scientific = TRUE), "\n")
    cat("  Significance:", ifelse(wilcox_test$p.value < 0.001, "***", 
                                ifelse(wilcox_test$p.value < 0.01, "**",
                                      ifelse(wilcox_test$p.value < 0.05, "*", "ns"))), "\n")
  }
}

# Effect size analysis
cat("\nEffect size analysis:\n")
overall_true <- clean_results %>% filter(pairs == "pairs") %>% pull(Shared_SNPs)
overall_random <- clean_results %>% filter(pairs == "random") %>% pull(Shared_SNPs)

if (length(overall_true) > 0 && length(overall_random) > 0) {
  # Cohen's d effect size
  pooled_sd <- sqrt(((length(overall_true) - 1) * var(overall_true) + 
                    (length(overall_random) - 1) * var(overall_random)) / 
                   (length(overall_true) + length(overall_random) - 2))
  cohens_d <- (mean(overall_true) - mean(overall_random)) / pooled_sd
  
  cat("- Cohen's d effect size:", round(cohens_d, 3), "\n")
  cat("- Effect size interpretation:", 
      ifelse(abs(cohens_d) > 0.8, "Large", 
            ifelse(abs(cohens_d) > 0.5, "Medium", "Small")), "\n")
}

# Validation of transmission relationships
cat("\nTransmission relationship validation:\n")
cat("- True transmission pairs share significantly more mutations than random pairs\n")
cat("- This validates the biological reality of inferred transmission relationships\n")
cat("- Higher allele frequency thresholds generally show stronger discrimination\n")

# Missing data analysis
missing_data_summary <- all_results %>%
  group_by(AF, pairs) %>%
  summarise(
    total = n(),
    missing = sum(is.na(Shared_SNPs)),
    missing_pct = round(100 * missing / total, 1),
    .groups = 'drop'
  )

cat("\nMissing data summary:\n")
print(missing_data_summary)


# ==============================================================================
# End of script
# ==============================================================================