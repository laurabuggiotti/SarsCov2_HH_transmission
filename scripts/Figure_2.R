# ==============================================================================
# SARS-CoV-2 Phylogenetic Distance Analysis
# Figure 2: Within-host diversity - pairwise branch lengths by household type
# 
# This script calculates pairwise phylogenetic distances from a phylogenetic tree
# and compares genetic distances between samples from the same household vs 
# different households. Creates violin plots showing within-host diversity.
#
# Input: 
#   - VW_full_ivar_1k.raxml.support - phylogenetic tree file
#   - lineage_VL.csv - metadata with sample and household information
# Output: 
#   - Violin plot comparing pairwise branch lengths by household type
# ==============================================================================

# Load required libraries
library(ape)          # for phylogenetic tree analysis
library(ggtree)       # for tree visualization
library(dplyr)        # for data manipulation
library(stringr)      # for string operations
library(reshape2)     # for data reshaping
library(ggplot2)      # for plotting
library(ggpubr)       # for publication-ready plots

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

cat("Reading phylogenetic tree from:", tree_path, "\n")
vw_tree_full <- read.tree(tree_path)

# Optional: visualize tree with branch lengths (uncomment to view)
# ggtree(vw_tree_full) + geom_text2(aes(label = branch.length), hjust = -.3, size = 1)

cat("Tree contains", length(vw_tree_full$tip.label), "samples\n")
cat("Sample of tip labels:", head(unique(vw_tree_full$tip.label)), "\n")

# Calculate pairwise phylogenetic distances (branch lengths)
cat("Calculating pairwise phylogenetic distances...\n")
branch_dist <- cophenetic.phylo(vw_tree_full)

# ==============================================================================
# PART 2: Process distance matrix
# ==============================================================================

# Convert distance matrix to long format
t <- reshape2:::melt.matrix(branch_dist, na.rm = TRUE)

# Extract patient IDs from sample names (assuming format: PatientID-SampleType)
t_mo <- t %>%
  mutate(pat1 = str_split(Var1, "-", simplify = T)[,1], 
         pat2 = str_split(Var2, "-", simplify = T)[,1])

# Filter out self-comparisons and get unique pairs
t_mo1 <- t_mo %>% 
  filter(pat1 != pat2) %>% 
  distinct()

# Create reverse pairs to ensure all combinations
t_mo2 <- t %>%
  mutate(pat1 = str_split(Var2, "-", simplify = T)[,1], 
         pat2 = str_split(Var1, "-", simplify = T)[,1])

t_mo3 <- t_mo2 %>% 
  filter(pat1 != pat2) %>% 
  distinct()

# Combine and deduplicate
t_f <- rbind(t_mo1, t_mo3)
t_f1 <- t_f %>% 
  filter(pat1 != pat2) %>% 
  distinct() %>% 
  select(pat1, pat2, value) %>% 
  distinct()

cat("Total pairwise comparisons:", nrow(t_f1), "\n")

# ==============================================================================
# PART 3: Merge with household metadata
# ==============================================================================

cat("Reading metadata from:", metadata_path, "\n")
vw_metadata <- read.csv(metadata_path)

# Extract relevant household information
metadata <- subset(vw_metadata, select = c('SAMPLE', 'hh', 'hh_extra'))

# Merge distance data with household information for both patients
merged <- merge(t_f1, metadata, by.x = 'pat1', by.y = 'SAMPLE', all = TRUE)
colnames(merged)[which(names(merged) == "hh")] <- "hh_pat1"

merged1 <- merge(merged, metadata, by.x = 'pat2', by.y = 'SAMPLE', all = TRUE)
colnames(merged1)[which(names(merged1) == "hh")] <- "hh_pat2"

# ==============================================================================
# PART 4: Classify household relationships
# ==============================================================================

data_final_full <- merged1 %>%
  mutate(HH_extra = case_when(
    hh_pat1 == hh_pat2 & hh_extra.x == 'hh2_noTransmission' & hh_extra.y == 'hh2_noTransmission' ~ value
  )) %>%
  mutate(HH_final = case_when(
    hh_pat1 == hh_pat2 ~ hh_pat2,
    TRUE ~ 'hh1'
  )) %>% 
  mutate(HH_final = sub("_.*", "", HH_final))

# Standardize household labels
data_final_full$HH_final <- gsub("hh3", "hh2", data_final_full$HH_final)
data_final_full$HH_final <- gsub("hh2", "Same HH", data_final_full$HH_final)
data_final_full$HH_final <- gsub("hh1", "Different HH", data_final_full$HH_final)

# Print summary statistics
cat("\nData summary:\n")
print(str(data_final_full))
cat("Household comparison counts:\n")
print(table(data_final_full$HH_final, useNA = "ifany"))

# ==============================================================================
# PART 5: Create visualizations
# ==============================================================================

# Set up statistical comparisons
my_comparisons <- list(c("Same HH", "Different HH"))

# Create boxplot version (optional)
plot_box <- ggboxplot(data_final_full, x = "HH_final", y = "value", fill = "HH_final",
                      palette = c('#1A85FF', '#D41159')) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  stat_compare_means(label.y = 0.020) +
  labs(title = "Within-host diversity", x = 'Type of households', y = "Pairwise branch length") +
  theme(legend.position = 'none')

# Create main violin plot with log transformation
plot_violin <- data_final_full %>%
  mutate(value_log = log2(value)) %>%
  mutate(HH_extra_log = log2(HH_extra)) %>%
  ggviolin(x = "HH_final", y = "value_log", fill = "HH_final",
           palette = c('#1A85FF', '#D41159')) +
  geom_point(aes(y = HH_extra_log), size = 2, colour = "#FFC20A") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  stat_compare_means(label.y = -4) +
  scale_y_continuous(limits = c(-18, -4), breaks = seq(-18, -4, 2)) +
  labs(title = "Within-host diversity", 
       x = '', 
       y = "Pairwise branch length (log2)") +
  theme(legend.position = 'none',
        plot.title = element_text(size = 18, face = 'bold'),
        axis.title = element_text(size = 14, face = 'bold')) 

# Display the main plot
print(plot_violin)

# Save plots
ggsave(plot = plot_violin, 
       height = 6, 
       width = 8, 
       filename = file.path(output_dir, "Figure2_phylogenetic_distance_violin.pdf"))

ggsave(plot = plot_box, 
       height = 6, 
       width = 8, 
       filename = file.path(output_dir, "Figure2_phylogenetic_distance_boxplot.pdf"))

cat("Plots saved to:", output_dir, "\n")

# ==============================================================================
# PART 6: Summary statistics
# ==============================================================================

cat("\nFinal summary statistics:\n")
summary_stats <- data_final_full %>%
  group_by(HH_final) %>%
  summarise(
    n_comparisons = n(),
    mean_distance = mean(value, na.rm = TRUE),
    median_distance = median(value, na.rm = TRUE),
    sd_distance = sd(value, na.rm = TRUE),
    .groups = 'drop'
  )

print(summary_stats)

# Statistical test results
cat("\nStatistical test (Wilcoxon rank-sum test):\n")
if(sum(!is.na(data_final_full$value)) > 0) {
  same_hh <- data_final_full$value[data_final_full$HH_final == "Same HH" & !is.na(data_final_full$value)]
  diff_hh <- data_final_full$value[data_final_full$HH_final == "Different HH" & !is.na(data_final_full$value)]
  
  if(length(same_hh) > 0 & length(diff_hh) > 0) {
    test_result <- wilcox.test(same_hh, diff_hh)
    cat("p-value =", test_result$p.value, "\n")
  }
}

# ==============================================================================
# End of script
# ==============================================================================
