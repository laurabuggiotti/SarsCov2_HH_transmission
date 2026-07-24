# ==============================================================================
# REPLACEMENT for Supplementary Figure 7
# Shared minority variants: household transmission pairs vs random pairs
#
# Changes from the original:
#   1. Uses the Sobel Leonard validated 3-50% minor allele frequency window
#      (the original used 5-95% and 10-95%, which include consensus-level
#      variants that two household members share simply by inheriting the
#      same lineage)
#   2. Applies the cohort-wide artifact blacklist (13 recurrent positions)
#      to remove amplicon-sequencing artifacts before counting shared variants
#   3. Matches variants by CHROM:POS:REF:ALT, not just CHROM:POS (so two
#      samples with different alt alleles at the same site aren't miscounted
#      as sharing)
#   4. Generates a control plot showing the BEFORE / AFTER artifact filter
#      comparison, demonstrating that filtering increases the household/random
#      gap rather than narrowing it
#
# Input:
#   - /SAN/.../bottleneck_out/pairs.json (household pairs)
#   - /SAN/.../bottleneck_out/random_pairs.json (random pairs)
#   - /SAN/.../VCF_noProbSitesnoMutPos_tsv/ (iVar TSV files)
#   - /SAN/.../variant_blacklist.tsv (artifact blacklist)
# Output:
#   - SuppFig7_shared_variants_filtered.pdf
#   - shared_variants_summary.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(jsonlite)
  library(readr)
})

# -------- USER PATHS — adjust to your cluster ---------
household_pairs_json <- "/SAN/breuerlab/LUNARstudy/VirusWatch_Arturo/bottleneck_out/pairs.json"
random_pairs_json    <- "/SAN/breuerlab/LUNARstudy/VirusWatch_Arturo/bottleneck_out/random_pairs.json"
tsv_dir              <- "/SAN/breuerlab/LUNARstudy/VirusWatch_Arturo/VW_VCF/VCF_noProbSitesnoMutPos_tsv/"
blacklist_path       <- "/SAN/breuerlab/LUNARstudy/VirusWatch_Arturo/variant_blacklist.tsv"
output_dir           <- "/SAN/breuerlab/LUNARstudy/VirusWatch_Arturo/figures_revised"
# ------------------------------------------------------

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Frequency window (Sobel Leonard validated)
AF_MIN <- 0.03
AF_MAX <- 0.50

# -------- Load blacklist (CHROM:POS:REF>ALT, minor allele convention) ---------
load_blacklist <- function(path) {
  if (!file.exists(path)) return(character(0))
  bl <- read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  bl$site_key
}
blacklist <- load_blacklist(blacklist_path)
cat("Loaded blacklist:", length(blacklist), "positions\n")

# -------- Read iVar TSV, return tibble of MINOR alleles in [AF_MIN, AF_MAX] ---
read_minor_variants <- function(tsv_path, af_min = AF_MIN, af_max = AF_MAX,
                                snvs_only = TRUE) {
  df <- tryCatch(
    read_tsv(tsv_path, show_col_types = FALSE),
    error = function(e) return(NULL)
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)
  names(df) <- trimws(names(df))
  if (!all(c("CHROM","POS","REF","ALT","ALT_FREQ","FILTER") %in% names(df))) return(NULL)

  df <- df %>%
    mutate(across(where(is.character), trimws)) %>%
    filter(FILTER == "PASS") %>%
    mutate(ALT_FREQ = as.numeric(ALT_FREQ),
           POS      = as.integer(POS)) %>%
    filter(!is.na(ALT_FREQ), !is.na(POS))

  if (snvs_only) {
    df <- df %>% filter(nchar(REF) == 1, nchar(ALT) == 1,
                        !startsWith(ALT, "+"), !startsWith(ALT, "-"))
  }
  if (nrow(df) == 0) return(NULL)

  # Compute minor freq and identify minor allele direction
  df <- df %>%
    mutate(minor_freq   = pmin(ALT_FREQ, 1 - ALT_FREQ),
           minor_is_alt = ALT_FREQ <= 0.5,
           # site_key in same format as the blacklist
           site_key = ifelse(minor_is_alt,
                             paste0(CHROM, ":", POS, ":", REF, ">", ALT),
                             paste0(CHROM, ":", POS, ":", ALT, ">", REF))) %>%
    filter(minor_freq >= af_min, minor_freq <= af_max)

  df
}

# -------- Count shared minority variants between two samples ----------------
count_shared <- function(sid1, sid2, samples, apply_filter = TRUE) {
  v1 <- samples[[sid1]]
  v2 <- samples[[sid2]]
  if (is.null(v1) || is.null(v2)) return(NA_integer_)
  if (nrow(v1) == 0 || nrow(v2) == 0) return(0L)
  k1 <- v1$site_key
  k2 <- v2$site_key
  if (apply_filter) {
    k1 <- setdiff(k1, blacklist)
    k2 <- setdiff(k2, blacklist)
  }
  length(intersect(k1, k2))
}

# -------- Load all sample data once -----------------------------------------
tsv_files <- list.files(tsv_dir, pattern = "\\.tsv$", full.names = TRUE)
cat("Loading", length(tsv_files), "TSV files...\n")

samples <- list()
for (i in seq_along(tsv_files)) {
  sid <- sub("\\.tsv$", "", basename(tsv_files[i]))
  samples[[sid]] <- read_minor_variants(tsv_files[i])
  if (i %% 50 == 0) cat("  loaded", i, "/", length(tsv_files), "\n")
}
cat("Loaded SNP data for", length(samples), "samples\n")

# -------- Load pair lists ----------------------------------------------------
pairs_from_json <- function(path) {
  raw <- fromJSON(path, simplifyVector = FALSE)
  do.call(rbind, lapply(raw, function(p) {
    data.frame(
      sample1 = sub("\\.tsv$", "", basename(p[[1]])),
      sample2 = sub("\\.tsv$", "", basename(p[[2]])),
      stringsAsFactors = FALSE
    )
  }))
}

hh_pairs   <- pairs_from_json(household_pairs_json)
rand_pairs <- pairs_from_json(random_pairs_json)
cat("Household pairs:", nrow(hh_pairs), " | Random pairs:", nrow(rand_pairs), "\n")

# -------- Count shared variants per pair, BEFORE and AFTER filter -----------
compute_pair_counts <- function(pairs_df, label) {
  result <- pairs_df %>%
    rowwise() %>%
    mutate(
      shared_before = count_shared(sample1, sample2, samples, apply_filter = FALSE),
      shared_after  = count_shared(sample1, sample2, samples, apply_filter = TRUE)
    ) %>%
    ungroup() %>%
    mutate(group = label)
  result
}

cat("Counting shared variants for household pairs...\n")
hh_counts <- compute_pair_counts(hh_pairs, "Household pairs")

cat("Counting shared variants for random pairs...\n")
rand_counts <- compute_pair_counts(rand_pairs, "Random pairs")

all_counts <- bind_rows(hh_counts, rand_counts) %>%
  mutate(group = factor(group, levels = c("Random pairs", "Household pairs")))

# -------- Stats --------------------------------------------------------------
summary_table <- all_counts %>%
  group_by(group) %>%
  summarise(n_pairs       = n(),
            median_before = median(shared_before, na.rm = TRUE),
            mean_before   = mean(shared_before, na.rm = TRUE),
            median_after  = median(shared_after, na.rm = TRUE),
            mean_after    = mean(shared_after, na.rm = TRUE),
            .groups = "drop")
print(summary_table)
write.csv(summary_table,
          file.path(output_dir, "shared_variants_summary.csv"),
          row.names = FALSE)

# Wilcoxon tests
wilcox_before <- wilcox.test(
  shared_before ~ group,
  data = all_counts,
  alternative = "less"  # random < household
)
wilcox_after <- wilcox.test(
  shared_after ~ group,
  data = all_counts,
  alternative = "less"
)
cat("\nWilcoxon (random < household):\n")
cat("  BEFORE filter: p =", format.pval(wilcox_before$p.value), "\n")
cat("  AFTER filter:  p =", format.pval(wilcox_after$p.value), "\n")

# -------- Pivot to long form for facetted plot ------------------------------
all_long <- all_counts %>%
  select(sample1, sample2, group, shared_before, shared_after) %>%
  tidyr::pivot_longer(c(shared_before, shared_after),
                      names_to = "filter", values_to = "shared") %>%
  mutate(filter = factor(filter,
                         levels = c("shared_before", "shared_after"),
                         labels = c("Before artifact filter",
                                    "After artifact filter")))

# -------- Plot --------------------------------------------------------------
p <- ggboxplot(all_long, x = "group", y = "shared",
               fill = "group", alpha = 0.7,
               outlier.shape = 21, outlier.size = 1) +
  scale_fill_manual(values = c("Random pairs" = "#FFB6C1",
                               "Household pairs" = "#1A85FF")) +
  facet_wrap(~ filter, scales = "free_y") +
  stat_compare_means(comparisons = list(c("Random pairs", "Household pairs")),
                     label = "p.signif",
                     method = "wilcox.test",
                     method.args = list(alternative = "less")) +
  labs(x = NULL,
       y = "Shared minority variants per pair\n(3-50% MAF window)",
       title = NULL) +
  theme_pubr() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 14, face = "bold"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", face = "bold", size = 14))

print(p)
ggsave(file.path(output_dir, "SuppFig7_shared_variants_filtered.pdf"),
       plot = p, width = 10, height = 6)
ggsave(file.path(output_dir, "SuppFig7_shared_variants_filtered.png"),
       plot = p, width = 10, height = 6, dpi = 300)

cat("\nDone. Plot saved to:", file.path(output_dir, "SuppFig7_shared_variants_filtered.pdf"), "\n")
