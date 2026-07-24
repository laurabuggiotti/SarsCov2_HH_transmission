# ==============================================================================
# Figure 4 — two-panel main figure
#
# Panel A: filtered per-pair Nb (artifact blacklist applied, 3-50% window)
# Panel B: TransPhylo posterior + filtered shared-variant proportions
#
# Both panels use hh2_X labels to match the original figure naming.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(cowplot)
  library(readr)
})

# -------- USER PATHS --------
pair_results_path <- "./pair_results_3to50_FILTERED.tsv"
posterior_path    <- "./PP_propShared_value_FILTERED.txt"
hh_mapping_path   <- "./hh_id_mapping.csv"
output_dir        <- "./figures_revised"
# ----------------------------

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

MIN_SITES <- 2   # min informative sites for Nb estimate

# ==============================================================================
# Load and merge data
# ==============================================================================

# Per-pair Nb results
bb <- read.delim(pair_results_path, sep = "\t", stringsAsFactors = FALSE)

# hh ID mapping: hh2_id, sample_a, sample_b
mapping <- read_csv(hh_mapping_path, show_col_types = FALSE)

# Function: extract donor and recipient sample IDs from pair label
parse_pair_label <- function(lbl) {
  parts <- strsplit(lbl, "__to__", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    return(c(parts[1], parts[2]))
  }
  return(c(NA, NA))
}

# Annotate per-pair results with hh2_id
bb <- bb %>%
  filter(pair != "bn_pooled_filt_result") %>%
  rowwise() %>%
  mutate(
    donor    = parse_pair_label(pair)[1],
    recipient = parse_pair_label(pair)[2]
  ) %>%
  ungroup()

# Lookup hh2_id by matching either sample
hh_by_sample <- mapping %>%
  pivot_longer(c(sample_a, sample_b), values_to = "sample") %>%
  select(hh2_id, sample)

bb <- bb %>%
  left_join(hh_by_sample, by = c("donor" = "sample")) %>%
  rename(hh2_id_donor = hh2_id) %>%
  left_join(hh_by_sample, by = c("recipient" = "sample")) %>%
  rename(hh2_id_recip = hh2_id) %>%
  mutate(hh2_id = coalesce(hh2_id_donor, hh2_id_recip)) %>%
  select(-hh2_id_donor, -hh2_id_recip)

cat("Pairs after annotation:", nrow(bb), "\n")

# Keep only pairs with valid Nb and enough sites
bb <- bb %>%
  mutate(
    Nb_num     = suppressWarnings(as.numeric(Nb)),
    CIl_num    = suppressWarnings(as.numeric(CI_low)),
    CIh_num    = suppressWarnings(as.numeric(CI_high)),
    n_sites    = as.integer(n_sites)
  ) %>%
  filter(!is.na(Nb_num), !is.na(hh2_id), n_sites >= MIN_SITES)

cat("Pairs with valid Nb (>=", MIN_SITES, "sites):", nrow(bb), "\n")
cat("Median Nb:", median(bb$Nb_num), "\n")
cat("IQR:", paste(quantile(bb$Nb_num, c(0.25, 0.75)), collapse = " - "), "\n")
cat("Range:", min(bb$Nb_num), "-", max(bb$Nb_num), "\n")
cat("Pairs at Nb=1:", sum(bb$Nb_num == 1), "/", nrow(bb), "\n")

# ==============================================================================
# Panel A — forest plot of per-pair Nb
# ==============================================================================

bb <- bb %>%
  mutate(
    hh2_id  = fct_reorder(hh2_id, Nb_num),
    Nb_log  = log10(Nb_num),
    CIl_log = log10(CIl_num),
    CIh_log = log10(CIh_num)
  )

panel_a <- ggplot(bb, aes(x = hh2_id, y = Nb_log, ymin = CIl_log, ymax = CIh_log)) +
  geom_linerange(colour = "olivedrab", linewidth = 0.5) +
  geom_point(colour = "olivedrab", size = 2.5, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_y_continuous(
    name = expression(bold("Transmission bottleneck size " * (N[b]))),
    breaks = c(0, 1, 2, 3),
    labels = c("1", "10", "100", "1000")
  ) +
  labs(x = NULL) +
  theme_bw() +
  theme(
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.y  = element_text(size = 12),
    axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank()
  )

print(panel_a)

# ==============================================================================
# Panel B — TransPhylo posterior + filtered shared-variant proportions
# ==============================================================================

# Load the filtered posterior file
df <- read.delim(posterior_path, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
df <- df %>% filter(!is.na(posteriorProb), posteriorProb != "")
df$posteriorProb <- as.numeric(df$posteriorProb)

# Match factor levels to Panel A ordering
panel_a_levels <- levels(bb$hh2_id)
df <- df %>% mutate(hh2_id = factor(hh2_type, levels = panel_a_levels))

# Pivot prop_shared_a/b to long for plotting
df_shared <- df %>%
  select(hh2_id, prop_shared_a, prop_shared_b) %>%
  mutate(
    prop_shared_a = suppressWarnings(as.numeric(prop_shared_a)),
    prop_shared_b = suppressWarnings(as.numeric(prop_shared_b))
  ) %>%
  pivot_longer(c(prop_shared_a, prop_shared_b),
               names_to = "type", values_to = "value") %>%
  filter(!is.na(hh2_id), !is.na(value))

panel_b <- ggplot() +
  geom_line(data = df_shared,
            aes(x = hh2_id, y = value, group = hh2_id),
            colour = "steelblue", linewidth = 0.5) +
  geom_point(data = df_shared,
             aes(x = hh2_id, y = value),
             colour = "steelblue", size = 2.5) +
  geom_point(data = df %>% filter(!is.na(hh2_id)),
             aes(x = hh2_id, y = posteriorProb),
             colour = "darkorange", size = 3, shape = 17) +
  ylim(0, 1) +
  labs(y = "Probability / Proportion", x = NULL) +
  theme_bw() +
  theme(
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.y  = element_text(size = 12),
    axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
    axis.ticks.x = element_blank()
  )

print(panel_b)

# ==============================================================================
# Combine
# ==============================================================================

final_plot <- plot_grid(
  panel_a, panel_b,
  ncol = 1,
  rel_heights = c(1, 0.55),
  align = "v",
  axis = "rlbt",
  labels = c("a", "b"),
  label_size = 18
)

ggsave(file.path(output_dir, "Figure_4.pdf"),
       plot = final_plot, width = 14, height = 10)
ggsave(file.path(output_dir, "Figure4_revised.png"),
       plot = final_plot, width = 14, height = 10, dpi = 300)

