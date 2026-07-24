# ==============================================================================
# Supplementary Figure X — Characterisation of recurrent artifact variant positions
#
# Three-panel figure documenting that the recurrent low-frequency variants
# excluded from the bottleneck analysis are spatially associated with
# ARTIC v4.1 primer binding sites and behave as primer-binding-site
# artifacts rather than classical cross-sample read misassignment.
#
# Panel a: Position of artifact variants overlaid on ARTIC v4.1 primer scheme
#          (genome track)
# Panel b: Table of 13 artifact positions with prevalence, median freq,
#          and primer-adjacency annotation
# Panel c: Lineage-stratified frequency distribution at position 14960
#          (the most prevalent artifact), showing variants distributed
#          across all major lineages at low frequency without ever
#          being fixed
#
# All inputs are small TSV/CSV files you copy from the cluster:
#   - variant_blacklist.tsv      from 06_build_variant_blacklist.py
#   - artifact_position_lineage_data.tsv   (generated below — needs cluster step)
# Plus a small inline data frame of primer-adjacency annotation we determined
# from your cluster run.
# ==============================================================================

# ----------- USER PATHS ----------------------------------------------------
input_dir  <- "~/UCL Dropbox/Laura Buggiotti/Mac/Desktop/VirusWatch/bottleneck_figures"
output_dir <- "~/UCL Dropbox/Laura Buggiotti/Mac/Desktop/VirusWatch/bottleneck_figures/output"

blacklist_path        <- file.path(input_dir, "variant_blacklist.tsv")
position_14960_path   <- file.path(input_dir, "position_14960_per_sample.tsv")
# (this is the file 08_check_index_hopping.py wrote out)
# ---------------------------------------------------------------------------


dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Install missing packages
required <- c("dplyr","tidyr","ggplot2","ggpubr","cowplot","forcats","scales",
              "gridExtra","grid","gtable")
missing <- required[!required %in% installed.packages()[,"Package"]]
if (length(missing) > 0) install.packages(missing, repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggpubr)
  library(cowplot); library(forcats); library(scales)
  library(gridExtra); library(grid); library(gtable)
})

# ============================================================================
# Inline data: artifact positions + primer adjacency (from your cluster run)
# ============================================================================

# 13 positions identified by the cohort-wide artifact screen (prevalence >5%)
# Columns derived from variant_blacklist.tsv and your primer-adjacency check
artifact_table <- tibble::tribble(
  ~position, ~ref, ~alt, ~n_samples, ~pct_samples, ~median_freq, ~sd_freq, ~near_primer, ~primer_name,
  14960,  "A", "T",  234,  97.9, 0.0948, 0.0629, TRUE,  "nCoV_50_LEFT",
  5744,  "T", "C",  200,  83.7, 0.0451, 0.0099, FALSE, NA_character_,
  15521,  "T", "A",  142,  59.4, 0.0632, 0.1052, TRUE,  "nCoV_52_LEFT",
  11291,  "G", "A",   49,  20.5, 0.0541, 0.0314, TRUE,  "nCoV_38_LEFT",
  8835,  "T", "C",   33,  13.8, 0.0632, 0.0791, FALSE, NA_character_,
  11956,  "T", "C",   29,  12.1, 0.1227, 0.0963, TRUE,  "nCoV_40_LEFT",
  28363,  "A", "T",   27,  11.3, 0.0833, 0.0319, FALSE, NA_character_,
  28370,  "A", "T",   27,  11.3, 0.1333, 0.0772, TRUE,  "nCoV_94_RIGHT",
  21639,  "C", "A",   24,  10.0, 0.0901, 0.0356, FALSE, NA_character_,
  18012,  "T", "C",   23,   9.6, 0.1194, 0.1267, TRUE,  "nCoV_59_RIGHT",
  1545,  "T", "C",   18,   7.5, 0.2126, 0.1566, TRUE,  "nCoV_6_LEFT",
  18014,  "C", "T",   17,   7.1, 0.0455, 0.1046, TRUE,  "nCoV_59_RIGHT",
  23118,  "A", "T",   12,   5.0, 0.0392, 0.0072, TRUE,  "nCoV_76_RIGHT_alt1"
) %>% mutate(
  variant_label = sprintf("%d %s>%s", position, ref, alt),
  variant_label = factor(variant_label, levels = variant_label[order(position)])
)

# ============================================================================
# Panel A: Genome track of artifact positions
# ============================================================================

# A schematic representation of the SARS-CoV-2 genome (29903 bp)
# with major gene boundaries marked, artifact positions plotted

genes <- tibble::tribble(
  ~gene,   ~start,  ~end,
  "ORF1ab",   266, 21555,
  "S",     21563, 25384,
  "ORF3a", 25393, 26220,
  "E",     26245, 26472,
  "M",     26523, 27191,
  "ORF6",  27202, 27387,
  "ORF7a", 27394, 27759,
  "ORF7b", 27756, 27887,
  "ORF8",  27894, 28259,
  "N",     28274, 29533
) %>% mutate(
  ymid = 0.5,
  fill = rep(c("#E8E8E8","#C8C8C8"), length.out = n())
)

# Only show major gene labels to avoid overcrowding in the 3' end
genes_to_label <- genes %>% filter(gene %in% c("ORF1ab", "S", "N"))
margin <- ggplot2::margin

panel_a <- ggplot() +
  # Genome scale bar
  annotate("segment", x = 0, xend = 29903, y = 0.5, yend = 0.5,
           colour = "grey70", linewidth = 0.6) +
  # Gene blocks (still draw all of them, just don't label the small ones)
  geom_rect(data = genes,
            aes(xmin = start, xmax = end, ymin = 0.35, ymax = 0.65,
                fill = fill), colour = "grey50", linewidth = 0.2) +
  geom_text(data = genes_to_label,
            aes(x = (start+end)/2, y = 0.5, label = gene),
            size = 3, vjust = 0.5, fontface = "bold") +
  scale_fill_identity() +
  # Artifact variants — colour-coded by primer adjacency
  geom_segment(data = artifact_table,
               aes(x = position, xend = position, y = 0.7, yend = 1.0,
                   colour = near_primer), linewidth = 0.8) +
  geom_point(data = artifact_table,
             aes(x = position, y = 1.05, colour = near_primer,
                 size = pct_samples), shape = 19) +
  scale_colour_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2"),
                      labels = c("TRUE" = "Within 30 bp of ARTIC primer",
                                 "FALSE" = "Not near primer"),
                      name = NULL) +
  scale_size_continuous(name = "% of samples",
                        range = c(2, 7),
                        breaks = c(5, 25, 50, 100)) +
  # Label the top artifacts (more headroom for the labels)
  geom_text(data = artifact_table %>% filter(pct_samples > 20),
            aes(x = position, y = 1.20, label = variant_label),
            angle = 90, hjust = 0, size = 2.8) +
  scale_x_continuous(
    name = "Position on SARS-CoV-2 genome (MN908947.3)",
    breaks = seq(0, 30000, 5000),
    labels = comma,
    limits = c(-500, 30500),
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(0.3, 2.0), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 11, face = "bold"),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.box = "horizontal",
    legend.text = element_text(size = 9),
    plot.margin = ggplot2::margin(t = 10, r = 15, b = 25, l = 15)
  )

print(panel_a)

# ============================================================================
# Panel B: Artifact table
# ============================================================================

# Render as a table grob
table_for_display <- artifact_table %>%
  mutate(
    Position = position,
    Variant = sprintf("%s>%s", ref, alt),
    `n (%) samples` = sprintf("%d (%.1f%%)", n_samples, pct_samples),
    `Median freq` = sprintf("%.3f", median_freq),
    `SD freq` = sprintf("%.3f", sd_freq),
    `ARTIC primer (within 30 bp)` = ifelse(near_primer, primer_name, "—")
  ) %>%
  select(Position, Variant, `n (%) samples`, `Median freq`, `SD freq`,
         `ARTIC primer (within 30 bp)`)

tab_theme <- ttheme_minimal(
  core    = list(fg_params = list(cex = 0.8),
                 bg_params = list(fill = c("white","grey95"), col = NA),
                 padding   = unit(c(4, 4), "mm")),
  colhead = list(fg_params = list(cex = 0.85, fontface = "bold"),
                 bg_params = list(fill = "grey85", col = NA),
                 padding   = unit(c(4, 4), "mm"))
)

panel_b <- tableGrob(table_for_display, rows = NULL, theme = tab_theme)

# Header annotation — split across two lines for readability
header_text <- textGrob(
  paste0("13 recurrent variant positions identified by cohort-wide screen ",
         "(>5% of 239 samples).\n",
         "9/13 (69%) within 30 bp of ARTIC v4.1 primer binding site",
         " \u2014 a ~5\u00d7 enrichment over expectation."),
  gp = gpar(fontsize = 10, fontface = "italic"),
  just = "centre"
)

# Tight stacking: header sits directly above the table with minimal gap
panel_b_full <- arrangeGrob(
  header_text,
  panel_b,
  ncol = 1,
  heights = unit.c(unit(2.2, "lines"), unit(1, "null")),
  padding = unit(0.1, "cm")
)

# ============================================================================
# Panel C: Lineage-stratified frequency distribution at position 14960
# ============================================================================

# Load the per-sample frequency data for position 14960
pos_data <- read.delim(position_14960_path, sep = "\t",
                       stringsAsFactors = FALSE)
# Columns: sample, freq, ref, alt, lineage

# Group rare lineages into "Other" for visual clarity
top_lineages <- pos_data %>%
  filter(!is.na(lineage), lineage != "NA") %>%
  count(lineage, sort = TRUE) %>%
  slice_head(n = 6) %>%
  pull(lineage)

pos_data <- pos_data %>%
  mutate(
    lineage_grouped = case_when(
      is.na(lineage) | lineage == "NA" ~ "Unassigned",
      lineage %in% top_lineages ~ lineage,
      TRUE ~ "Other"
    ),
    lineage_grouped = factor(lineage_grouped,
                             levels = c(top_lineages, "Other", "Unassigned"))
  )

panel_c <- ggplot(pos_data, aes(x = freq, y = lineage_grouped,
                                colour = lineage_grouped)) +
  # Shaded artifact band
  annotate("rect", xmin = 0.03, xmax = 0.15, ymin = -Inf, ymax = Inf,
           fill = "#FFC4C4", alpha = 0.3) +
  # Shaded fixation band
  annotate("rect", xmin = 0.90, xmax = 1.00, ymin = -Inf, ymax = Inf,
           fill = "#C4E0C4", alpha = 0.3) +
  geom_jitter(width = 0, height = 0.25, alpha = 0.7, size = 1.5) +
  scale_colour_brewer(palette = "Set2", guide = "none") +
  scale_x_continuous(
    name = sprintf("ALT_FREQ at position 14960 (A>T)"),
    breaks = c(0, 0.03, 0.15, 0.5, 0.9, 1.0),
    labels = c("0", "0.03", "0.15", "0.5", "0.9", "1.0"),
    limits = c(-0.02, 1.02),
    expand = c(0, 0)
  ) +
  labs(y = "Pangolin lineage") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    panel.grid.major.y = element_line(colour = "grey90"),
    panel.grid.minor = element_blank(),
    plot.margin = ggplot2::margin(t = 25, r = 15, b = 5, l = 15)
  ) +
  # Place labels OUTSIDE the panel using coord_cartesian + clip = "off"
  coord_cartesian(clip = "off") +
  annotate("text", x = 0.09, y = length(levels(pos_data$lineage_grouped)) + 0.7,
           label = "Artifact band\n(3-15%)",
           size = 3.2, colour = "#9C2424", hjust = 0.5, fontface = "italic") +
  annotate("text", x = 0.95, y = length(levels(pos_data$lineage_grouped)) + 0.7,
           label = "Fixation\n(\u226590%)",
           size = 3.2, colour = "#244C24", hjust = 0.5, fontface = "italic")

print(panel_c)

# ============================================================================
# Assemble multi-panel figure
# ============================================================================

# Panel B is a grob; wrap it in plot_grid-compatible form
top_row <- plot_grid(panel_a, labels = "a", label_size = 18)
mid_row <- plot_grid(panel_b_full, labels = "b", label_size = 18)
bot_row <- plot_grid(panel_c, labels = "c", label_size = 18)

final_fig <- plot_grid(
  top_row, mid_row, bot_row,
  ncol = 1,
  rel_heights = c(0.65, 0.85, 0.95),
  align = "v",
  axis = "lr"
)

ggsave(file.path(output_dir, "SuppFigX_artifact_characterisation.pdf"),
       plot = final_fig, width = 11, height = 13)
ggsave(file.path(output_dir, "SuppFigX_artifact_characterisation.png"),
       plot = final_fig, width = 11, height = 13, dpi = 300)

cat("Figure saved to:", file.path(output_dir, "SuppFigX_artifact_characterisation.pdf"), "\n")
