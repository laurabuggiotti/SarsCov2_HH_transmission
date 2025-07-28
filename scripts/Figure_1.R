# ==============================================================================
# SARS-CoV-2 Lineage Distribution Analysis
# Figure 1: Participants by viral lineage and variant of concern
# 
# This script generates a horizontal bar chart showing the distribution of 
# SARS-CoV-2 lineages among study participants, grouped and colored by 
# variant of concern (VOC/VUI).
#
# Input: lineage_VL.csv - contains participant lineage and VOC data
# Output: lineage_participants.pdf - horizontal bar chart
# ==============================================================================

# Load required libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(forcats)

# Set file paths (adjust as needed for your directory structure)
data_path <- "data/lineage_VL.csv"
output_dir <- "figures"
output_file <- "lineage_participants.pdf"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read and explore data
cat("Reading data from:", data_path, "\n")
vw_metadata <- read.csv(data_path)

# Data exploration (optional - comment out for production)
cat("Unique samples:", length(unique(vw_metadata$SAMPLE)), "\n")
cat("Data structure:\n")
str(vw_metadata)
cat("Unique VOC/VUI categories:", paste(unique(vw_metadata$VOC_VUI), collapse = ", "), "\n")

# Process data: filter out missing VOC data and count participants by lineage
No_lineage <- vw_metadata %>%
  filter(VOC_VUI != 'NA') %>% 
  group_by(VOC_VUI, LINEAGE) %>%
  tally() %>%
  ungroup()

# Set factor levels for consistent ordering
No_lineage$VOC_VUI <- factor(No_lineage$VOC_VUI, 
                            levels = c('Omicron', 'Omicron (BA.2-like)', 'Omicron (BA.5-like)'))

No_lineage$LINEAGE <- factor(No_lineage$LINEAGE, 
                            levels = unique(No_lineage$LINEAGE[order(No_lineage$VOC_VUI)]))

# Create the plot
plot <- No_lineage %>%    
  ggplot(aes(x = n, y = fct_inorder(LINEAGE), fill = VOC_VUI)) +   
  labs(title = "", 
       y = "Lineage", 
       x = "# participants",
       fill = "VOC/VUI") +   
  geom_bar(stat = "identity") +    
  scale_fill_manual(values = c('#1A85FF', '#5D3A9B', '#E66100')) +   
  theme_pubr() +   
  theme(     
    legend.position = "right",     
    axis.title.x = element_text(size = 14, face = 'bold'),     
    axis.title.y = element_text(size = 14, face = 'bold'),     
    axis.text.x.bottom = element_text(vjust = 0.5, size = 12),     
    axis.text.y = element_text(size = 12),     
    strip.text.x = element_text(face = "bold", colour = "white", size = 12),     
    strip.text.y = element_text(face = "bold", colour = "white", size = 12),     
    strip.background = element_rect(colour = "black", fill = "black")   
  ) 

# Display plot
print(plot)

# Save plot
output_path <- file.path(output_dir, output_file)
ggsave(plot = plot, 
       height = 5.5, 
       width = 7, 
       filename = output_path,
       device = "pdf")

cat("Plot saved to:", output_path, "\n")

# Print summary statistics
cat("\nSummary:\n")
cat("Total participants included:", sum(No_lineage$n), "\n")
cat("Number of unique lineages:", length(unique(No_lineage$LINEAGE)), "\n")
cat("Participants by VOC/VUI:\n")
print(No_lineage %>% group_by(VOC_VUI) %>% summarise(total = sum(n)))

# ==============================================================================
# End of script
# ==============================================================================
