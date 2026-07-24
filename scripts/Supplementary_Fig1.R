# ==============================================================================
# SARS-CoV-2 Geographic Distribution Analysis - Supplementary Material
# Supplementary Figure 2: Geographic distribution of study participants across UK
# 
# This script creates a map visualization showing the geographic distribution
# of study participants across UK regions, with points colored by region and
# shaped by household type. This provides important context for understanding
# the spatial representation of the study cohort.
#
# Input: 
#   - UK administrative boundaries shapefile (Local Authority Districts)
#   - LocalAuthorithy_lineage.csv - Study participant geographic data
# Output: 
#   - UK map with participant locations marked by region and household type
# ==============================================================================

# Load required libraries
library(osmdata)      # for OpenStreetMap data access
library(sf)           # for spatial data handling
library(mapview)      # for interactive mapping (optional)
library(rgdal)        # for spatial data reading
library(ggplot2)      # for plotting
library(dplyr)        # for data manipulation
library(giscoR)       # for country boundaries
library(maps)         # for world cities data
library(HYPEtools)    # for distinct color palettes
library(pals)         # for polychrome color palettes

# Set file paths (adjust as needed for your directory structure)
shapefile_path <- "data/LAD_MAY_2024_UK_BFC.shp"
participant_data_path <- "data/LocalAuthorithy_lineage.csv"
output_dir <- "figures"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==============================================================================
# PART 1: Load and process UK administrative boundaries
# ==============================================================================

cat("Loading UK administrative boundaries...\n")

# Method 1: Using giscoR for country boundaries (recommended)
cat("Downloading UK country boundaries from giscoR...\n")
UK <- gisco_get_countries(country = "UK", resolution = 1)

# Method 2: Load local shapefile for detailed boundaries (if available)
if (file.exists(shapefile_path)) {
  cat("Loading detailed administrative boundaries from:", shapefile_path, "\n")
  region <- read_sf(shapefile_path)
  cat("Loaded", nrow(region), "administrative districts\n")
  
  # Extract district coordinates
  district_vw <- region %>%
    select(LAD24CD, LONG, LAT) %>%
    distinct()
  
  cat("Extracted coordinates for", nrow(district_vw), "unique districts\n")
} else {
  cat("Warning: Detailed shapefile not found at:", shapefile_path, "\n")
  cat("Will use basic UK boundaries only\n")
  district_vw <- NULL
}

# ==============================================================================
# PART 2: Load and process participant geographic data
# ==============================================================================

cat("Loading participant geographic data from:", participant_data_path, "\n")
vw_district_lineage <- read.csv(participant_data_path)

cat("Participant data contains", nrow(vw_district_lineage), "records\n")
cat("Data structure:\n")
str(vw_district_lineage)

# Check for required columns
required_cols <- c("laua..local.authority.")
missing_cols <- setdiff(required_cols, colnames(vw_district_lineage))

if (length(missing_cols) > 0) {
  cat("Warning: Missing expected columns:", paste(missing_cols, collapse = ", "), "\n")
  cat("Available columns:", paste(colnames(vw_district_lineage), collapse = ", "), "\n")
}

# ==============================================================================
# PART 3: Merge participant data with geographic coordinates
# ==============================================================================

if (!is.null(district_vw)) {
  cat("Merging participant data with district coordinates...\n")
  
  metadata_lineage <- merge(vw_district_lineage, district_vw, 
                           by.x = 'laua..local.authority.', 
                           by.y = 'LAD24CD', 
                           all = FALSE)
  
  cat("Successfully merged", nrow(metadata_lineage), "records with coordinates\n")
  cat("Merged data structure:\n")
  str(metadata_lineage)
  
  # Use merged data for plotting
  metadata <- metadata_lineage
} else {
  cat("Using participant data without detailed coordinates\n")
  # You may need to add default coordinates or use alternative data source
  metadata <- vw_district_lineage
}

# ==============================================================================
# PART 4: Data exploration and validation
# ==============================================================================

cat("\n=== DATA EXPLORATION ===\n")

# Check for required plotting variables
required_plot_vars <- c("Region", "HH_type", "LONG", "LAT")
available_vars <- intersect(required_plot_vars, colnames(metadata))
missing_vars <- setdiff(required_plot_vars, colnames(metadata))

cat("Available plotting variables:", paste(available_vars, collapse = ", "), "\n")
if (length(missing_vars) > 0) {
  cat("Missing plotting variables:", paste(missing_vars, collapse = ", "), "\n")
}

# Summary of regions and household types
if ("Region" %in% colnames(metadata)) {
  cat("\nUnique regions (", length(unique(metadata$Region)), "):\n")
  region_counts <- table(metadata$Region)
  print(region_counts)
}

if ("HH_type" %in% colnames(metadata)) {
  cat("\nHousehold types:\n")
  hh_counts <- table(metadata$HH_type)
  print(hh_counts)
}

# Geographic range
if (all(c("LONG", "LAT") %in% colnames(metadata))) {
  cat("\nGeographic range:\n")
  cat("Longitude:", round(min(metadata$LONG, na.rm = TRUE), 3), "to", 
      round(max(metadata$LONG, na.rm = TRUE), 3), "\n")
  cat("Latitude:", round(min(metadata$LAT, na.rm = TRUE), 3), "to", 
      round(max(metadata$LAT, na.rm = TRUE), 3), "\n")
}

# ==============================================================================
# PART 5: Create color palette for regions
# ==============================================================================

if ("Region" %in% colnames(metadata)) {
  n_regions <- length(unique(metadata$Region))
  cat("Creating color palette for", n_regions, "regions...\n")
  
  # Generate distinct colors
  if (n_regions <= 33) {
    # Use polychrome palette for up to 33 colors
    region_colors <- unname(polychrome(n_regions))
  } else {
    # For more regions, use distinctColorPalette
    if (require(HYPEtools, quietly = TRUE)) {
      region_colors <- distinctColorPalette(count = n_regions, seed = NULL, darken = 0)
    } else {
      # Fallback to rainbow colors
      region_colors <- rainbow(n_regions)
    }
  }
} else {
  cat("No region variable found - using default colors\n")
  region_colors <- c("#1f77b4")  # Default blue
}

# ==============================================================================
# PART 6: Create the map visualization
# ==============================================================================

cat("Creating UK geographic distribution map...\n")

# Check if we have all required variables for the full plot
can_create_full_plot <- all(c("LONG", "LAT", "Region", "HH_type") %in% colnames(metadata))

if (can_create_full_plot) {
  uk_plot_region <- ggplot() +
    geom_sf(data = UK, fill = 'gray', alpha = 0.1) +
    geom_point(data = metadata, 
               aes(x = LONG, y = LAT, shape = HH_type, color = Region), 
               size = 2.5) +
    scale_size_continuous(range = c(1, 12)) +
    scale_color_manual(values = region_colors, name = "Region") +
    scale_shape_manual(values = c(16, 17, 15), name = "Household\nType") +  # Adjust shapes as needed
    theme_void() +
    ylim(50, 59) +  # Focus on UK mainland
    labs(title = "Geographic Distribution of Study Participants",
         subtitle = "Points colored by region, shaped by household type") +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10, face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
} else {
  # Create simplified plot with available data
  cat("Creating simplified map with available data...\n")
  uk_plot_region <- ggplot() +
    geom_sf(data = UK, fill = 'gray', alpha = 0.3) +
    theme_void() +
    ylim(50, 59) +
    labs(title = "UK Study Area",
         subtitle = "Geographic context for SARS-CoV-2 transmission study") +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  
  cat("Note: Limited geographic data available for detailed participant mapping\n")
}

# Display the plot
print(uk_plot_region)

# ==============================================================================
# PART 7: Save the map
# ==============================================================================

output_path <- file.path(output_dir, "SupplementaryFigure2_geographic_distribution.pdf")
ggsave(plot = uk_plot_region, height = 14, width = 12, filename = output_path)

cat("Supplementary Figure 2 saved to:", output_path, "\n")

# Also save as PNG for presentations
png_path <- file.path(output_dir, "SupplementaryFigure2_geographic_distribution.png")
ggsave(plot = uk_plot_region, height = 14, width = 12, filename = png_path, dpi = 300)

cat("PNG version saved to:", png_path, "\n")

# ==============================================================================
# PART 8: Summary statistics and interpretation
# ==============================================================================

cat("\n=== SUPPLEMENTARY FIGURE 2 ANALYSIS SUMMARY ===\n")

if (can_create_full_plot) {
  # Geographic coverage analysis
  cat("Geographic coverage analysis:\n")
  cat("- Total participants plotted:", nrow(metadata), "\n")
  cat("- Number of regions represented:", length(unique(metadata$Region)), "\n")
  cat("- Number of household types:", length(unique(metadata$HH_type)), "\n")
  
  # Regional distribution
  if ("Region" %in% colnames(metadata)) {
    regional_summary <- metadata %>%
      group_by(Region) %>%
      summarise(n_participants = n(), .groups = 'drop') %>%
      arrange(desc(n_participants))
    
    cat("\nTop 5 regions by participant count:\n")
    print(head(regional_summary, 5))
  }
  
  # Household type distribution
  if ("HH_type" %in% colnames(metadata)) {
    hh_summary <- metadata %>%
      group_by(HH_type) %>%
      summarise(n_participants = n(), .groups = 'drop')
    
    cat("\nParticipants by household type:\n")
    print(hh_summary)
  }
  
  # Geographic spread metrics
  if (all(c("LONG", "LAT") %in% colnames(metadata))) {
    geographic_summary <- metadata %>%
      summarise(
        long_range = max(LONG, na.rm = TRUE) - min(LONG, na.rm = TRUE),
        lat_range = max(LAT, na.rm = TRUE) - min(LAT, na.rm = TRUE),
        centroid_long = mean(LONG, na.rm = TRUE),
        centroid_lat = mean(LAT, na.rm = TRUE)
      )
    
    cat("\nGeographic spread:\n")
    cat("- Longitude range:", round(geographic_summary$long_range, 3), "degrees\n")
    cat("- Latitude range:", round(geographic_summary$lat_range, 3), "degrees\n")
    cat("- Study centroid:", round(geographic_summary$centroid_long, 3), ",", 
        round(geographic_summary$centroid_lat, 3), "\n")
  }
} else {
  cat("Limited data available for detailed geographic analysis\n")
  cat("Map shows UK administrative boundaries for context\n")
}

cat("\nInterpretation notes:\n")
cat("- This supplementary figure provides geographic context for the study\n")
cat("- Point locations show spatial distribution of participating households\n")
cat("- Different colors represent different UK regions\n")
cat("- Different shapes represent different household types\n")
cat("- Geographic spread indicates study representativeness across UK\n")

# ==============================================================================
# End of script
# ==============================================================================
