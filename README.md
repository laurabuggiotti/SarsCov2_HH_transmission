# SARS-CoV-2 Household Transmission Analysis

**Reconstruction of SARS-CoV-2 transmissibility within households in the UK Virus Watch Study**

## Overview

Analysis pipeline for studying SARS-CoV-2 transmission dynamics within UK households using genomic, phylogenetic, and epidemiological approaches. Each script generates a specific figure from the publication.

**Key objectives:**
- Estimate household transmission bottleneck sizes
- Reconstruct transmission chains using phylogenetics
- Analyze geographic distribution across UK regions
- Validate transmission relationships

## Repository Structure

```
SarsCov2_HH_transmission/
├── scripts/
│   ├── Figure1.R          # Viral lineage distribution
│   ├── Figure2.R          # Pairwise branch distance
│   ├── Figure3.R          # Maximum likelihood trees
│   ├── Figure4.R          # Transmission bottlenecks
│   ├── Supplementary_Fig1.R    # Bottleneck density plots
│   ├── Supplementary_Fig2.R    # Household distribution in the UK
│   ├── Supplementary_Fig4.R    # Tree comparison
│   ├── Supplementary_Fig6.R    # Transmission inference (posterior probabilities)
│   └── Supplementary_Fig7.R    # Inferred transmission pairs vs random pairs
└── README.md
```

## Scripts Overview

### Main Figures

| Script | Purpose | Methods |
|--------|---------|---------|
| `Figure1.R` | Viral lineage distribution | Bar charts, factor ordering |
| `Figure2.R` | Pairwise branch distances | Distance matrices, violin plots |
| `Figure3.R` | Maximum likelihood trees | Tree visualization, bootstrap support |
| `Figure4.R` | Transmission bottlenecks | Bottleneck inference, confidence intervals |

### Supplementary Figures

| Script | Purpose | Methods |
|--------|---------|---------|
| `Supplementary_Fig1.R` | Bottleneck validation | Density plots, sensitivity analysis |
| `Supplementary_Fig2.R` | Household distribution in the UK | Spatial mapping, UK boundaries |
| `Supplementary_Fig4.R` | Tree comparison | Tanglegram visualization |
| `Supplementary_Fig6.R` | Transmission inference | BactDating, TransPhylo, MCMC |
| `Supplementary_Fig7.R` | Inferred transmission pairs vs random pairs | Shared mutation analysis |

## Quick Start

### Install Dependencies
```r
required_packages <- c(
  "tidyverse", "ggplot2", "ggpubr", "ape", "ggtree", 
  "BactDating", "TransPhylo", "sf", "cowplot"
)
install.packages(required_packages)
```

### Run Analysis
```bash
# Clone repository
git clone https://github.com/laurabuggiotti/SarsCov2_HH_transmission.git
cd SarsCov2_HH_transmission

# Run specific figures
Rscript scripts/Figure1.R
Rscript scripts/Figure4.R
Rscript scripts/Supplementary_Fig7.R
```

## Data Requirements

⚠️ **Analysis code only** - Raw data not included due to privacy constraints.

## Key Dependencies

**External Methods:**
- Minority variant calling: [arturotorreso/scov2_withinHost](https://github.com/arturotorreso/scov2_withinHost)
- Molecular dating: [BactDating](https://github.com/xavierdidelot/BactDating)
- Transmission inference: [TransPhylo](https://github.com/xavierdidelot/TransPhylo)
- Bottleneck inference: [Bottleneck](https://github.com/weissmanlab/BB_bottleneck)

## Citation

```bibtex
@article{buggiotti2025household,
  title={Reconstruction of SARS-CoV-2 transmissibility within households in the UK Virus Watch Study},
  author={Buggiotti, Laura and Torres Ortiz, Arturo and Didelot, Xavier and others},
  journal={medRxiv},
  year={2025},
  doi={10.1101/2025.06.12.25329406}
}
```

**Additional citations required:**
- [arturotorreso/scov2_withinHost](https://github.com/arturotorreso/scov2_withinHost) for minority variant methods

---
Part of the [UK Virus Watch Study](https://www.ucl.ac.uk/epidemiology-health-care/research/epidemiology-and-public-health/research/virus-watch)