# SARS-CoV-2 Household Transmission Analysis

**Reconstruction of SARS-CoV-2 transmissibility within households in the UK Virus Watch Study**

## Overview

Analysis pipeline for studying SARS-CoV-2 transmission dynamics within UK households using genomic, phylogenetic, and epidemiological approaches. Each script generates a specific figure from the publication.

**Key objectives:**

- Reconstruct transmission chains using phylogenetic analysis incorporating within-host diversity
- Estimate household transmission bottleneck sizes with rigorous artefact controls
- Validate transmission relationships against random non-household pair negative controls
- Characterise the geographic distribution across UK regions

## Repository Structure

```
SarsCov2_HH_transmission/
├── scripts/
│   ├── Figure1.R                  # Viral lineage distribution
│   ├── Figure2.R                  # Pairwise branch distance (within-host vs consensus)
│   ├── Figure3.R                  # Within-host maximum likelihood tree with bootstrap support
│   ├── Figure4.R                  # Per-pair transmission bottlenecks and TransPhylo posteriors
│   ├── Supplementary_Fig1.R       # Household distribution in the UK
│   ├── Supplementary_Fig2.R       # UK-wide lineage surveillance context
│   ├── Supplementary_Fig3.R       # Within-host vs consensus tree tanglegram
│   ├── Supplementary_Fig4.R       # Global phylogeny with Nextstrain background sequences
│   ├── Supplementary_Fig5.R       # iSNV diagnostics (counts, Ct/VL relationships, SFS)
│   ├── Supplementary_Fig6.R       # Cohort-wide artefact characterisation
│   └── Supplementary_Fig7.R       # Household vs random pair shared variants
├── source_data/
│   ├── README.md                  # Description of source data files
│   ├── Figure1_lineage_counts.csv
│   ├── Figure2_pairwise_distances.csv
│   ├── Figure3_within_host_tree.treefile
│   └── Figure4_bottleneck_and_transmission.csv
└── README.md
```

## Scripts Overview

### Main Figures

| Script      | Purpose                                                         | Methods                                                            |
| ----------- | --------------------------------------------------------------- | ------------------------------------------------------------------ |
| `Figure1.R` | Viral lineage distribution across cohort participants           | Bar chart of lineage counts grouped by VOC/VUI                     |
| `Figure2.R` | Pairwise branch distances                                       | Violin plots comparing same/different household pairs, Wilcoxon    |
| `Figure3.R` | Within-host phylogenetic tree with variant counts               | Circular tree with bootstrap-coloured nodes and household rings    |
| `Figure4.R` | Per-pair transmission bottleneck and TransPhylo posteriors      | Forest plot of Nb estimates; TransPhylo posteriors + prop_shared   |

### Supplementary Figures

| Script                 | Purpose                                                                          | Methods                                        |
| ---------------------- | -------------------------------------------------------------------------------- | ---------------------------------------------- |
| `Supplementary_Fig1.R` | Household distribution across UK regions                                         | Spatial mapping, UK boundaries                 |
| `Supplementary_Fig2.R` | UK-wide contemporaneous lineage context                                          | COG-UK / Nextstrain surveillance data          |
| `Supplementary_Fig3.R` | Tanglegram comparing within-host vs consensus tree topology (HH2 pairs)          | Cladogram tanglegram with concordance lines    |
| `Supplementary_Fig4.R` | Global phylogeny with 4,219 UK Nextstrain sequences                              | RAxML tree with categorical tip annotations    |
| `Supplementary_Fig5.R` | Within-host variant characterisation                                             | iSNV counts, Ct/VL correlations, SFS           |
| `Supplementary_Fig6.R` | Cohort-wide artefact screen — 13 recurrent low-frequency positions               | Genome track, primer overlap, lineage strat.   |
| `Supplementary_Fig7.R` | Household vs random pair shared variants (before/after artefact filter)          | Wilcoxon comparison, box/scatter plots         |

## Source Data

The `source_data/` folder contains the numerical values underlying the graphs and charts in each main figure. See `source_data/README.md` for detailed column descriptions.

## Quick Start

### Install Dependencies

```r
required_packages <- c(
  "tidyverse", "ggplot2", "ggpubr", "ape", "ggtree", "ggtreeExtra",
  "phytools", "BactDating", "TransPhylo", "sf", "cowplot", "ggnewscale"
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

**Analysis code only** — raw sequencing data and individual-level metadata are not included due to Virus Watch data-sharing agreements. Aggregated source data underlying the figures is provided in `source_data/`.

Consensus sequences generated by this study are available under ENA accession [pending].

## Key Dependencies

**External Methods:**

- Minority variant calling: [arturotorreso/scov2_withinHost](https://github.com/arturotorreso/scov2_withinHost)
- Molecular dating: [BactDating](https://github.com/xavierdidelot/BactDating)
- Transmission inference: [TransPhylo](https://github.com/xavierdidelot/TransPhylo)
- Bottleneck inference: [Bottleneck](https://github.com/weissmanlab/BB_bottleneck)

## Citation

```
@article{buggiotti2025household,
  title={Reconstruction of SARS-CoV-2 transmissibility within households in the UK Virus Watch Study},
  author={Buggiotti L., Torres Ortiz A., Didelot X., Yavlinsky A., Geismar C., Kovar J., Miller C., Martin Bernal L. M., Williams R., Abubakar I., Aldridge R.W., Breuer J.},
  journal={medRxiv},
  year={2025},
  doi={10.1101/2025.06.12.25329406}
}
```

**Additional citations required:**

- [arturotorreso/scov2_withinHost](https://github.com/arturotorreso/scov2_withinHost) for minority variant methods

---

Part of the [UK Virus Watch Study](https://www.ucl.ac.uk/epidemiology-health-care/research/epidemiology-and-public-health/research/virus-watch)
