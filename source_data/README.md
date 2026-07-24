# Source Data

This folder contains the source data underlying the graphs and charts in the main figures of:

**Reconstruction of SARS-CoV-2 transmissibility within households in the UK Virus Watch Study**  
Buggiotti L. et al., *Communications Medicine* (2025)  
Preprint: medRxiv [10.1101/2025.06.12.25329406](https://doi.org/10.1101/2025.06.12.25329406)

## Files

| File | Figure | Description |
|------|--------|-------------|
| `Figure1_lineage_counts.csv` | Figure 1 | Number of Virus Watch participants assigned to each SARS-CoV-2 Pango lineage, grouped by broader VOC/VUI category (Omicron, Omicron BA.2-like, Omicron BA.5-like) |
| `Figure2_pairwise_distances.csv` | Figure 2 | Pairwise phylogenetic branch distances between SARS-CoV-2 samples, categorised by tree type (within-host vs consensus) and pair type (same household vs different household). Same-household pairs are further annotated as included in transmission analyses or excluded as inferred community introductions |
| `Figure3_within_host_tree.treefile` | Figure 3 | RAxML maximum-likelihood phylogenetic tree in Newick format, inferred from the within-host diversity alignment using the 16-state substitution model with 1,000 bootstrap replicates |
| `Figure4_bottleneck_and_transmission.csv` | Figure 4 | For each of the 62 inferred household transmission pairs: per-pair transmission bottleneck estimate (Nb) and 95% confidence interval from the exact beta-binomial method, number of informative sites, TransPhylo posterior probability of direct transmission, and proportion of shared minority variants in each direction |

## File formats

- All tabular files are CSV format (comma-separated, UTF-8 encoded)
- The phylogenetic tree is in Newick format, compatible with standard viewers (FigTree, iTOL, `ape` in R, `ete3` in Python)

## Column descriptions

### Figure1_lineage_counts.csv

| Column | Description |
|--------|-------------|
| `LINEAGE` | Pango lineage assignment (e.g. XBB.1.5, BA.2.1) |
| `VOC_VUI` | UK Health Security Agency variant grouping (Omicron, Omicron BA.2-like, Omicron BA.5-like) |
| `n_participants` | Number of Virus Watch participants assigned to that lineage |

### Figure2_pairwise_distances.csv

| Column | Description |
|--------|-------------|
| `tree_type` | Phylogeny used: `within-host` (16-state model, IUPAC-coded minority variants) or `consensus` (GTR+G) |
| `pair_type` | `same_household` or `different_household` |
| `excluded_status` | For same-household pairs only: `included_transmission` (used in analyses) or `excluded_no_transmission` (community introduction) |
| `pairwise_distance` | Pairwise branch distance between the two samples (log₂ scale for plotting) |

### Figure4_bottleneck_and_transmission.csv

| Column | Description |
|--------|-------------|
| `hh2_id` | Anonymised household identifier (hh2_1 through hh2_62) |
| `n_sites` | Number of informative sites in the 3–50% minor allele frequency window after artefact filtering |
| `Nb` | Maximum-likelihood transmission bottleneck size estimate (exact beta-binomial method) |
| `CI_low` | Lower bound of 95% confidence interval on Nb |
| `CI_high` | Upper bound of 95% confidence interval on Nb |
| `posteriorProb` | TransPhylo posterior probability of direct transmission between the two household members |
| `prop_shared_a` | Proportion of donor minority variants also present in recipient |
| `prop_shared_b` | Proportion of recipient minority variants also present in donor |
| `n_a` | Number of donor minority variants |
| `n_b` | Number of recipient minority variants |
| `n_shared` | Number of minority variants shared between donor and recipient |

## Note on missing values

For pairs with fewer than two informative sites after artefact filtering, transmission bottleneck estimation is not possible; these rows contain `NA` values in the `Nb`, `CI_low`, and `CI_high` columns but retain valid TransPhylo posterior probabilities.

## Data privacy

To protect participant privacy, this source data folder contains only the numerical values plotted in each figure. Personal identifiers, sampling dates, Ct values, viral loads, and other clinical metadata have been excluded. Sample-level metadata is available upon reasonable request to the corresponding author, subject to the Virus Watch data-sharing agreements.

## Citation

If you use these data, please cite:

Buggiotti L., Torres Ortiz A., Didelot X., Yavlinsky A., Geismar C., Kovar J., Miller C., Martin Bernal L. M., Williams R., Abubakar I., Aldridge R.W., Breuer J. Reconstruction of SARS-CoV-2 transmissibility within households in the UK Virus Watch Study. *Communications Medicine* (2026).

## Contact

For questions about the data, please contact:  
Laura Buggiotti  
UCL Great Ormond Street Institute of Child Health  
laura.buggiotti@ucl.ac.uk