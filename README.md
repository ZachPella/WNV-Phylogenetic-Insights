# West Nile Virus Phylogenetic Insights: Nebraska 2023 & Global Context

This repository presents a phylogenetic analysis of West Nile Virus (WNV) genomes, with a particular focus on sequences collected in **Nebraska during 2023**. It also includes a global overview of WNV genome availability over time, providing context for the dataset used in this study.

The primary goals of this project are:
1.  **Visualize the distribution of WNV genome submissions over time**, highlighting the rationale for focusing on data up to a specific year.
2.  **Generate a comprehensive phylogenetic tree**, emphasizing the genomic diversity of WNV strains from Nebraska in 2023 within a broader global context.

---

## Project Components & Visualizations

### 1. WNV Genome Counts by Year

This section analyzes the annual submission rates of WNV genomes to public databases. The generated bar graph demonstrates a significant drop-off in recent years, which informed our decision to set a **2019 cutoff** for the broader phylogenetic analysis. This ensures that our global tree is built upon a robust and relatively complete dataset, while still allowing for the detailed inclusion of our newly sequenced 2023 Nebraska samples.

* **Rationale:** The visualization clearly illustrates the decline in available data post-2019, justifying why a full phylogenetic analysis incorporating all subsequent years would be less representative of global diversity.
* **Key Insight:** There's a notable disparity in genome availability pre- and post-2019.

---

### 2. Phylogenetic Tree of WNV Genomes (Nebraska 2023 Highlighted)

A large-scale phylogenetic tree constructed using `baltic` is presented. This tree includes a global collection of WNV genomes, with a special emphasis on samples collected in Nebraska during 2023.

* **Nebraska 2023 Samples:** Tips corresponding to Nebraska 2023 WNV genomes are distinctly colored by their specific **Nebraska region** (Central, West, East) and are made visually prominent (larger points, bold colors).
* **Global Context:** All other WNV genomes in the tree are displayed in muted tones (e.g., light grey) to provide essential phylogenetic context without overshadowing the Nebraska 2023 samples. Branches leading to highlighted Nebraska samples are also subtly emphasized to trace their lineage.
* **Methodology:** The tree leverages absolute time (decimal years) for its x-axis, providing a temporal perspective on WNV evolution.

---

## Data

* `tree_2025.nwk`: The Newick format phylogenetic tree file used for visualization.
* `updated_metadata.tsv`: Tab-separated file containing metadata for all WNV genomes, including collection dates and geographic regions, which are crucial for coloring and filtering the tree and bar graphs.

---

## Setup and Usage

To replicate the visualizations and analysis from this repository, you'll need a Python environment with the following libraries installed:

* `baltic` (for phylogenetic tree visualization)
* `pandas` (for data manipulation)
* `matplotlib` (for plotting)
* `numpy` (for numerical operations)

You can install these dependencies using pip:

```bash
pip install baltic pandas matplotlib numpy
