# Multi-Omics Data Integration Project

## Description
This project integrates multi-omics data using two tools: **DIABLO** and **MOFA**. The goal is to identify key patterns and potential biomarkers related to cancer progression.

## Tools
- **DIABLO**: Used for supervised integration and classification.
- **MOFA**: Used for unsupervised integration of multi-omics data.

## Data
- **DIABLO Data**: Located in `DIABLO/data/`.
  - `raw/`: Contains raw data files for DIABLO.
  - `processed/`: Contains preprocessed data ready for DIABLO analysis.
- **MOFA Data**: Located in `MOFA/data/`.
  - `raw/`: Contains raw data files for MOFA.
  - `processed/`: Contains preprocessed data ready for MOFA analysis.

## Structure
- **DIABLO/**: Contains scripts, results, and figures for DIABLO analysis.
- **MOFA/**: Contains scripts, results, and figures for MOFA analysis.

## How to Run
1. Clone this repository.
2. Install the required R packages.
3. Run the scripts in the following order:
   - `DIABLO/scripts/DIABLO.Rmd`
   - `MOFA/scripts/MOFA.Rmd`

## Results
- Integrated datasets and differential expression results are stored in the `results/` folders.
- Visualizations are saved in the `figures/` folders.

## References
- [DIABLO Documentation](https://mixomics.org/mixdiablo/)
- [MOFA Documentation](https://biofam.github.io/MOFA2/)
