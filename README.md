# Data and Code for: Phytosociological analysis of Juniperus woodlands, NW Algeria

This repository contains the raw data matrices and the R script used for the statistical analysis presented in the manuscript submitted to the *Journal of Vegetation Science*.

## Study Context
The analysis focuses on floristic differentiation and beta diversity patterns in coastal and montane *Juniperus* woodlands of northwestern Algeria.

## Repository Structure

### 1. Data Files
All data matrices are structured with **Sites (Relevés)** as rows and **Variables/Species** as columns.

| Filename | Description |
| :--- | :--- |
| `env_final.xlsx` | **Environmental Matrix**: Contains topographic, edaphic, and bioclimatic variables for each sampled stand. |
| `ad_final.xlsx` | **Abundance-Dominance Matrix**: Floristic data with abundance coefficients (transformed). Used for quantitative analysis. |
| `pa_final.xlsx` | **Presence-Absence Matrix**: Binary floristic data (0/1). Used for beta diversity and richness calculations. |
| `inv_final.xlsx` | **Inventory Matrix**: Complete raw floristic inventory of the studied relevés. |

### 2. Analysis Script
* **File:** `script_analysis.R` (Renamed from original upload for clarity)
* **Language:** R Statistical Software
* **Key Packages Required:** `vegan`, `ade4`, `cluster`, `labdsv` (and dependencies).

## Reproducibility Workflow

To reproduce the analysis and figures presented in the article:

1.  Clone this repository or download the files.
2.  Set your R working directory to the folder containing the downloaded files.
3.  Ensure the required packages are installed:
    ```R
    install.packages(c("vegan", "ade4", "cluster", "labdsv", "readxl"))
    ```
4.  Run the script `script_analysis.R`.

## License
This data is shared under the **Creative Commons Attribution 4.0 International (CC BY 4.0)** license, allowing reuse with appropriate citation of the original article.

## Contact
For any queries regarding the data or code, please open an issue in this repository.
