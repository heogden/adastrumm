# Fat-data analysis

This directory contains the R script used to reproduce the figures from the body-fat data analysis in the AdaStruMM paper.

The analysis uses the `fat` dataset from the [`ALA`](https://rdrr.io/rforge/ALA/) package. It fits the AdaStruMM model and the comparison model described in the paper, then produces the fitted-trajectory, population-average and derivative plots.

## Running the analysis

The script should be run from the root of the `adastrumm` repository. First install the current version of the package and its dependencies, then run the code in `reproduce_fat_analysis.R`.

The analysis requires the following R packages:

- `adastrumm`
- `ALA`
- `dplyr`
- `ggplot2`
- `nlme`
- `tidyr`
- `tikzDevice`

The model fitting may take some time.

## Output

The generated TikZ files used in the paper are included in this
repository, in `paper/fat_analysis/figures`.  Running the script will
overwrite them.
