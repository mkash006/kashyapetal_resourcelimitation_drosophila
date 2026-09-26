# Resource limitation shapes the evolution of sexual dimorphism and condition dependence in *Drosophila melanogaster*

This repository contains the data and R code associated with our manuscript, **“Resource Limitation Shapes the Evolution of Sexual Dimorphism and Condition Dependence in *Drosophila melanogaster*.”**

Sexual dimorphism often depends on the resources available during development. In this study, we ask whether this relationship is only an immediate plastic response to resource limitation, or whether it can itself evolve after long-term adaptation to a resource-limited environment. We compare four replicate *D. melanogaster* populations that have evolved for more than 250 generations under larval crowding with their four ancestral control populations. Both selection regimes were reared in a common garden across high, intermediate, and low larval density.

The analysis focuses on two related questions. First, how does larval resource availability affect sexual dimorphism in adult body size? Second, does long-term adaptation to larval crowding alter this relationship? We analyse composite body size using a principal components analysis and linear mixed-effects models, and then examine individual morphological traits after correcting for body size. Standardized effect sizes are used to describe both sexual dimorphism and the response of each sex to larval density.

## Repository contents

- `data/cdsd_dat.csv` — morphological dataset used for the analyses.
- `R/analysis.R` — complete analysis script.
- `outputs/figures/` — directory where figures are written when the script is run.
- `outputs/tables/` — directory where model tables and effect-size estimates are written when the script is run.

The data contain individual measurements of femur, tibia, thorax, left wing (LW), right wing (RW), and average wing length, together with sex, larval density treatment, selection history, and ancestry block. In the original data file the treatment column was named `Treatement`; the analysis script renames this to `Treatment` immediately after reading the data, but otherwise preserves the original coding. `CU` denotes the crowding-evolved populations and `MB` denotes the control populations. Larval density treatments are `HD`, `ID`, and `LD` for high, intermediate, and low density.

## Running the analysis

The revised manuscript was analysed in R 4.4.1. To reproduce the analyses, clone or download this repository, set the repository root as the working directory, and run:

```r
source("R/analysis.R")
```

The script checks that the required packages are installed before starting. These are `car`, `DHARMa`, `dplyr`, `emmeans`, `ggplot2`, `lme4`, `patchwork`, `reshape2`, `tibble`, and `tidyr`.

The script writes all generated tables and figures to `outputs/`. It also saves `sessionInfo()` so that the package versions used in a particular run can be recorded.

## Analysis overview

We first use a PCA to combine femur, tibia, thorax, and wing lengths into a composite measure of body size. The PCA is fitted across all individuals, sexes, selection regimes, and density treatments after centering and scaling the traits. PC1 is then analysed with a linear mixed-effects model containing the full `Treatment × Sex × Selection` factorial structure and ancestry block as a random effect.

The full model is compared with reduced models using maximum likelihood because AIC values from models with different fixed effects should not be compared after REML fitting. The final inferential model is then refitted using REML, matching the analysis described in the manuscript. Model assumptions are checked using standard residual plots and simulation-based diagnostics in `DHARMa`. The script also identifies observations with standardized residuals greater than three in absolute value and repeats the model as a sensitivity analysis without those observations; the full dataset remains the primary analysis.

We use estimated marginal means and Cohen's *d* to quantify sexual dimorphism within each density treatment and selection regime, and to quantify the response of females and males across the density gradient. This separates changes in the magnitude of dimorphism from differences that arise only through residual variance or statistical power.

For the trait-wise analyses, each focal trait is corrected for body size using a PCA constructed without that focal trait. This avoids using a trait to define the body-size axis against which the same trait is corrected. As a supplementary check, femur, tibia, and thorax are also corrected using wing length as an independent body-size proxy.

Finally, the script compares the magnitude of sexual dimorphism among size-corrected traits with their sensitivity to larval density using a Spearman rank correlation. As in the manuscript, this comparison is treated as descriptive because it contains only six traits.


## Citation

If you use these data or code, please cite the associated manuscript and the archived Zenodo release. The Zenodo DOI can be added here after the GitHub repository has been connected to Zenodo and the first release has been archived.

## Authors

Mayank Kashyap and co-authors.

For questions about the analysis, please open an issue in this repository.
