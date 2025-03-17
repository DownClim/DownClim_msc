# DownClim manuscript
Sylvain Schmitt
Mar 17, 2025

[![](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![lint](https://github.com/DownClim/DownClim_msc/workflows/lint/badge.svg)](https://github.com/DownClim/DownClim_msc/actions?query=workflow%3Alint)

**DownClim_msc** is a sub-project of **DownClim** developping analyses
for the manuscript in preparation entitled: “*DownClim*: downscaling
climate models to improve ensemble forecasting in ecology”.

## Usage

**DownClim_msc** takes advantage of
[snakemake](https://snakemake.readthedocs.io/en/stable/) and
[conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
workflow for now but will be replaced by a call to **DownClim** python
package in a short future. Then all analyses rely on the quarto
documents (`files.qmd`) that can be run with R and associated
environment defined with
[renv](https://rstudio.github.io/renv/articles/renv.html).

## Project

**DownClim_msc** includes:

- [snakemake](https://snakemake.readthedocs.io/en/stable/) and
  [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
  workflow for automatic data retrieval:
  - Workflow definition in `Snakefile`
  - Worflow scripts in `scripts/`
  - Workflow configuration in `config/`
  - Worflow results in `results/`
  - conda environments in `envs/`
- Analyse of the data with associated documentation and figures:
  - Reproductive analyses in `files.qmd`
  - Resulting pages in `docs/`
  - Document structure definition in `_quarto.yml`
- All data in `data/`
- Intermediary files in `outputs/`
- Figures in `figures/`
- R environment definition with
  [renv](https://rstudio.github.io/renv/articles/renv.html) in `renv/`
  and `renv/lock`
- R files (`.Rbuildignore` , `.Rdata` , `.Rprofile` , `.Rhistory`)
- Git and GitHub files (`.gitignore` , `.github/`)
- Project documentation (`README.qmd` , `README.md` , `NEWS.md` ,
  `LICENSE`)

## Contribution

You can contribute to the project by forking the repository on github
and cloning the fork to your machine using several options, including
GitHub desktop GUI.

## Help

Please preferentially create an issue on GitHub for any questions, bugs
or help needed regarding **DownClim_msc**.

## Poeple

- Sylvain Schmitt (sylvain.schmitt@cirad.fr)
- Ghislain Vieilledent (ghislain.vieilledent@cirad.fr)
- Thomas Arsouze (thomas.arsouze@cirad.fr)
- Achille Mauri (mauri.achille@gmail.com)
