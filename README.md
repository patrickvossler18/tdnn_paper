Optimal Nonparametric Inference with Two-Scale Distributional Nearest Neighbors Reproduction Materials
================

This GitHub repository contains the code needed to reproduce the analyses, visualizations, and
tables for ["Optimal Nonparametric Inference with Two-Scale Distributional Nearest Neighbors"](https://arxiv.org/abs/1808.08469)
by Emre Demirkaya, Yingying Fan, Lan Gao, Jinchi Lv, Patrick Vossler, Jingbo Wang.

## Code

This folder contains the R code for the simulation settings as well as for the real data analysis.

* `setting_1.R` contains the code needed to reproduce the results in Tables 1 and 5 of the paper.
* `setting_2.R` contains the code needed to reproduce the results in Table 2 of the paper.
* `setting_3.R` contains the code needed to reproduce the results in Table 3 of the paper.
* `figures_1_2.R` contains the code needed to reproduce Figures 1 and 2 of the paper
* `real_data_analysis.R` contains the code needed to reproduce the results in Table 4 of the paper.

**Installation instructions for installing dependencies and running the replication code can be found in `code/README.md`.**

## Data
This folder contains the abalone data set from the UCI data repository used in the real data section of the paper. The data set can also be downloaded at https://archive.ics.uci.edu/ml/datasets/abalone.

## Output 

By default, all functions save output to the `output` folder. This can be changed by modifying the `OUTPUT_DIR` variable in each of the files.


