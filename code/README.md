# Replication code
## Installing dependencies
This repository uses [`renv`](https://rstudio.github.io/renv/index.html) to manage R dependencies. For those not familiar with `renv`, we recommend reading through the [introductory vignette](https://rstudio.github.io/renv/articles/renv.html). The package can be installed with the command

``` r 
install.packages("renv")
```

To install the required packages, we suggest opening the R project using RStudio because this will cause `renv` to automatically bootstrap itself, thereby downloading and installing the appropriate version of `renv` into the project library. After the project has been launched, run `renv::restore()` to install of the required packages to the project library.

For users using other IDEs, run the `renv::restore()` command in the `tdnn_paper` directory after installing the `renv` package.

### TDNN R package
Our implementation of the algorithms described in our paper are available in the `tdnn_package` repository and the code version of the code used to generate our results is available in source and binary versions at https://github.com/patrickvossler18/tdnn_package/releases/tag/v0.1.0

Note: The TDNN R package is included in the `renv` dependencies and will be installed along with other required packages by default when installing the dependencies using `renv`.

## Code execution options

The files for generating the tables in our paper (`setting_1.R`, `setting_2.R`, `setting_3.R`) can be run on the command line using the `Rscript` command and by providing command line arguments. For example, to generate the results from the first simulation setting using normally distributed covariates we can run:

``` bash
Rscript setting_1.R --data_type normal
```

We can return a list of command line arguments by running

``` bash
Rscript setting_1.R --help
```

We also provide bash scripts for each setting which loop through the parameter combinations reported in the paper. As an example, the bash script for generating the results for normally distributed covariates and uniformly distributed covariates can be executed in the command line with the command

``` bash
bash setting_1.sh
```

Finally, the code can be run interactively by replacing the variables set by the command line argument parser with desired values. For the first simulation setting this would mean changing the `data_type <- opt$data_type` to one of the two possible values `unif` or `normal`.

## Real data analysis replication

We repeat the data splitting procedure implemented in `real_data_analysis.R` for 50 times and obtain averaged prediction errors over those 50 repetitions. To increase computing efficiency, we separate the 50 repetitions into 50 parallel computation, and the random seeds are chosen from `seq(from=1234, to = 1724, by = 10)`. In our computation, the following procedure was run in parallel in 50 different R sessions and used a different seed from seq(from=1234, to = 1724, by = 10) for different R sessions, then we take average of the results from 50 parallel R sessions.

We provide a bash script `real_data_analysis.R` for automating the method described above by looping over the sequence of random seed values, calculating the estimate for each of the repetitions and saving the result for each random seed value to a CSV file `output/real_data_analysis_results.csv`. Finally, the results are averaged together in `real_data_analysis_results.R`.