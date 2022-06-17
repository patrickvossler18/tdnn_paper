#!/bin/bash
SETTING=real_data_analysis.R

for seed_val in {1234..1724..10}; 
do
  Rscript $SETTING --seed_val $seed_val >&1 | tee "../output/real_data_analysis_seed_val_${seed_val}_$(date +"%Y_%m_%d_%I_%M_%p").log"
done

Rscript real_data_analysis_results.R