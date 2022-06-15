#!/bin/bash
SETTING=setting_1.R
SETTING_NUM=1

declare -a data_types=("normal" "unif")

for DATA in "${data_types[@]}"
do
    Rscript $SETTING --data_type $DATA >&1 | tee "setting_${SETTING_NUM}_p_${DATA}_$(date +"%Y_%m_%d_%I_%M_%p").log"
done