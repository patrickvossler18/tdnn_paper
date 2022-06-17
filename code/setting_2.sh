#!/bin/bash
SETTING=setting_2.R
SETTING_NUM=2

for DIM in 3 5 10 15 20;
do
    Rscript $SETTING --dimension $DIM >&1 | tee "../output/setting_${SETTING_NUM}_p_${DIM}_$(date +"%Y_%m_%d_%I_%M_%p").log"
done