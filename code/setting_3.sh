#!/bin/bash
SETTING=setting_3.R
SETTING_NUM=3

for DIM in in 3 5 10 15 20;
do
    Rscript $SETTING --dimension $DIM >&1 | tee "setting_${SETTING_NUM}_p_${DIM}_$(date +"%Y_%m_%d_%I_%M_%p").log"
done