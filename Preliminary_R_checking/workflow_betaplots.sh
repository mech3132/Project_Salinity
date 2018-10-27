#!/bin/bash

##### Preliminary beta div plots to check where all samples lay ####

fraser16=../raw_otu_and_mf/Fraser_16S/OUTPUT_FILES/otu_table_nochloromito_wtaxa.biom
fraser18=../raw_otu_and_mf/Fraser_18S/OUTPUT_FILES/otu_table_nocon.biom
plotscript=../R/PCOA_plots.R

filter_samples_from_otu_table.py -i $fraser16 -n 1000 -o fraser16_otutable_min1000.biom
filter_samples_from_otu_table.py -i $fraser18 -n 1000 -o fraser18_otutable_min1000.biom

beta_diversity.py -i fraser16_otutable_min1000.biom -m bray_curtis -o beta_div_fr16
beta_diversity.py -i fraser18_otutable_min1000.biom -m bray_curtis -o beta_div_fr18

Rscript $plotscript