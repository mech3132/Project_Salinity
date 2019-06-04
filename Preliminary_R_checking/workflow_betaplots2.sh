#!/bin/bash

##### Preliminary beta div plots to check where all samples lay ####

## Paths for McKenzie Computer
fraser16=Fraser_16S/otu_table_nochloromito_wtaxa.biom
fraser18=Fraser_18S/otu_table_nocon.biom
tree16=/Users/mckenzielab/Documents/Melissa/Project_salinity/trees/SILVA_132_16S/NR99_otus_16S.tre
tree18=melissa_18s_epa_placement/raxml_epa/RAxML_labelledTree.EPA_placement_18s.figtree.tre

## Paths for MC Computer
fraser16=../raw_otu_and_mf/Fraser_16S/OUTPUT_FILES/otu_table_nochloromito_wtaxa.biom
fraser18=../raw_otu_and_mf/Fraser_18S/OUTPUT_FILES/otu_table_nocon.biom
tree16=../SILVA_132_files/16s/NR99_otus_16S.tre
tree18=../SILVA_132_files/melissa_18s_epa_placement/raxml_epa/RAxML_labelledTree.EPA_placement_18s.figtree.tre

## Paths for ACU
fraser16=Fraser_16S/otu_table_nochloromito_wtaxa.biom
fraser18=Fraser_18S/otu_table_nocon.biom
tree16=~/parfreylab/shared/databases/SILVA_132_by_Stilian/trees/16S_only/NR99/NR99_otus_16S.tre 
tree18=melissa_18s_epa_placement/raxml_epa/RAxML_labelledTree.EPA_placement_18s.figtree.tre


plotscript=PCOA_plots.R

sed -e 's/QUERY___//g' $tree18 > tree18.tre
# sed -e 's/QUERY___//g' melissa_18s_epa_placement/raxml_epa/RAxML_labelledTree.EPA_placement_18s.figtree.tre > tree18.tre

filter_samples_from_otu_table.py -i $fraser16 -n 1000 -o fraser16_otutable_min1000.biom
filter_samples_from_otu_table.py -i $fraser18 -n 1000 -o fraser18_otutable_min1000.biom

single_rarefaction.py -i fraser16_otutable_min1000.biom -d 1000 -o fraser16_otutable_rare1000.biom
single_rarefaction.py -i fraser18_otutable_min1000.biom -d 1000 -o fraser18_otutable_rare1000.biom


beta_diversity.py -i fraser16_otutable_min1000.biom -t $tree16 -m bray_curtis,weighted_unifrac,unweighted_unifrac -o beta_div_fr16_min
beta_diversity.py -i fraser16_otutable_rare1000.biom -t $tree16 -m bray_curtis,weighted_unifrac,unweighted_unifrac -o beta_div_fr16_rare

beta_diversity.py -i fraser18_otutable_min1000.biom -t tree18.tre -m bray_curtis,weighted_unifrac,unweighted_unifrac -o beta_div_fr18_min
beta_diversity.py -i fraser18_otutable_rare1000.biom -t tree18.tre -m bray_curtis,weighted_unifrac,unweighted_unifrac -o beta_div_fr18_rare

Rscript $plotscript
