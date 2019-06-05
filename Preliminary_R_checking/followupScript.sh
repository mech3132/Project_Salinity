#!bin/bash sh

## Paths for ZCU
otu16=Fraser_16S/otu_table_nochloromito_col_wtaxa.biom
otu18=Fraser_18S/otu_table_col_wtaxa.biom
mf16=Fraser_16S/MF_16sFraser_noConCOL.txt
mf18=Fraser_18S/MF_18sFraser_noConCOL.txt
tree16=~/parfreylab/shared/databases/SILVA_132_by_Stilian/trees/16S_only/NR99/NR99_otus_16S.tre 
tree18=tree18.tre

##Paths for McKenzie Computer or colorado microbe
# fraser16=Fraser_16S/otu_table_nochloromito_col_wtaxa.biom
# fraser18=Fraser_18S/otu_table_col_wtaxa.biom
otu16=Fraser_16S/otu_table_nochloromito_col_wtaxa.biom
otu18=Fraser_18S/otu_table_col_wtaxa.biom
mf16=Fraser_16S/MF_16sFraser_noConCOL.txt
mf18=Fraser_18S/MF_18sFraser_noConCOL.txt
# tree16=/Users/mckenzielab/Documents/Melissa/Project_salinity/trees/SILVA_132_16S/NR99_otus_16S.tre
tree16=Fraser_16S/NR99_otus_16S.tre
tree18=tree18.tre

## Paths for MC Computer
# fraser16=../raw_otu_and_mf/Fraser_16S/OUTPUT_FILES/otu_table_nochloromito_col_wtaxa.biom
# fraser18=../raw_otu_and_mf/Fraser_16S/OUTPUT_FILES/otu_table_col_wtaxa.biom
otu16=../raw_otu_and_mf/Fraser_16S/OUTPUT_FILES/otu_table_nochloromito_col_wtaxa.biom
otu18=../raw_otu_and_mf/Fraser_18S/OUTPUT_FILES/otu_table_col_wtaxa.biom
mf16=../raw_otu_and_mf/Fraser_16S/OUTPUT_FILES/MF_16sFraser_noConCOL.txt
mf18=../raw_otu_and_mf/Fraser_18S/OUTPUT_FILES/MF_18sFraser_noConCOL.txt
tree16=../SILVA_132_files/16s/NR99_otus_16S.tre
tree18=tree18.tre

single_rarefaction.py -i $otu16 -o fraser16_OTUTable_rare1000_col.biom -d 1000
single_rarefaction.py -i $otu18 -o fraser18_OTUTable_rare1000_col.biom -d 1000


beta_diversity.py -i fraser16_OTUTable_rare1000_col.biom -m weighted_unifrac,unweighted_unifrac,bray_curtis -t $tree16 -o beta_div_COLLAPSED_16S
beta_diversity.py -i fraser18_OTUTable_rare1000_col.biom -m weighted_unifrac,unweighted_unifrac,bray_curtis -t $tree18 -o beta_div_COLLAPSED_18S