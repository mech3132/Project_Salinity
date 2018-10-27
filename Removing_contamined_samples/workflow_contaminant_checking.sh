#!/bin/bash 

badsampleIDscript=../R/ID_bad_samples_7aug2018.R
scanScript=../R/Scan_for_contam_VProjSal.R
otu16=../raw_otu_and_mf/Fraser_16S/MANUAL_INPUT_FILES/fraser16S_otu_table_raw.biom
otu18=../raw_otu_and_mf/Fraser_18S/MANUAL_INPUT_FILES/OTU_Table_wtaxa.biom
mf16=../raw_otu_and_mf/Fraser_16S/MANUAL_INPUT_FILES/Fraser16s_mappingfile_merged.txt
mf18=../raw_otu_and_mf/Fraser_18S/MANUAL_INPUT_FILES/MF_Fraser18s_control_readded.txt
# mf18=
#taxonomy
taxonomy16=../SILVA_132_files/16s/taxonomy_7_levels_wheaders.txt
taxonomy18=../SILVA_132_files/18s/taxonomy_7_levelsw_headers.txt

filter_samples_from_otu_table.py -i $otu16 -n 1000 -o fraser16_otutable_min1000.biom
filter_samples_from_otu_table.py -i $otu18 -n 1000 -o fraser18_otutable_min1000.biom

beta_diversity.py -i fraser16_otutable_min1000.biom -m bray_curtis -o beta_div_fr16
beta_diversity.py -i fraser18_otutable_min1000.biom -m bray_curtis -o beta_div_fr18

Rscript $badsampleIDscript

# First, add the taxonomy onto the otu table
biom add-metadata --observation-metadata-fp $taxonomy16 --sc-separated taxonomy -i $otu16 -o otu_table_wtaxa_16.biom
biom add-metadata --observation-metadata-fp $taxonomy18 --sc-separated taxonomy -i $otu18 -o otu_table_wtaxa_18.biom

# change into text format
biom convert -i otu_table_wtaxa_16.biom --header-key taxonomy --to-tsv -o otu_table_wtaxa_16_text.txt
biom convert -i otu_table_wtaxa_18.biom --header-key taxonomy --to-tsv -o otu_table_wtaxa_18_text.txt

# First, let's identify some potential contaminants
Rscript $scanScript -i otu_table_wtaxa_16_text.txt -m $mf16 -b control:y --filter_low_abund_per_sample 5 -o filter_contam_16 --dashes TRUE --numberstart FALSE -s Year:2014,2015 --filter_low_abund_overall 10 --delim ';'

# Then, let's delete low-abundance OTUs
# t is min reads in sample for OTU to be kept; otherwise set to 0
# n is min reads in sample for sample to be kept.
# python /Users/melissachen/Documents/Masters/UNIVERSALCODE_git/MSC/Delete_low_abundance_v2.py -i OTU_Table_wtaxa.biom -t 3 -n 2000 -o otu_table_wtaxa_16S_filtered_t3T1n2000

