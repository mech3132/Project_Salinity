#!/bin/bash

#### Executable #####
#### For generating downstream outputs from raw OTU tables ######

# Filters out chloro + mito 
# Adds re-adds metadata to otu table

biomBALTIC=MANUAL_INPUT_FILES/otu_table_absolute_unfiltered.biom
balticMF=MANUAL_INPUT_FILES/metadata_table.tsv
taxonomy16=../../SILVA_132_files/16s/taxonomy_7_levels_wheaders.txt

mkdir OUTPUT_FILES

# convert to text format
biom convert -i $biomBALTIC --to-tsv --header-key taxonomy -o OUTPUT_FILES/otu_table_absolute_unfiltered.txt

grep "Chloroplast" OUTPUT_FILES/otu_table_absolute_unfiltered.txt > OUTPUT_FILES/toDelete.txt
grep "Mitochond" OUTPUT_FILES/otu_table_absolute_unfiltered.txt >> OUTPUT_FILES/toDelete.txt
grep "Eukary" OUTPUT_FILES/otu_table_absolute_unfiltered.txt >> OUTPUT_FILES/toDelete.txt

filter_otus_from_otu_table.py -i $biomBALTIC -e OUTPUT_FILES/toDelete.txt -o OUTPUT_FILES/otu_table_nochloromito.biom

# Also, re-add names
biom add-metadata --observation-metadata-fp $taxonomy16 --sc-separated taxonomy -i OUTPUT_FILES/otu_table_nochloromito.biom -o OUTPUT_FILES/otu_table_nochloromito_newtax.biom