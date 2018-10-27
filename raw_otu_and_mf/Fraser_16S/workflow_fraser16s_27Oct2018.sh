#!/bin/bash


#### Executable #####
#### For generating downstream outputs from raw OTU tables ######

# Filters out chloro + mito 
# Filters out some potential contaminants
# Adds re-adds metadata to otu table


biomFRASER=MANUAL_INPUT_FILES/fraser16S_otu_table_raw.biom
fraserMF=MANUAL_INPUT_FILES/Fraser16s_mappingfile_merged.txt
taxonomy16=../../SILVA_132_files/16s/taxonomy_7_levels_wheaders.txt

tokeep=tokeep

mkdir OUTPUT_FILES

# First, let's filter with old MF to get rid  of controls etc
filter_samples_from_otu_table.py -i $biomFRASER -o OUTPUT_FILES/otu_table_nocon.biom -m $fraserMF -s ${tokeep}:Y -n 1000
# filter out chloroplast stuff
# filter_samples_from_otu_table.py -i $biomBALTIC -m $balticMF -n 1000 -o ./OTUTableandMappingFile/OTU_Table_baltic.biom --output_mapping_fp ./OTUTableandMappingFile/MF_baltic.txt
biom add-metadata --observation-metadata-fp $taxonomy16 --sc-separated taxonomy -i OUTPUT_FILES/otu_table_nocon.biom -o OUTPUT_FILES/otu_table_nocon_wtaxa.biom
biom convert -i OUTPUT_FILES/otu_table_nocon_wtaxa.biom --table-type="OTU table" --to-tsv --header-key taxonomy -o OUTPUT_FILES/OTUTable_text.txt

grep "Chloroplast" OUTPUT_FILES/OTUTable_text.txt > OUTPUT_FILES/toDelete.txt
grep "Mitochond" OUTPUT_FILES/OTUTable_text.txt >> OUTPUT_FILES/toDelete.txt
grep "Eukary" OUTPUT_FILES/OTUTable_text.txt >> OUTPUT_FILES/toDelete.txt

## Manual removal of contaminants; contaminants from contam analysis. ##
grep "GQ350231.1.1372" OUTPUT_FILES/OTUTable_text.txt >> OUTPUT_FILES/toDelete.txt
grep "JX525782.1.1434" OUTPUT_FILES/OTUTable_text.txt >> OUTPUT_FILES/toDelete.txt
grep "KY608105.1.1200" OUTPUT_FILES/OTUTable_text.txt >> OUTPUT_FILES/toDelete.txt


filter_otus_from_otu_table.py -i OUTPUT_FILES/otu_table_nocon_wtaxa.biom -e OUTPUT_FILES/toDelete.txt -o OUTPUT_FILES/otu_table_nochloromito_wtaxa.biom

# Now, let's collapse samples
collapse_samples.py -b OUTPUT_FILES/otu_table_nochloromito_wtaxa.biom -m $fraserMF --output_biom_fp OUTPUT_FILES/otu_table_nochloromito_col_wtaxa.biom --output_mapping_fp OUTPUT_FILES/MF_16sFraser_noConCOL.txt --collapse_fields 'ColRep'

