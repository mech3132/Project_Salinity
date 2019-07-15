#!/bin/bash


#### Executable #####
#### For generating downstream outputs from raw OTU tables ######

# Filters out chloro + mito 
# Filters out some potential contaminants
# Adds re-adds metadata to otu table

# biomFRASER=MANUAL_INPUT_FILES/OTU_Table_wtaxa.biom
biomFRASER=MANUAL_INPUT_FILES/OTU_Table_wtaxa_unassigned.biom
fraserMF=MANUAL_INPUT_FILES/MF_Fraser18s_control_readded.txt
tokeep=tokeep

mkdir OUTPUT_FILES

# First, let's filter with old MF to get rid  of controls etc
filter_samples_from_otu_table.py -i $biomFRASER -o OUTPUT_FILES/otu_table_nocon_temp.biom --sample_id_fp $fraserMF 
filter_samples_from_otu_table.py -i OUTPUT_FILES/otu_table_nocon_temp.biom -o OUTPUT_FILES/otu_table_nocon1.biom -m $fraserMF -s ${tokeep}:Y -n 1000

biom convert -i OUTPUT_FILES/otu_table_nocon1.biom --to-tsv --header-key taxonomy -o OUTPUT_FILES/otu_table_temp.txt

touch OUTPUT_FILES/toDelete.txt
grep "Unassigned" OUTPUT_FILES/otu_table_temp.txt > OUTPUT_FILES/toDelete.txt
grep "70124" OUTPUT_FILES/otu_table_temp.txt >> OUTPUT_FILES/toDelete.txt
grep "59829" OUTPUT_FILES/otu_table_temp.txt >> OUTPUT_FILES/toDelete.txt

filter_otus_from_otu_table.py -i OUTPUT_FILES/otu_table_nocon1.biom -e OUTPUT_FILES/toDelete.txt -o OUTPUT_FILES/otu_table_nocon.biom
rm OUTPUT_FILES/otu_table_nocon1.biom
rm OUTPUT_FILES/otu_table_nocon_temp.biom
rm OUTPUT_FILES/otu_table_temp.txt


# Now, let's collapse samples
collapse_samples.py -b OUTPUT_FILES/otu_table_nocon.biom -m $fraserMF --output_biom_fp OUTPUT_FILES/otu_table_col_wtaxa.biom --output_mapping_fp OUTPUT_FILES/MF_18sFraser_noConCOL.txt --collapse_fields 'ColRep'

biom convert -i OUTPUT_FILES/otu_table_col_wtaxa.biom --table-type="OTU table" --to-tsv --header-key taxonomy -o OUTPUT_FILES/OTUTable_text.txt

