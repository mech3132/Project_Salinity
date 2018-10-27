#!/bin/bash

# Workflow to run all of the SALBIN at once

GradBinFP=../QTAG/QTAG_11aug2018.py

biomBALTIC=../raw_otu_and_mf/Baltic_16S/OUTPUT_FILES/otu_table_nochloromito_newtax.biom
MFBALTIC=../raw_otu_and_mf/Baltic_16S/MANUAL_INPUT_FILES/metadata_table.tsv
outputBALTIC=16sBALTIC
biomFRASER16=../raw_otu_and_mf/Fraser_16s/OUTPUT_FILES/otu_table_nochloromito_col_wtaxa.biom
MFFRASER16=../raw_otu_and_mf/Fraser_16s/OUTPUT_FILES/MF_16sFraser_noConCOL.txt
outputFRASER16=16sFRASER
biomFRASER18=../raw_otu_and_mf/Fraser_18s/OUTPUT_FILES/otu_table_col_wtaxa.biom
MFFRASER18=../raw_otu_and_mf/Fraser_18s/OUTPUT_FILES/MF_18sFraser_noConCOL.txt
outputFRASER18=18sFraser
tree16=../SILVA_132_files/16s/NR99_otus_16S.tre

gradient='fresh,brackish,marine'
gradientNAME=SalinityEnviron

RscriptFP=../QTAG/QTAG_graphing_15May2018.R

typetaxaFP=../R/TypesAcrossTaxonomy_V5_7june2018.R
taxatabFP=../R/RworkflowComparisons_16sfraserbaltic.R
basicmb=../R/Basic_MB_stats.R
SalPD=../R/Salinity_PD_betweentypes_9July2018.R

getOTUs=../R/make_list_all_OTUs.R
################

dobaltic='True'
dofraser16='True'
dofraser18='True'
doall='True'


if [ $dobaltic == 'True' ] 
then
	# This is for Baltic 16s

	# Make sure is \n
	tr '\r\n' '\n' < $MFBALTIC > temp.txt
	tr '\n\n' '\n' < temp.txt > $MFBALTIC
	rm temp.txt

	python $GradBinFP -t $biomBALTIC -m $MFBALTIC -M $gradientNAME --gradient $gradient -o $outputBALTIC -R ../$RscriptFP
# 	python $GradBinFP -t $biomBALTIC -m $MFBALTIC -M $gradientNAME --gradient $gradient -o TESTBALTIC_hiabundonly -R $RscriptFP --minCountTable 150

# 	cd $outputBALTIC
# # 	python $otuExtract
# 	cd ..
fi


if [ $dofraser16 == 'True' ]
then

	# This is for Fraser 16s

	# Make sure is \n
	tr '\r\n' '\n' < $MFFRASER16 > temp.txt
	tr '\n\n' '\n' < temp.txt > $MFFRASER16
	rm temp.txt

	python $GradBinFP -t $biomFRASER16 -m $MFFRASER16 -M $gradientNAME --gradient $gradient -o $outputFRASER16 -R ../$RscriptFP
# 	python $GradBinFP -t $biomFRASER16 -m $MFFRASER16 -M $gradientNAME --gradient $gradient -o TESTFRASER_hiabundonly -R $RscriptFP --minCountTable 2000

# 	cd $outputFRASER16
# # 	python $otuExtract
# 	cd ..
fi


if [ $dofraser18 == 'True' ]
then

# This is for Fraser 18s
	# Make sure is \n
	tr '\r\n' '\n' < $MFFRASER18 > temp.txt
	tr '\n\n' '\n' < temp.txt > $MFFRASER18
	rm temp.txt

	python $GradBinFP -t $biomFRASER18 -m $MFFRASER18 -M $gradientNAME --gradient $gradient -o $outputFRASER18 -R ../$RscriptFP
# 	cd $outputFRASER18
# # 	python $otuExtract
# 	cd ..
fi 

# This is to compareFraser with Baltic
# mkdir COMPARISONS
# cd COMPARISONS

# python $comparingGradBin -n Baltic,Fraser -m ../${outputBALTIC}/modelBoundaries_type.txt,../${outputFRASER16}/modelBoundaries_type.txt -N ../${outputBALTIC}/gradient.txt -t $taxamap16



if [ $doall == 'True' ] 
then
	### Now, get list of otus
	
	Rscript $getOTUs
	
	### make trees
	
	# For Baltic
	filter_tree.py -i $tree16 -t 16sBALTIC/OTUs.txt -o tree_B16_filt.tre

	# For Fraser
	filter_tree.py -i $tree16 -t 16sFRASER/OTUs.txt -o tree_F16_filt.tre

	# For Fraser 18s
	# filter_tree.py -i /Users/melissachen/Documents/Masters/Project_Environmental/FromBotaClust_1feb2018/trees/SILVA_132/18s/*.tre -t /Users/melissachen/Documents/Masters/Project_Environmental/FromBotaClust_1feb2018/SALBIN_7June2018_vsearch_silva132/18sFraser_7june2018/OTUs.txt -o tree_F18_filt.tre
	# first, add "BACTERIA" to list of OTUs to keep
	echo "\r\n" >> 18sFraser/OTUs.txt
	echo "BACTERIA" >> 18sFraser/OTUs.txt

	filter_tree.py -i ../raw_otu_and_mf/Fraser_18S/MANUAL_INPUT_FILES/tmp_withbact.tre -t 18sFraser/OTUs.txt -o tree_F18_filt.tre


	# For combo baltic and fraser 16s; **** must already have made non-duplicate OTU list from both datasets
	filter_tree.py -i $tree16 -o tree_BF16_filt.tre -t allOTUs_combo_BF16.txt


	# Run remaining scripts
	Rscript $basicmb
	Rscript $taxatabFP
	Rscript $typetaxaFP
	Rscript $SalPD

fi



