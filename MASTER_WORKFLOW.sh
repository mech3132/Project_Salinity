#!/bin/bash

# Master workflow

# First, run all biom adjustment scripts, pre-processing, collapsing, etc

sh /Users/melissachen/Documents/Masters/Project_Environmental/Project_Salinity/raw_otu_and_mf/Baltic_16S/workflow_baltic16s_27Oct2018.sh

sh /Users/melissachen/Documents/Masters/Project_Environmental/Project_Salinity/raw_otu_and_mf/Fraser_16S/workflow_fraser16s_27Oct2018.sh

sh /Users/melissachen/Documents/Masters/Project_Environmental/Project_Salinity/raw_otu_and_mf/Fraser_18S/workflow_fraser18s_27oct2018.sh

# Now, run some preliminary data checks and make PCOA plots

sh /Users/melissachen/Documents/Masters/Project_Environmental/Project_Salinity/Removing_contamined_samples/workflow_contaminant_checking.sh

sh /Users/melissachen/Documents/Masters/Project_Environmental/Project_Salinity/Preliminary_R_checking/workflow_betaplots.sh

# Now, run SALBIN

sh /Users/melissachen/Documents/Masters/Project_Environmental/Project_Salinity/SALBIN/workflow_GradBin_7June2018.sh