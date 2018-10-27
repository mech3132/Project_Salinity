#!/bin/bash

#### Get unique list of OTUs
baltic <- "16sBALTIC/OTUs.txt"
fraser <- "16sFRASER/OTUs.txt"

botus <- read.delim(baltic, header=FALSE, stringsAsFactors = FALSE)
fotus <- read.delim(fraser, header=FALSE, stringsAsFactors = FALSE)

allotus <- unique(c(botus$V1,fotus$V1))

write.table(allotus, file="allOTUs_combo_BF16.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
