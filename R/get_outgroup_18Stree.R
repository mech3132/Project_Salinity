#!bin/bash RScript
#### GET OUTGROUP ####


library(phytools)
library(picante)
library(ape)

tree <- read.tree("/Users/melissachen/Documents/Masters/Project_QTAG_writing/SALBIN/tree_F18_filt.tre")
taxaleg <- read.delim("/Users/melissachen/Documents/Masters/Project_QTAG_writing/SALBIN/18sFraser/taxaIDLegend.txt", header = F)
colnames(taxaleg) <-c("OTUID","taxonomy")

taxaleg <- taxaleg[taxaleg$OTUID %in% tree$tip.label,]

taxa.split <- sapply(taxaleg$taxonomy, FUN = function(x) {strsplit(as.character(x), split="; __", fixed = TRUE)})

class_list2 <- vector(length=length(taxa.split))
for ( i in 1:length(taxa.split) ) {
    if ( is.na(taxa.split[[i]][2]) ) {
        class_list2[i] <- taxa.split[[i]][1]
    } else {
        class_list2[i] <- taxa.split[[i]][2]
    }
}

names(class_list2) <- as.character(taxaleg$OTUID)
# unique(class_list2)
opistikonts <- names(class_list2)[grep(pattern="Opisthokonta", x=class_list2)]

opist_node <- findMRCA(tree, tips=opistikonts)



