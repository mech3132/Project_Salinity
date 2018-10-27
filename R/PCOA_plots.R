#!/bin/bash
library(MASS) # for NMDS plotting (isoMDS)
##### This script takes otutables an stuff from salinity project and plots them to see if there are differences in extr/seq methods #######

#### Set FP #####
# setwd("~/Documents/Masters/Project_Environmental/Project_Salinity/Preliminary_R_checking/")
dm16sFP <- "beta_div_fr16/bray_curtis_fraser16_otutable_min1000.txt"
dm18sFP <- "beta_div_fr18/bray_curtis_fraser18_otutable_min1000.txt"
mf16FP <- "../raw_otu_and_mf/Fraser_16S/MANUAL_INPUT_FILES/Fraser16s_mappingfile_merged.txt"
mf18FP <- "../raw_otu_and_mf/Fraser_18S/MANUAL_INPUT_FILES/MF_Fraser18s_control_readded.txt"

#############

dm16 <- read.delim(dm16sFP, header=TRUE, stringsAsFactors = FALSE, row.names=1)
dm18 <- read.delim(dm18sFP, header=TRUE, stringsAsFactors = FALSE, row.names=1)
mf16 <- read.delim(mf16FP, header=TRUE, stringsAsFactors = FALSE, row.names=1)
mf18 <- read.delim(mf18FP, header=TRUE, stringsAsFactors = FALSE, row.names=1)

mf16 <- mf16[match(rownames(dm16), rownames(mf16)),]
mf18 <- mf18[match(rownames(dm18), rownames(mf18)),]

set.seed(1234012934)
nmds16 <- isoMDS(dist(dm16), k=2)
nmds18 <- isoMDS(dist(dm18), k=2)

# Plot via colors

salGrad <- colorRampPalette(colors = c("white","blue","darkblue"))
colorSal <- salGrad(34)
salRound <- round(as.numeric(mf16$SalinityEnviron))
salRound18s <- round(as.numeric(mf18$SalinityEnviron))

#### NMDS 16s #####
pdf("Fr16s_NMDS.pdf",7,7)
par(mfrow=c(2,2), mar=c(4.1,4.1,2.1,2.1))
plot(nmds16$points, xlab="NMDS1", ylab="NMDS2"
     , pch=21
     , bg=colorSal[factor(salRound)]
)
legend("bottomright",legend=c("Salinity=0","Salinity=33"),pch=21, pt.bg=c("white","darkblue"), cex=0.7)
plot(nmds16$points, xlab="NMDS1", ylab="NMDS2"
     , pch=21
     , bg=c("purple","green")[factor(mf16$Extrmethod)]
)
legend("bottomright",legend=c("Single Tube","Plate"),pch=21, pt.bg=c("green","purple"), cex=0.7)
plot(nmds16$points, xlab="NMDS1", ylab="NMDS2"
     , pch=21
     , bg=c("red","blue","orange")[factor(mf16$Year)]
)
legend("bottomright",legend=c("Fr2014-MC_FN","Fr2015-FN","Fr2016-MC"),pch=21, pt.bg=c("red","blue","orange"), cex=0.7)
plot(nmds16$points, xlab="NMDS1", ylab="NMDS2"
     , pch=21
     , bg=c("grey","yellow")[factor(mf16$Polymerase)]
)
legend("bottomright",legend=c("Phusion","5Prime"),pch=21, pt.bg=c("grey","yellow"), cex=0.7)
dev.off()

#### NMDS 18s ####
pdf("Fr18s_NMDS.pdf",7,7)
par(mfrow=c(2,2), mar=c(4.1,4.1,2.1,2.1))
plot(nmds18$points, xlab="NMDS1", ylab="NMDS2"
     , pch=21
     , bg=colorSal[factor(salRound18s)]
)
legend("bottomright",legend=c("Salinity=0","Salinity=33"),pch=21, pt.bg=c("white","darkblue"), cex=0.7)
plot(nmds18$points, xlab="NMDS1", ylab="NMDS2"
     , pch=21
     , bg=c("purple","green")[factor(mf18$Extrmethod)]
)
legend("bottomright",legend=c("Single Tube","Plate"),pch=21, pt.bg=c("green","purple"), cex=0.7)
plot(nmds18$points, xlab="NMDS1", ylab="NMDS2"
     , pch=21
     , bg=c("red","blue","orange")[factor(mf18$Year)]
)
legend("bottomright",legend=c("Fr2014-MC_FN","Fr2015-FN","Fr2016-MC"),pch=21, pt.bg=c("red","blue","orange"), cex=0.7)
dev.off()


###### CHECKING DELETED SAMPLES #######
