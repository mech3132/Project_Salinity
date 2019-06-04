#!/bin/bash
library(MASS) # for NMDS plotting (isoMDS)
library(vegan) # for adonis
library(chron) # for date
##### This script takes otutables an stuff from salinity project and plots them to see if there are differences in extr/seq methods #######

#### Set FP #####
# setwd("~/Documents/Masters/Project_Environmental/Project_Salinity/Preliminary_R_checking/")
dm16sFP_min_BC <- "beta_div_fr16_min/bray_curtis_fraser16_otutable_min1000.txt"
dm16sFP_min_UWU <- "beta_div_fr16_min/unweighted_unifrac_fraser16_otutable_min1000.txt"
dm16sFP_min_WU <- "beta_div_fr16_min/weighted_unifrac_fraser16_otutable_min1000.txt"

dm16sFP_rare_BC <- "beta_div_fr16_rare/bray_curtis_fraser16_otutable_rare1000.txt"
dm16sFP_rare_UWU <- "beta_div_fr16_rare/unweighted_unifrac_fraser16_otutable_rare1000.txt"
dm16sFP_rare_WU <- "beta_div_fr16_rare/weighted_unifrac_fraser16_otutable_rare1000.txt"

dm18sFP_min_BC <- "beta_div_fr18_min/bray_curtis_fraser18_otutable_min1000.txt"
dm18sFP_min_UWU <- "beta_div_fr18_min/unweighted_unifrac_fraser18_otutable_min1000.txt"
dm18sFP_min_WU <- "beta_div_fr18_min/weighted_unifrac_fraser18_otutable_min1000.txt"

dm18sFP_rare_BC <- "beta_div_fr18_rare/bray_curtis_fraser18_otutable_rare1000.txt"
dm18sFP_rare_UWU <- "beta_div_fr18_rare/unweighted_unifrac_fraser18_otutable_rare1000.txt"
dm18sFP_rare_WU <- "beta_div_fr18_rare/weighted_unifrac_fraser18_otutable_rare1000.txt"

mf16FP <- "../raw_otu_and_mf/Fraser_16S/MANUAL_INPUT_FILES/Fraser16s_mappingfile_merged.txt"
mf18FP <- "../raw_otu_and_mf/Fraser_18S/MANUAL_INPUT_FILES/MF_Fraser18s_control_readded.txt"

#############
# setwd("/Users/melissachen/Documents/Masters/Project_QTAG_writing/Preliminary_R_checking") # For troubleshootin

for ( type in c("min","rare")) {
    for ( metric in c("BC","UWU","WU")) {
        dm16 <- read.delim(get(paste0("dm16sFP_",type,"_",metric)), header=TRUE, stringsAsFactors = FALSE, row.names=1)
        dm18 <- read.delim(get(paste0("dm18sFP_",type,"_",metric)), header=TRUE, stringsAsFactors = FALSE, row.names=1)
          
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
        pdf(paste0("Fr16s_NMDS_",type,"_",metric,".pdf"),7,7)
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
        
        adonis(dm16 ~ SalinityEnviron + Month + Month:Day, data=mf16, method="bray")
        mf16$juliandate <- julian(x=mf16$Month,d=mf16$Day,y=0, origin. = c(1,0,0))
        adonis(dm16 ~ SalinityEnviron*juliandate, data=mf16, method="bray")
        
        #### NMDS 18s ####
        pdf(paste0("Fr18s_NMDS_",type,"_",metric,".pdf"),7,7)
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
        
        mf16$juliandate <- julian(x=mf16$Month,d=mf16$Day,y=0, origin. = c(1,0,0))
        capture.output(
            adonis(dist(dm16) ~ SalinityEnviron*juliandate, data=mf16)
            , file=paste0("adonis_fr16_",type,"_",metric,".txt"))
        
        mf18$juliandate <- julian(x=mf18$Month,d=mf18$Day,y=0, origin. = c(1,0,0))
        capture.output(
            adonis(dist(dm18) ~ SalinityEnviron*juliandate, data=mf18)
            , file=paste0("adonis_fr18_",type,"_",metric,".txt"))
        
    }
    
}


# dm16 <- read.delim(dm16sFP, header=TRUE, stringsAsFactors = FALSE, row.names=1)
# dm18_BC <- read.delim(dm18sFP_BC, header=TRUE, stringsAsFactors = FALSE, row.names=1)
# dm18_UWU <- read.delim(dm18sFP_UWU, header=TRUE, stringsAsFactors = FALSE, row.names=1)
# dm18_WU <- read.delim(dm18sFP_WU, header=TRUE, stringsAsFactors = FALSE, row.names=1)
# 
# mf16 <- read.delim(mf16FP, header=TRUE, stringsAsFactors = FALSE, row.names=1)
# mf18 <- read.delim(mf18FP, header=TRUE, stringsAsFactors = FALSE, row.names=1)
# 
# mf16 <- mf16[match(rownames(dm16), rownames(mf16)),]
# 
# mf18 <- mf18[match(rownames(dm18_BC), rownames(mf18)),]
# mf18 <- mf18[match(rownames(dm18_UWU), rownames(mf18)),]
# mf18 <- mf18[match(rownames(dm18_WU), rownames(mf18)),]
# 
# set.seed(1234012934)
# nmds16 <- isoMDS(dist(dm16), k=2)
# 
# nmds18_BC <- isoMDS(dist(dm18_BC), k=2)
# nmds18_UWU <- isoMDS(dist(dm18_UWU), k=2)
# nmds18_WU <- isoMDS(dist(dm18_WU), k=2)
# 
# 
# # Plot via colors
# 
# salGrad <- colorRampPalette(colors = c("white","blue","darkblue"))
# colorSal <- salGrad(34)
# salRound <- round(as.numeric(mf16$SalinityEnviron))
# salRound18s <- round(as.numeric(mf18$SalinityEnviron))
# 
# #### NMDS 16s #####
# pdf("Fr16s_NMDS.pdf",7,7)
# par(mfrow=c(2,2), mar=c(4.1,4.1,2.1,2.1))
# plot(nmds16$points, xlab="NMDS1", ylab="NMDS2"
#      , pch=21
#      , bg=colorSal[factor(salRound)]
# )
# legend("bottomright",legend=c("Salinity=0","Salinity=33"),pch=21, pt.bg=c("white","darkblue"), cex=0.7)
# plot(nmds16$points, xlab="NMDS1", ylab="NMDS2"
#      , pch=21
#      , bg=c("purple","green")[factor(mf16$Extrmethod)]
# )
# legend("bottomright",legend=c("Single Tube","Plate"),pch=21, pt.bg=c("green","purple"), cex=0.7)
# plot(nmds16$points, xlab="NMDS1", ylab="NMDS2"
#      , pch=21
#      , bg=c("red","blue","orange")[factor(mf16$Year)]
# )
# legend("bottomright",legend=c("Fr2014-MC_FN","Fr2015-FN","Fr2016-MC"),pch=21, pt.bg=c("red","blue","orange"), cex=0.7)
# plot(nmds16$points, xlab="NMDS1", ylab="NMDS2"
#      , pch=21
#      , bg=c("grey","yellow")[factor(mf16$Polymerase)]
# )
# legend("bottomright",legend=c("Phusion","5Prime"),pch=21, pt.bg=c("grey","yellow"), cex=0.7)
# dev.off()
# 
# adonis(dm16 ~ SalinityEnviron + Month + Month:Day, data=mf16, method="bray")
# mf16$juliandate <- julian(x=mf16$Month,d=mf16$Day,y=0, origin. = c(1,0,0))
# adonis(dm16 ~ SalinityEnviron*juliandate, data=mf16, method="bray")
# 
# #### NMDS 18s ####
# pdf("Fr18s_NMDS_BC.pdf",7,7)
# par(mfrow=c(2,2), mar=c(4.1,4.1,2.1,2.1))
# plot(nmds18_BC$points, xlab="NMDS1", ylab="NMDS2"
#      , pch=21
#      , bg=colorSal[factor(salRound18s)]
# )
# legend("bottomright",legend=c("Salinity=0","Salinity=33"),pch=21, pt.bg=c("white","darkblue"), cex=0.7)
# plot(nmds18_BC$points, xlab="NMDS1", ylab="NMDS2"
#      , pch=21
#      , bg=c("purple","green")[factor(mf18$Extrmethod)]
# )
# legend("bottomright",legend=c("Single Tube","Plate"),pch=21, pt.bg=c("green","purple"), cex=0.7)
# plot(nmds18_BC$points, xlab="NMDS1", ylab="NMDS2"
#      , pch=21
#      , bg=c("red","blue","orange")[factor(mf18$Year)]
# )
# legend("bottomright",legend=c("Fr2014-MC_FN","Fr2015-FN","Fr2016-MC"),pch=21, pt.bg=c("red","blue","orange"), cex=0.7)
# dev.off()
# 
# pdf("Fr18s_NMDS_UWU.pdf",7,7)
# par(mfrow=c(2,2), mar=c(4.1,4.1,2.1,2.1))
# plot(nmds18_UWU$points, xlab="NMDS1", ylab="NMDS2"
#      , pch=21
#      , bg=colorSal[factor(salRound18s)]
# )
# legend("bottomright",legend=c("Salinity=0","Salinity=33"),pch=21, pt.bg=c("white","darkblue"), cex=0.7)
# plot(nmds18_UWU$points, xlab="NMDS1", ylab="NMDS2"
#      , pch=21
#      , bg=c("purple","green")[factor(mf18$Extrmethod)]
# )
# legend("bottomright",legend=c("Single Tube","Plate"),pch=21, pt.bg=c("green","purple"), cex=0.7)
# plot(nmds18_UWU$points, xlab="NMDS1", ylab="NMDS2"
#      , pch=21
#      , bg=c("red","blue","orange")[factor(mf18$Year)]
# )
# legend("bottomright",legend=c("Fr2014-MC_FN","Fr2015-FN","Fr2016-MC"),pch=21, pt.bg=c("red","blue","orange"), cex=0.7)
# dev.off()
# 
# pdf("Fr18s_NMDS_WU.pdf",7,7)
# par(mfrow=c(2,2), mar=c(4.1,4.1,2.1,2.1))
# plot(nmds18_WU$points, xlab="NMDS1", ylab="NMDS2"
#      , pch=21
#      , bg=colorSal[factor(salRound18s)]
# )
# legend("bottomright",legend=c("Salinity=0","Salinity=33"),pch=21, pt.bg=c("white","darkblue"), cex=0.7)
# plot(nmds18_WU$points, xlab="NMDS1", ylab="NMDS2"
#      , pch=21
#      , bg=c("purple","green")[factor(mf18$Extrmethod)]
# )
# legend("bottomright",legend=c("Single Tube","Plate"),pch=21, pt.bg=c("green","purple"), cex=0.7)
# plot(nmds18_WU$points, xlab="NMDS1", ylab="NMDS2"
#      , pch=21
#      , bg=c("red","blue","orange")[factor(mf18$Year)]
# )
# legend("bottomright",legend=c("Fr2014-MC_FN","Fr2015-FN","Fr2016-MC"),pch=21, pt.bg=c("red","blue","orange"), cex=0.7)
# dev.off()
# 
# 
# adonis(dm18 ~ SalinityEnviron + Year + Year:Month + Year:Month:Day, data=mf18, method="bray")
# mf18$juliandate <- julian(x=mf18$Month,d=mf18$Day,y=0, origin. = c(1,0,0))
# adonis(dm18 ~ SalinityEnviron*juliandate, data=mf18, method="bray")

###### CHECKING DELETED SAMPLES #######
