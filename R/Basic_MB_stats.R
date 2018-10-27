#!/bin/bash Rscript

##### Get some basic stats from the modelBoundaries files #####
setwd('../SALBIN')

mbB16FP="16sBALTIC/modelBoundaries_type.txt"
mbF16FP="16sFRASER/modelBoundaries_type.txt"
mbF18FP="18sFraser/modelBoundaries_type.txt"

mbB16 = read.delim(mbB16FP, header=TRUE, row.names=1, stringsAsFactors = FALSE)
mbF16 = read.delim(mbF16FP, header=TRUE, row.names=1, stringsAsFactors = FALSE)
mbF18 = read.delim(mbF18FP, header=TRUE, row.names=1, stringsAsFactors = FALSE)

# get rid of noclass
mbB16.filt = mbB16[mbB16$typeSimple!='noclass',]
mbF16.filt = mbF16[mbF16$typeSimple!='noclass',]
mbF18.filt = mbF18[mbF18$typeSimple!='noclass',]

# B16
B16.ranges = matrix(ncol=2,nrow=nrow(mbB16.filt),dimnames = list(rownames(mbB16.filt),c("first","second")))
for ( r in 1:nrow(B16.ranges) ) {
    temprange = mbB16.filt[r,c("boundaries","boundariestwo")]
    temptype = mbB16.filt[r,"typeSimple"]
    if ( temptype == "freshRestricted" ) {
        B16.ranges[r,1] = 0
        B16.ranges[r,2] = as.numeric(temprange[1])
    } else if ( temptype == "marineRestricted" ) {
        B16.ranges[r,1] = as.numeric(temprange[1])
        B16.ranges[r,2] = 33
    } else if (temptype == "brackishRestricted" ) {
        B16.ranges[r,] = as.numeric(temprange)
    }
}
B16.ranges = as.data.frame(B16.ranges)
B16.ranges$rlen = abs(B16.ranges$first-B16.ranges$second)

# F16
F16.ranges = matrix(ncol=2,nrow=nrow(mbF16.filt),dimnames = list(rownames(mbF16.filt),c("first","second")))
for ( r in 1:nrow(F16.ranges) ) {
    temprange = mbF16.filt[r,c("boundaries","boundariestwo")]
    temptype = mbF16.filt[r,"typeSimple"]
    if ( temptype == "freshRestricted" ) {
        F16.ranges[r,1] = 0
        F16.ranges[r,2] = as.numeric(temprange[1])
    } else if ( temptype == "marineRestricted" ) {
        F16.ranges[r,1] = as.numeric(temprange[1])
        F16.ranges[r,2] = 33
    } else if (temptype == "brackishRestricted" ) {
        F16.ranges[r,] = as.numeric(temprange)
    }
}
F16.ranges = as.data.frame(F16.ranges)
F16.ranges$rlen = abs(F16.ranges$first-F16.ranges$second)

# F18
F18.ranges = matrix(ncol=2,nrow=nrow(mbF18.filt),dimnames = list(rownames(mbF18.filt),c("first","second")))
for ( r in 1:nrow(F18.ranges) ) {
    temprange = mbF18.filt[r,c("boundaries","boundariestwo")]
    temptype = mbF18.filt[r,"typeSimple"]
    if ( temptype == "freshRestricted" ) {
        F18.ranges[r,1] = 0
        F18.ranges[r,2] = as.numeric(temprange[1])
    } else if ( temptype == "marineRestricted" ) {
        F18.ranges[r,1] = as.numeric(temprange[1])
        F18.ranges[r,2] = 33
    } else if (temptype == "brackishRestricted" ) {
        F18.ranges[r,] = as.numeric(temprange)
    }
}
F18.ranges = as.data.frame(F18.ranges)
F18.ranges$rlen = abs(F18.ranges$first-F18.ranges$second)

maxrange = max(c(B16.ranges$rlen, F16.ranges$rlen, F18.ranges$rlen))

#baltic range limits
print(c(min(B16.ranges$rlen),max(B16.ranges$rlen)))
print(median(B16.ranges$rlen))
#fraser range limits
print(c(min(F16.ranges$rlen),max(F16.ranges$rlen)))
print(median(F16.ranges$rlen))
#fraser18 range limits
print(c(min(F18.ranges$rlen), max(F18.ranges$rlen)))
print(median(F18.ranges$rlen))



pdf("rangesizes.pdf",5,10)
# quartz(,5,10)
par(mfrow=c(3,1))
hist(B16.ranges$rlen, xlim=c(0,maxrange), main="(A) Baltic Sea (16S) salinity tolerance range", xlab="Salinity Tolerance Range")
abline(v=median(B16.ranges$rlen), col="red", lty=3, sub=paste0("Median: ",median(B16.ranges$rlen)))
hist(F16.ranges$rlen, xlim=c(0,maxrange), main="(B) Fraser River Estuary (16S) salinity tolerance range", xlab="Salinity Tolerance Range")
abline(v=median(F16.ranges$rlen), col="red", lty=3, sub=paste0("Median: ",median(F16.ranges$rlen)))
hist(F18.ranges$rlen, xlim=c(0,maxrange), main="(C) Fraser River Estuary (18S) salinity tolerance range", xlab="Salinity Tolerance Range")
abline(v=median(F18.ranges$rlen), col="red", lty=3, sub=paste0("Median: ",median(F18.ranges$rlen)))
dev.off()

wilcox.test(B16.ranges$rlen, F16.ranges$rlen)
wilcox.test(F18.ranges$rlen, F16.ranges$rlen)
wilcox.test(B16.ranges$rlen, F18.ranges$rlen)

#count number in each
sum(mbB16$typeSimple=="brackishRestricted")/nrow(mbB16) # brackish Baltic
sum(mbB16$typeSimple=="freshRestricted")/nrow(mbB16) # fresh Baltic
sum(mbB16$typeSimple=="marineRestricted")/nrow(mbB16) # marine Baltic
sum(mbB16$typeSimple=="noclass")/nrow(mbB16) # marine Baltic
nrow(mbB16)


sum(mbF16$typeSimple=="brackishRestricted")/nrow(mbF16) # brackish Fraser
sum(mbF16$typeSimple=="freshRestricted")/nrow(mbF16) # fresh Fraser
sum(mbF16$typeSimple=="marineRestricted")/nrow(mbF16) # marine Fraser
sum(mbF16$typeSimple=="noclass")/nrow(mbF16) # noclass Baltic
nrow(mbF16)

sum(mbF18$typeSimple=="brackishRestricted")/nrow(mbF18) # brackish Fraser
sum(mbF18$typeSimple=="freshRestricted")/nrow(mbF18) # fresh Fraser
sum(mbF18$typeSimple=="marineRestricted")/nrow(mbF18) # marine Fraser
sum(mbF18$typeSimple=="noclass")/nrow(mbF18) # noclass Fraser
nrow(mbF18)


## Finding mid-point in tolerance range

F16.ranges$mid = apply(F16.ranges[,c(1,2)], MARGIN = 1, FUN =mean)
B16.ranges$mid = apply(B16.ranges[,c(1,2)], MARGIN = 1, FUN =mean)
F18.ranges$mid = apply(F18.ranges[,c(1,2)], MARGIN = 1, FUN =mean)

maxmid = max(c(F16.ranges$mid,B16.ranges$mid,F18.ranges$mid))+2

pdf("midpoints.pdf",5,10)
# quartz(,5,10)
par(mfrow=c(3,1))
hist(B16.ranges$mid, xlim=c(0,maxmid), main="Baltic 16s range midpoint", xlab="Midpoint of tolerance range", breaks = 20)
hist(F16.ranges$mid, xlim=c(0,maxmid), main="Fraser 16s range midpoint", xlab="Midpoint of tolerance range", breaks = 20)
hist(F18.ranges$mid, xlim=c(0,maxmid), main="Fraser 18s range midpoint", xlab="Midpoint of tolerance range", breaks = 20)
dev.off()


##### Comparing shared OTUs ######

####### Shared otus by classification (complex) #####
commonOTUs <- rownames(mbB16[na.omit(match(rownames(mbF16), rownames(mbB16))),])

combinedmb <- matrix(ncol=6, nrow=length(commonOTUs),dimnames=list(commonOTUs,c("btype","ftype","bb1","bb2","fb1","fb2")))
for ( otu in commonOTUs ) {
    ftemp <- as.vector(mbF16[otu,c('type','boundaries','boundariestwo')])
    btemp <- as.vector(mbB16[otu,c('type','boundaries','boundariestwo')])
    combinedmb[otu,] <- c(as.character(btemp[1]),as.character(ftemp[1]),as.numeric(btemp[c(2,3)]), as.numeric(ftemp[c(2,3)]))
}

sum(combinedmb[,c("btype")] == combinedmb[,c("ftype")])/nrow(combinedmb)

categories <- c("freshBloom" # 1
                ,"freshRestricted" #1
                ,"freshPeak" #2
                ,"brackishPeakLoToler" #4
                ,"brackishBloom" #5
                ,"brackishRestricted" #5
                ,"brackishPeakHiToler" #6
                ,"marinePeak" #8
                ,"marineRestricted" #9
                ,"marineBloom") #9
posCat <- c(1,1,2,4,5,5,6,8,9,9)

# lcat <- length(categories)
lcat <- max(posCat)

fpos <- sapply(combinedmb[,c("ftype")],function(x) {
    posCat[match(x,categories)]
})
bpos <- sapply(combinedmb[,c("btype")],function(x) {
    posCat[match(x,categories)]
})

combpos <- as.data.frame(cbind(bpos,fpos, comb=paste0(bpos,"-",fpos)))
# combpos[which(combpos$comb == '9-4'),]
counts_combpos <- table(combpos[,'comb'])
maxCount <- max(counts_combpos)
minCount <- min(counts_combpos)
dividend <- (maxCount)/40
thickness_combpos <- counts_combpos/dividend

pdf("shared_otus_lineplot_complex.pdf")
# quartz()
par(mar=c(4.2,10.2,2.2,10.2))
plot(NULL, xlim=c(0,1), ylim=c(0,lcat+2), xaxt="n", yaxt="n", xlab="",ylab="", bty="n")
for ( i in unique(posCat) ) {
    for (j in unique(posCat) ) {
        combpos_temp <- paste0(i,"-",j)
        thick_temp <- thickness_combpos[combpos_temp]
        if ( !is.na(thick_temp)) {
            # print(paste0(i,"-",j))
            lines(x=c(0,1), y=c(i,j), col=rgb(0.5,0.5,0.5,0.3), lwd=thick_temp)
        }
    }
}
axis(1, at = c(0,1), labels = c("Baltic","Fraser"), las=1, line=-1, cex.axis=2, tick = FALSE)
axis(side = 2, at=unique(posCat), labels=c("Fresh","Fresh/ marinetolerant","Brack/ freshtolerant","Brackish","Brack / marinetolerant","Marine/ freshtolerant","Marine"), las=2)
axis(side = 4, at=unique(posCat), labels=c("Fresh","Fresh/ marinetolerant","Brack/ freshtolerant","Brackish","Brack / marinetolerant","Marine/ freshtolerant","Marine"), las=2)
legend("top", legend=c(paste0("       ",maxCount," OTUs"),"",paste0("       ",minCount," OTU")), lwd=c(max(thickness_combpos),NA,min(thickness_combpos)), bty="n", col=rgb(0.5,0.5,0.5,0.3))
dev.off()

####### Shared otus by classification (simple) #####

combinedmb.simple <- matrix(ncol=6, nrow=length(commonOTUs),dimnames=list(commonOTUs,c("btype","ftype","bb1","bb2","fb1","fb2")))
for ( otu in commonOTUs ) {
    ftemp <- as.vector(mbF16[otu,c('typeSimple','boundaries','boundariestwo')])
    btemp <- as.vector(mbB16[otu,c('typeSimple','boundaries','boundariestwo')])
    combinedmb.simple[otu,] <- c(as.character(btemp[1]),as.character(ftemp[1]),as.numeric(btemp[c(2,3)]), as.numeric(ftemp[c(2,3)]))
}

categories.simple <- c("freshRestricted","brackishRestricted","marineRestricted")
lcat <- length(categories.simple)

fpos <- sapply(combinedmb.simple[,c("ftype")],function(x) {
    match(x,categories.simple)
})
bpos <- sapply(combinedmb.simple[,c("btype")],function(x) {
    match(x,categories.simple)
})

combpos <- as.data.frame(cbind(bpos,fpos,comb=paste0(bpos, "-",fpos)))
counts_combpos <- table(combpos[,'comb'])
maxCount <- max(counts_combpos)
minCount <- min(counts_combpos)
dividend <- (maxCount)/40
thickness_combpos <- counts_combpos/dividend

pdf("shared_otus_lineplot_simple.pdf")
par(mar=c(4.2,10.2,2.2,10.2))
plot(NULL, xlim=c(0,1), ylim=c(0,lcat+1), xaxt="n", yaxt="n", xlab="",ylab="")
for ( i in 1:lcat ) {
    for (j in 1:lcat) {
        combpos_temp <- paste0(i,"-",j)
        thick_temp <- thickness_combpos[combpos_temp]
        if (!is.na(thick_temp)) {
            lines(x=c(0,1), y=c(i,j), col=rgb(0.5,0.5,0.5,0.3), lwd=thick_temp)
        }
    }
}
axis(1, at = c(0,1), labels = c("Baltic","Fraser"), las=2)
axis(side = 2, at=seq(1,lcat), labels=categories.simple, las=2)
axis(side = 4, at=seq(1,lcat), labels=categories.simple, las=2)
dev.off()

########## shared otus by peaking location #########
blim <- c(0,36)
flim <- c(0,34)

categories <- c("freshbloom","freshRestricted","freshPeak","brackishPeakLoToler","brackishBloom","BrackishPeakAllToler","brackishRestricted","brackishPeakHiToler","marinePeak","marineRestricted","marineBloom")
lcat <- length(categories)

peaksTable <- matrix(ncol=2,nrow=length(commonOTUs), dimnames=list(commonOTUs,c("b","f")))
for ( i in commonOTUs)  {
    #BALTIC
    if ( combinedmb.simple[i,'btype'] == 'freshRestricted' ) {
        peakpos <- (as.numeric(combinedmb.simple[i,'bb1']) - flim[1])/2
    } else if ( combinedmb.simple[i,'btype'] == 'brackishRestricted' ) {
        peakpos <- (as.numeric(combinedmb.simple[i,'bb2']) + as.numeric(combinedmb.simple[i,'bb1']))/2
    } else if ( combinedmb.simple[i,'btype'] == 'marineRestricted' ) {
        peakpos <- (flim[2] - as.numeric(combinedmb.simple[i,'bb1']))/2
    }
    peaksTable[i,"b"] <- peakpos
    
    #FRASER
    if ( combinedmb.simple[i,'ftype'] == 'freshRestricted' ) {
        peakpos <- (as.numeric(combinedmb.simple[i,'fb1']) - flim[1])/2
    } else if ( combinedmb.simple[i,'ftype'] == 'brackishRestricted' ) {
        peakpos <- (as.numeric(combinedmb.simple[i,'fb2']) + as.numeric(combinedmb.simple[i,'fb1']))/2
    } else if ( combinedmb.simple[i,'ftype'] == 'marineRestricted' ) {
        peakpos <- (flim[2] - as.numeric(combinedmb.simple[i,'fb1']))/2
    }
    peaksTable[i,"f"] <- peakpos
}
# 
# 
# combpos <- as.data.frame(cbind(fpos,bpos, comb=paste0(fpos,"-",bpos)))
# counts_combpos <- table(combpos[,'comb'])
# maxCount <- max(counts_combpos)
# minCount <- min(counts_combpos)
# dividend <- (maxCount)/40
# thickness_combpos <- counts_combpos/dividend

pdf("shared_otus_lineplot_bysal.pdf")
par(mar=c(4.2,4.2,2.2,4.2))
plot(NULL, xlim=c(0,1), ylim=c(0,36), xaxt="n", yaxt="n", xlab="",ylab="")
for ( r in 1:nrow(peaksTable)) {
    lines(x=c(0,1),y=peaksTable[r,], col=rgb(0.5,0.5,0.5,0.2))
}
axis(1, at = c(0,1), labels = c("Baltic","Fraser"), las=2)
axis(side = 2, at=seq(0,35,by=5), labels=seq(0,35,by=5), las=1)
axis(side = 4, at=seq(0,35,by=5), labels=seq(0,35,by=5), las=1)
mtext(side=2, text="Salinity", line=3)
mtext(side=4, text="Salinity", line=3)
dev.off()

####### COMPOSITION OF SPECIALIST TYPES #########

taxaRefB16FP = '16sBALTIC/taxaIDLegend.txt'
taxaRefF16FP = '16sFRASER/taxaIDLegend.txt'
taxaRefF18FP = '18sFraser/taxaIDLegend.txt'

taxaRefB16 = read.delim(taxaRefB16FP,header=FALSE, row.names=1)
taxaRefF16 = read.delim(taxaRefF16FP,header=FALSE, row.names=1)
taxaRefF18 = read.delim(taxaRefF18FP,header=FALSE, row.names=1)

### CLASS LEVEL FOR BALTIC
B16.fresh.taxa = as.character(taxaRefB16[rownames(mbB16.filt[mbB16.filt$typeSimple == "freshRestricted",]),])
B16.marine.taxa = as.character(taxaRefB16[rownames(mbB16.filt[mbB16.filt$typeSimple == "marineRestricted",]),])
B16.brackish.taxa = as.character(taxaRefB16[rownames(mbB16.filt[mbB16.filt$typeSimple == "brackishRestricted",]),])

# Breakdown by class
# FRESH
B16.fresh.taxa.split = sapply(1:length(B16.fresh.taxa), function(x) {strsplit(B16.fresh.taxa[x], split="; ", fixed=TRUE)})
round(table(sapply(B16.fresh.taxa.split, FUN="[[", INDEX=3))/length(B16.fresh.taxa.split),2)
# MARINE
B16.marine.taxa.split = sapply(1:length(B16.marine.taxa), function(x) {strsplit(B16.marine.taxa[x], split="; ", fixed=TRUE)})
round(table(sapply(B16.marine.taxa.split, FUN="[[", INDEX=3))/length(B16.marine.taxa.split),2)
# BRACKISH
B16.brackish.taxa.split = sapply(1:length(B16.brackish.taxa), function(x) {strsplit(B16.brackish.taxa[x], split="; ", fixed=TRUE)})
round(table(sapply(B16.brackish.taxa.split, FUN="[[", INDEX=3))/length(B16.brackish.taxa.split),2)



### CLASS LEVEL FOR FRASER
F16.fresh.taxa = as.character(taxaRefF16[rownames(mbF16.filt[mbF16.filt$typeSimple == "freshRestricted",]),])
F16.marine.taxa = as.character(taxaRefF16[rownames(mbF16.filt[mbF16.filt$typeSimple == "marineRestricted",]),])
F16.brackish.taxa = as.character(taxaRefF16[rownames(mbF16.filt[mbF16.filt$typeSimple == "brackishRestricted",]),])

# Breakdown by class
# FRESH
F16.fresh.taxa.split = sapply(1:length(F16.fresh.taxa), function(x) {strsplit(F16.fresh.taxa[x], split="; ", fixed=TRUE)})
round(table(sapply(F16.fresh.taxa.split, FUN="[[", INDEX=3))/length(F16.fresh.taxa.split),2)
# MARINE
F16.marine.taxa.split = sapply(1:length(F16.marine.taxa), function(x) {strsplit(F16.marine.taxa[x], split="; ", fixed=TRUE)})
round(table(sapply(F16.marine.taxa.split, FUN="[[", INDEX=3))/length(F16.marine.taxa.split),2)
# BRACKISH
F16.brackish.taxa.split = sapply(1:length(F16.brackish.taxa), function(x) {strsplit(F16.brackish.taxa[x], split="; ", fixed=TRUE)})
round(table(sapply(F16.brackish.taxa.split, FUN="[[", INDEX=3))/length(F16.brackish.taxa.split),2)



### CLASS LEVEL FOR FRASER
F18.fresh.taxa = as.character(taxaRefF18[rownames(mbF18.filt[mbF18.filt$typeSimple == "freshRestricted",]),])
F18.marine.taxa = as.character(taxaRefF18[rownames(mbF18.filt[mbF18.filt$typeSimple == "marineRestricted",]),])
F18.brackish.taxa = as.character(taxaRefF18[rownames(mbF18.filt[mbF18.filt$typeSimple == "brackishRestricted",]),])

# Breakdown by class
# FRESH
F18.fresh.taxa.split = sapply(1:length(F18.fresh.taxa), function(x) {strsplit(F18.fresh.taxa[x], split="; ", fixed=TRUE)})
# add command because 18s has some "Unassigned"
F18.fresh.taxa.split <- lapply(F18.fresh.taxa.split, `length<-`, max(lengths(F18.fresh.taxa.split)))
round(table(sapply(F18.fresh.taxa.split, FUN="[[", INDEX=3))/length(F18.fresh.taxa.split),2)

# MARINE
F18.marine.taxa.split = sapply(1:length(F18.marine.taxa), function(x) {strsplit(F18.marine.taxa[x], split="; ", fixed=TRUE)})
F18.marine.taxa.split <- lapply(F18.marine.taxa.split, `length<-`, max(lengths(F18.marine.taxa.split)))
round(table(sapply(F18.marine.taxa.split, FUN="[[", INDEX=3))/length(F18.marine.taxa.split),2)

# BRACKISH
F18.brackish.taxa.split = sapply(1:length(F18.brackish.taxa), function(x) {strsplit(F18.brackish.taxa[x], split="; ", fixed=TRUE)})
F18.brackish.taxa.split <- lapply(F18.brackish.taxa.split, `length<-`, max(lengths(F18.brackish.taxa.split)))
round(table(sapply(F18.brackish.taxa.split, FUN="[[", INDEX=3))/length(F18.brackish.taxa.split),2)





##### FIND OTUS that are classified as opposites #######
# Number of classifications that are EXACT same
sum(combinedmb[,'btype'] == combinedmb[,'ftype'])
nrow(combinedmb)

# Number of classifications that are EXACT same; noclass ommitted
combinedmb.filt <- combinedmb[-grep("noclass", combinedmb),]
sum(combinedmb.filt[,'btype'] == combinedmb.filt[,'ftype'])
nrow(combinedmb.filt)
65/214

# Number of classifications that are same broad type
sum(combinedmb.simple[,'btype'] == combinedmb.simple[,'ftype'])
nrow(combinedmb.simple)
128/303

# Number of classifications that are same broad type; noclass ommitted
combinedmb.simple.filt <- combinedmb.simple[-grep("noclass", combinedmb.simple),]
sum(combinedmb.simple.filt[,'btype'] == combinedmb.simple.filt[,'ftype'])
nrow(combinedmb.simple.filt)
99/214

# find the one that is opposite
oppositeClass <- c()
oppositeClass <- c(oppositeClass, names(which((combinedmb[,'btype'] == 'brackishPeakLoToler') & (combinedmb[,'ftype'] == 'marineRestricted'))))
names(which((combinedmb[,'btype'] == 'brackishPeakHiToler') & (combinedmb[,'ftype'] == 'brackishPeakHiToler')))


names(which((combinedmb[,'btype'] == 'brackishRestricted') & (combinedmb[,'ftype'] == 'brackishRestricted')))


# 
# ##### Comparing shared OTUs WITH HIABUND ONLY ######
# 
# ####### Shared otus by classification (complex) #####
# commonOTUs <- rownames(mbB16.hiabund[na.omit(match(rownames(mbF16.hiabund), rownames(mbB16.hiabund))),])
# 
# combinedmb <- matrix(ncol=6, nrow=length(commonOTUs),dimnames=list(commonOTUs,c("btype","ftype","bb1","bb2","fb1","fb2")))
# for ( otu in commonOTUs ) {
#     ftemp <- as.vector(mbF16.hiabund[otu,c('type','boundaries','boundariestwo')])
#     btemp <- as.vector(mbB16.hiabund[otu,c('type','boundaries','boundariestwo')])
#     combinedmb[otu,] <- c(as.character(btemp[1]),as.character(ftemp[1]),as.numeric(btemp[c(2,3)]), as.numeric(ftemp[c(2,3)]))
# }
# 
# categories <- c("freshBloom" # 1
#                 ,"freshRestricted" #1
#                 ,"freshPeak" #2
#                 ,"brackishPeakLoToler" #4
#                 ,"brackishBloom" #5
#                 ,"brackishRestricted" #5
#                 ,"brackishPeakHiToler" #6
#                 ,"marinePeak" #8
#                 ,"marineRestricted" #9
#                 ,"marineBloom") #9
# posCat <- c(1,1,2,4,5,5,6,8,9,9)
# 
# # lcat <- length(categories)
# lcat <- max(posCat)
# 
# fpos <- sapply(combinedmb[,c("ftype")],function(x) {
#     posCat[match(x,categories)]
# })
# bpos <- sapply(combinedmb[,c("btype")],function(x) {
#     posCat[match(x,categories)]
# })
# 
# combpos <- as.data.frame(cbind(bpos,fpos, comb=paste0(bpos,"-",fpos)))
# # combpos[which(combpos$comb == '9-4'),]
# counts_combpos <- table(combpos[,'comb'])
# maxCount <- max(counts_combpos)
# minCount <- min(counts_combpos)
# dividend <- (maxCount)/40
# thickness_combpos <- counts_combpos/dividend
# 
# pdf("shared_otus_lineplot_complex_hiabund.pdf")
# # quartz()
# par(mar=c(4.2,10.2,2.2,10.2))
# plot(NULL, xlim=c(0,1), ylim=c(0,lcat+2), xaxt="n", yaxt="n", xlab="",ylab="", bty="n")
# for ( i in unique(posCat) ) {
#     for (j in unique(posCat) ) {
#         combpos_temp <- paste0(i,"-",j)
#         thick_temp <- thickness_combpos[combpos_temp]
#         if ( !is.na(thick_temp)) {
#             # print(paste0(i,"-",j))
#             lines(x=c(0,1), y=c(i,j), col=rgb(0.5,0.5,0.5,0.3), lwd=thick_temp)
#         }
#     }
# }
# axis(1, at = c(0,1), labels = c("Baltic","Fraser"), las=1, line=-1, cex.axis=2, tick = FALSE)
# axis(side = 2, at=unique(posCat), labels=c("Fresh","Fresh/ marinetolerant","Brack/ freshtolerant","Brackish","Brack / marinetolerant","Marine/ freshtolerant","Marine"), las=2)
# axis(side = 4, at=unique(posCat), labels=c("Fresh","Fresh/ marinetolerant","Brack/ freshtolerant","Brackish","Brack / marinetolerant","Marine/ freshtolerant","Marine"), las=2)
# legend("top", legend=c(paste0("       ",maxCount," OTUs"),"",paste0("       ",minCount," OTU")), lwd=c(max(thickness_combpos),NA,min(thickness_combpos)), bty="n", col=rgb(0.5,0.5,0.5,0.3))
# dev.off()

