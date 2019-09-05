#!bin/bash Rscript

############# GOAL ##############
# This script takes modelBoundaries_type.txt and the OTU Table (txt) and plots tolerance ranges
# for individual OTUs, grouped by phylogenetic group.
# Working directory should be SALBIN

###### Packages #######
library("tidyverse")
library("ggplot2")
library("gridExtra")
library("gtable")

############### Set pathways ##################

# Baltic 16S
B16_mb_pwd <- "./16sBaltic/modelBoundaries_type.txt"
B16_otu_pwd <- "./16sBaltic/OTUTableText.txt"
B16_abund_pwd <- "./16sBaltic/taxa_abundances_across_gradient.txt"

# Fraser 16S
F16_mb_pwd <- "./16sFraser/modelBoundaries_type.txt"
F16_otu_pwd <- "./16sFraser/OTUTableText.txt"
F16_abund_pwd <- "./16sFraser/taxa_abundances_across_gradient.txt"

# Fraser 18S
F18_mb_pwd <- "./18sFraser/modelBoundaries_type.txt"
F18_otu_pwd <- "./18sFraser/OTUTableText.txt"
F18_abund_pwd <- "./18sFraser/taxa_abundances_across_gradient.txt"

############### Load files ##################

# Baltic 16S
B16_mb <- read.delim(B16_mb_pwd, header=TRUE, stringsAsFactors=FALSE)
B16_otu <- read.delim(B16_otu_pwd, header=TRUE)
B16_abund <- read.delim(B16_abund_pwd, header=FALSE, stringsAsFactors = FALSE, row.names = 1)

# Fraser 16S
F16_mb <- read.delim(F16_mb_pwd, header=TRUE, stringsAsFactors=FALSE)
F16_otu <- read.delim(F16_otu_pwd, header=TRUE)
F16_abund <- read.delim(F16_abund_pwd, header=FALSE, stringsAsFactors = FALSE, row.names = 1)

# Fraser 18S
F18_mb <- read.delim(F18_mb_pwd, header=TRUE, stringsAsFactors=FALSE)
F18_otu <- read.delim(F18_otu_pwd, header=TRUE)
F18_abund <- read.delim(F18_abund_pwd, header=FALSE, stringsAsFactors = FALSE, row.names = 1)


############# changing otu tables ###########
# Baltic 16
B16_otu_tax <- B16_otu[,c(1,ncol(B16_otu))]
B16_otu_tax$taxonomy <- gsub("D_[0-9]*__","",B16_otu_tax$taxonomy)
# B16_otu_main <- B16_otu[,2:(ncol(B16_otu)-1)]
# rownames(B16_otu_main) <- B16_otu[,1]

# Fraser 16
F16_otu_tax <- F16_otu[,c(1,ncol(F16_otu))]
F16_otu_tax$taxonomy <- gsub("D_[0-9]*__","",F16_otu_tax$taxonomy)
# F16_otu_main <- F16_otu[,2:(ncol(F16_otu)-1)]
# rownames(F16_otu_main) <- F16_otu[,1]

# Fraser 18
F18_otu_tax <- F18_otu[,c(1,ncol(F18_otu))]
F18_otu_tax$taxonomy <- gsub("__","",F18_otu_tax$taxonomy)
# F18_otu_main <- F18_otu[,2:(ncol(F18_otu)-1)]
# rownames(F18_otu_main) <- F18_otu[,1]

##### Make a separated taxonomy #########
B16_otu_tax <- B16_otu_tax %>%
    as_tibble() %>%
    separate(col=taxonomy, into=c("D0","D1","D2","D3","D4","D5","D6"), sep = "; ", remove=FALSE,fill="left")
F16_otu_tax <- F16_otu_tax %>%
    as_tibble() %>%
    separate(col=taxonomy, into=c("D0","D1","D2","D3","D4","D5","D6"), sep = "; ", remove=FALSE,fill="left")
F18_otu_tax <- F18_otu_tax %>%
    as_tibble() %>%
    separate(col=taxonomy, into=c("D0","D1","D2","D3","D4","D5","D6","D7"), sep = "; ", remove=FALSE,fill="left")
#### Adjust abundance table #####
# colnames(B16_abund) <- gsub("X","",colnames(B16_abund))
# colnames(F16_abund) <- gsub("X","",colnames(F16_abund))
# colnames(F18_abund) <- gsub("X","",colnames(F18_abund))

B16_abund_t <- as.data.frame(t(B16_abund))
F16_abund_t <- as.data.frame(t(F16_abund))
F18_abund_t <- as.data.frame(t(F18_abund))


###### Group by family, genus ########

dir.create("Tolerance_plots_by_family")
for ( d in c("B16","F16","F18") ) {
    dir.create(paste0("Tolerance_plots_by_family/",d))
    if ( d %in% c("B16","F16")) {
        lvl="D4"
    } else {
        lvl="D5"
    }
    # Get rid of "no class" individuals
    keep_taxa <- get(paste0(d,"_mb")) %>%
        as_tibble() %>%
        filter(typeSimple!="noclass") %>%
        pull(taxa)
    temp_otu_tax <- get(paste0(d,"_otu_tax")) %>%
        filter(X.OTUID %in% keep_taxa)
    
    # Get list of families
    f_all <- unique(unlist(temp_otu_tax[,lvl]))
    
    for ( f in f_all) {
        current.fam <- unlist(temp_otu_tax[temp_otu_tax[,lvl] == f,"X.OTUID"])
        nas <- which(is.na(current.fam))
        if (length(nas)>0) {
            current.fam <- current.fam[-nas]
        }
        current.abund <- get(paste0(d,"_abund_t"))[,c(as.character(current.fam), "Gradient")]
        
        maxSal <- max(get(paste0(d,"_abund_t"))[,"Gradient"])
        
        gg1<- current.abund %>%
            as_tibble() %>%
            gather(key=OTU, value=RelAbund, -Gradient) %>%
            ggplot() +
            geom_segment(aes(x=Gradient,xend=Gradient, y=0, yend=RelAbund, col=OTU), alpha=0.2, lwd=2) +
            xlim(0,maxSal) +
            scale_y_continuous(position = "right") +
            ylab("Relative Abundance") +
            xlab("Salinity")
        
        gg2 <- get(paste0(d,"_mb")) %>%
            as_tibble() %>%
            filter(typeSimple!="noclass") %>%
            filter(taxa %in% current.fam) %>%
            select(taxa,typeSimple, boundaries, boundariestwo) %>%
            mutate(min=ifelse(typeSimple=="marineRestricted", boundaries, ifelse(typeSimple=="brackishRestricted",boundaries, 0))
                   , max=ifelse(typeSimple=="freshRestricted",boundaries, ifelse(typeSimple=="brackishRestricted",boundariestwo, maxSal)))%>%
            ggplot() +
            geom_segment(aes(x=min,xend=max, y=taxa,yend=taxa,col=as.factor(taxa))) +
            scale_y_discrete(position="right")+
            xlab("Salinity") +
            xlim(0,maxSal)
        
        g1 <- ggplotGrob(gg1)
        g2 <- ggplotGrob(gg2)
        g <- rbind(g1, g2, size = "first")
        g$widths <- grid::unit.pmax(g1$widths, g2$widths)
        grid::grid.newpage()
        pdf(paste0("Tolerance_plots_by_family/",d,"/",f,".pdf"))
        grid::grid.draw(g)
        dev.off()
    }
    
}


