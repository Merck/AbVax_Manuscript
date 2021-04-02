#!/usr/bin/env Rscript
# Copyright © 2021 Merck Sharp & Dohme Corp. a subsidiary of Merck & Co., Inc., Kenilworth, NJ, USA.
# Geoffrey Hannigan edited by Begüm Topçuoglu
# Systems Biology Research Group
# Merck Exploratory Science Center, Cambridge


#################
# Set Libraries #
#################

library("vegan")
library("stringr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("viridis")
library("cowplot")
library("gtools")
library("parallel")
library("data.table")
library("ggrepel")


###############
# DATA IMPORT #
###############

# Import the master file for correlation analysis
master_taxa <- as.data.frame(fread(
	file = "./data/master_taxa.tsv",
	header=TRUE,
	sep="\t"))

# Import the master metadata files for study set
master_meta <- as.data.frame(fread(
	file = "./data/master_metadata.tsv",
	header=TRUE,
	sep="\t"))%>%
	filter(StudyID %in% c("MOUSE_002", "MOUSE_004")) %>%
	filter(PostDose_2 != 0)%>%
	filter(Antibiotic != "Water")

########################
# CORRELATION ANALYSIS #
########################

# Can bring together the two tables as well
mpmerge <- inner_join(master_taxa, master_meta, by = c("SampleID" = "SampleID"))%>%
	filter(StudyID == "MOUSE_002")

mpmerge$TaxaID <- sub('_.*', '', mpmerge$TaxaID)

# GENUS
corlist <- lapply(unique(mpmerge$StudyID), function (i) {
		tdf <- mpmerge %>%
			group_by(TaxaID) %>%
			summarize(
				corval = cor.test(relabund, PostDose_2, method = "spearman", adjust = "BH")$estimate,
				pvalue = cor.test(relabund, PostDose_2, method = "spearman", adjust = "BH")$p.value
			) %>%
			drop_na() %>%
			as.data.frame()
		tdf$StudyID <- i
		return(tdf)
	})

corlisto <- do.call(rbind, corlist) %>%
	mutate(level = case_when(pvalue < 0.01 ~ "FDR-corrected significance", pvalue >= 0.01 ~ "Not significant"))

corlisto$TaxaID[which(corlisto$level == "Not significant")] <- NA




volcanop <- ggplot(corlisto, aes(x = corval, y = -log10(pvalue), color=level, label=TaxaID)) +
	theme_classic() +
	geom_point() +
	scale_colour_manual(values = c("red", "black")) +
	geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
	geom_vline(xintercept = 0.5, linetype = "dashed") +
	geom_vline(xintercept = -0.5, linetype = "dashed") +
	xlim(c(-1,1)) +
	xlab("Correlation Coefficient") +
	ylab("FDR Corrected P-Value (-log10)") +
	theme(legend.position="none",
				axis.title.x = element_text(size = 24),
  			axis.title.y = element_text(size = 24),
				axis.text = element_text(size = 24, color = "black"),
				axis.ticks.length = unit(0.5, "cm")) +
	geom_label_repel(size = 7)




pdf(file = "./figures/genus_002_volcano.pdf", width = 6, height = 6)

	volcanop

dev.off()







mpmerge <- inner_join(master_taxa, master_meta, by = c("SampleID" = "SampleID"))%>%
 filter(StudyID == "MOUSE_004")

mpmerge$TaxaID <- sub('_.*', '', mpmerge$TaxaID)

# GENUS
corlist <- lapply(unique(mpmerge$StudyID), function (i) {
	 tdf <- mpmerge %>%
		 group_by(TaxaID) %>%
		 summarize(
			 corval = cor.test(relabund, PostDose_2, method = "spearman", adjust = "BH")$estimate,
			 pvalue = cor.test(relabund, PostDose_2, method = "spearman", adjust = "BH")$p.value
		 ) %>%
		 drop_na() %>%
		 as.data.frame()
	 tdf$StudyID <- i
	 return(tdf)
 })

corlisto <- do.call(rbind, corlist) %>%
	mutate(level = case_when(pvalue < 0.01 ~ "FDR-corrected significance", pvalue >= 0.01 ~ "Not significant"))

corlisto$TaxaID[which(corlisto$level == "Not significant")] <- NA



	volcanop <- ggplot(corlisto, aes(x = corval, y = -log10(pvalue), color=level, label=TaxaID)) +
	theme_classic() +
	geom_point() +
	scale_colour_manual(values = c("red", "black")) +
	geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
	geom_vline(xintercept = 0.5, linetype = "dashed") +
	geom_vline(xintercept = -0.5, linetype = "dashed") +
	xlim(c(-1,1)) +
	xlab("Correlation Coefficient") +
	ylab("FDR Corrected P-Value (-log10)") +
	theme(legend.position="none",
				axis.title.x = element_text(size = 24),
  			axis.title.y = element_text(size = 24),
				axis.text = element_text(size = 24, color = "black"),
				axis.ticks.length = unit(0.5, "cm")) +
	geom_label_repel(size = 7)



pdf(file = "./figures/genus_004_volcano.pdf", width = 6, height = 6)

	volcanop

dev.off()
