#!/usr/bin/env Rscript
# Copyright © 2021 Merck Sharp & Dohme Corp. a subsidiary of Merck & Co., Inc., Kenilworth, NJ, USA.
# Geoffrey Hannigan
# Systems Biology Research Group
# Merck Exploratory Science Center, Cambridge

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(tibble)
library(gage)
library(cowplot)
library(gridExtra)
library(caret)


######################
# IMPORT DATA TABLES #
######################

# Import the taxonomic abundance table
iko <- read.delim(
	file = "./mAbx-Vax-002/data/ko/NormRelAbundMatrix.tsv",
	header = TRUE,
	sep = "\t")

# Import the metadata tables
imet <- read.delim(
	file = "./mAbx-Vax-002/data/metadata/abx002Metadata.tsv",
	header = TRUE,
	sep = "\t")

# Import the metadata tables
antit <- read.delim(
	file = "./mAbx-Vax-002/data/metadata/TiterTable.tsv",
	header = TRUE,
	sep = "\t")

antit$Antibiotic <- gsub("_Enrofloxacin", "", antit$Antibiotic)
antit <- antit %>%
	filter(Target != "Fungi") %>%
	filter(Antibiotic != "Ampicillin")


################
# GAGE Analyis #
################

# Here we want to run the analysis assuming differences in the same directions
# Also look only at the non-redundant pathways with significance


gageall <- lapply(as.vector(unique(antit$Antibiotic)[2:7]), function(i) {
	write(i, stdout())

	# Prep the pathway data
	kgko <- kegg.gsets("ko", id.type = "kegg")
	kgkogs <- kgko$kg.sets[ kgko$sigmet.idx ]

	waterdf <- iko %>%
		dplyr::select(contains("Watrer"))

	vancodf <- iko %>%
		dplyr::select(contains(i))

	ikoid <- iko %>%
		dplyr::select("KO")

	mcontdf <- cbind(ikoid, waterdf, vancodf)

	abxlength <- ncol(vancodf)

	# Add kegg ids
	row.names(mcontdf) <- mcontdf$KO
	mcontdf <- mcontdf[,-1]

	res <- gage(mcontdf,
		gsets= kgkogs,
		ref= c(1:5),
		samp= c(6:10),
		compare= "unpaired",
		same.dir= TRUE)

	ress <- as.data.frame(res$stat)

	ress$KOID <- row.names(ress)
	row.names(ress) <- NULL

	rko <- ress %>%
		gather(key = "SampleID", value = "Enrichment", -KOID) %>%
		filter(SampleID != "stat.mean") %>%
		filter(is.finite(Enrichment))

	# Cut down the sample names to match metadata
	rko$SampleID <- gsub("\\.1\\.[A-Za-z].+", "", rko$SampleID, perl = TRUE)
	rko$SampleID <- gsub("\\.", "_", rko$SampleID, perl = TRUE)

	irko <- left_join(rko, antit, by = c("SampleID" = "SampleID"))

	return(irko)
})

gagedf <- do.call(rbind, gageall)

gagedf <- gagedf[complete.cases(gagedf), ]


gagetable <- gagedf %>%
	group_by(KOID) %>%
	summarize(correlation = cor.test(Enrichment, OVA_IgG_titer_ng_ml)$estimate, fdr_pvalue = p.adjust(cor.test(Enrichment, OVA_IgG_titer_ng_ml)$p.value, method = "fdr", n = length(Enrichment))) %>%
	arrange(correlation) %>%
	filter(fdr_pvalue < 0.005) %>%
	as.data.frame()

#########################
#abby modify appearances#
#########################

sigplotlist2 <- lapply(c(gagetable$KOID), function(i) {
  write(i, stdout())

  oplot <- gagedf %>%
    filter(KOID == i) %>%
    ggplot(aes(x = Enrichment, y = OVA_IgG_titer_ng_ml)) +
    theme_classic() +
    geom_smooth(method = "lm", size=0.2, colour="red") +
    geom_point(shape=21, size=1.5, aes(fill = Antibiotic), stroke=.3) +
    #geom_point(shape=21, size=5, aes(fill = Antibiotic), stroke=.3) +
    theme(axis.text = element_text(family="sans", size=8),
          text = element_text(size = 8))+
          #legend.position="none")
    theme(legend.position="none")+
    scale_fill_manual(values=c("#FF7F00","#596A98","#FFE528","#449B75","#AC5782","#E41A1C")) +
    xlab("Metagenomic Pathway Enrichment") +
    ylab(expression(paste("anti-OVA IgG titer (", mu, "g/ml)")))

  return(oplot)
})
sigplotlist2[[1]]

pdf(
  file = "./mAbx-Vax-002/figures/titer_corr_ubiq_manuscript.pdf",
  width = 3,
  height = 1.7)
sigplotlist2[[1]]
dev.off()

pdf(
  file = "./mAbx-Vax-002/figures/titer_corr_sulf_manuscript.pdf",
  width = 3,
  height = 1.7)
sigplotlist2[[2]]
dev.off()

pdf(
  file = "./mAbx-Vax-002/figures/titer_corr_tryp_manuscript.pdf",
  width = 3,
  height = 1.7)
sigplotlist2[[3]]
dev.off()

pdf(
  file = "./mAbx-Vax-002/figures/titer_corr_sel_manuscript.pdf",
  width = 3,
  height = 1.7)
sigplotlist2[[4]]
dev.off()

##########################
#end abby's modifications#
#########################
