#!/usr/bin/env Rscript
# Copyright © 2021 Merck Sharp & Dohme Corp. a subsidiary of Merck & Co., Inc., Kenilworth, NJ, USA.
# Geoffrey Hannigan
# Systems Biology Research Group
# Merck Exploratory Science Center, Cambridge

#################
# Set Libraries #
#################

library("vegan")
library("dplyr")
library("tidyr")
library("ggplot2")
library("cowplot")
library("viridis")
library("stringr")
library("data.table")

######################
# Import Data Tables #
######################

# Import the OTU abundance table
# Pull in the large table quickly with fread
otucounts <- as.data.frame(fread(
	file="./mAbx-Vax-004/data/16S/renamed/OtuCounts.tsv",
	header=TRUE,
	sep="\t"))

# Import the OTU taxonomic assignments (down to genus level)
otutaxonomy <- read.delim(
	file="./mAbx-Vax-004/data/16S/renamed/OtuTaxonomy.tsv",
	header=TRUE,
	sep="\t")
# Format the IDs by parsing the ID column
otutaxonomy$Taxonomy <- as.character(gsub("\\(\\d+\\)", "", otutaxonomy$Taxonomy, perl = TRUE))

splittax <- as.data.frame(gsub(";", "", str_split_fixed(otutaxonomy$Taxonomy, ";", 6)))
colnames(splittax) <- c(
	"Kingdom",
	"Phylum",
	"Class",
	"Order",
	"Family",
	"Genus")

oparse <- cbind(otutaxonomy[,-3], splittax)

metadata <- read.delim(
	file="./mAbx-Vax-004/data/metadata/004Metadata.tsv",
	header=TRUE,
	sep="\t")

animalData <- read.delim(
	file = "./mAbx-Vax-004//data/metadata/004Animals.tsv",
	header = TRUE,
	sep = "\t")

metadata <- left_join(metadata, animalData, by = c("AnimalID" = "AnimalID"))

# Make a column that matched the taxonomy IDs
metadata$DiversigenSampleID_Fixed <- gsub("\\.", "_", metadata$DiversigenSampleID)



################
# Run Analysis #
################

# ALPHA

# Create a rarefaction curve to visualize depth
formotu <- otucounts
row.names(formotu) <- otucounts[,2]
formotu <- formotu[,-c(1:3)]

# There is just way too much excess OTU information here
# Cut those out with very few counts overall
formotu <- formotu[,colSums(formotu) > 10]

rsumotu <- as.data.frame(rowSums(formotu))
colnames(rsumotu) <- "counts"
rsumotu$sampleid <- rownames(rsumotu)
rownames(rsumotu) <- NULL

countplot <- ggplot(rsumotu, aes(x = sampleid, y = counts)) +
	theme_classic() +
	geom_bar(stat = "identity", width = 1) +
	coord_flip() +
	geom_hline(yintercept = 10000, linetype = "dashed") +
	theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# Get richness and Shannon diversity
# Make sure to rarefy dataset
rotu <- rrarefy(formotu, sample = 10000)
length(rownames(rotu))
rotu <- rotu[c(rowSums(rotu) >= 10000),]
rotu[ rotu < 2 ] <- 0
length(rownames(rotu))

# Correlate the resulting titers with the distances from water like in 002

sdist <- as.data.frame(as.matrix(vegdist(rotu, method = "bray")))
sdist$samplecompare <- row.names(sdist)
row.names(sdist) <- NULL
sdistg <- gather(sdist, sampleid, dist, -samplecompare)
sdistg <- left_join(sdistg, metadata, by = c("samplecompare" = "DiversigenSampleID_Fixed"))
sdistg <- left_join(sdistg, metadata, by = c("sampleid" = "DiversigenSampleID_Fixed"))

sf <- sdistg %>%
	filter(Antibiotics.y == "Water") %>%
	filter(TreatmentDays.y == "Water") %>%
	filter(TreatmentDays.x == 7 | TreatmentDays.x == "Water") %>%
	filter(Day.x == -1 & Day.y == -1) %>%
	filter(dist != 0) %>%
	select(samplecompare, sampleid, dist, Antibiotics.x, Antibiotics.y, WashoutDays.x, WashoutDays.y, Day.x, Day.y, OVA.IgG.2.x)

sf <- sf %>%
	group_by(samplecompare) %>%
	summarize(mean = mean(dist)) %>%
	as.data.frame()

sf <- left_join(sf, metadata, by = c("samplecompare" = "DiversigenSampleID_Fixed")) %>%
	select(samplecompare, mean, Antibiotics, WashoutDays, Day, OVA.IgG.2)

beta_corr <- ggplot(sf, aes(x = mean, y = OVA.IgG.2)) +
	theme_classic() +
	geom_point(aes(colour = WashoutDays)) +
	xlab("Dist From Water") +
	ylab("Vaccine Titer") +
	geom_smooth(method = "lm")

pdf(
	file = "./mAbx-Vax-004/figures/Beta_Washout_Correlation.pdf",
	width = 6,
	height = 4)

	beta_corr

dev.off()
