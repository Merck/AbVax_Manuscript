#!/usr/bin/env Rscript
# Copyright © 2021 Merck Sharp & Dohme Corp. a subsidiary of Merck & Co., Inc., Kenilworth, NJ, USA.
# Geoffrey Hannigan
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
library("RColorBrewer")
library("Hmisc")
library("cowplot")



######################
# Import Data Tables #
######################

# Import the OTU abundance table
otucounts <- read.delim(
	file="./mAbx-Vax-002/data/16S/renamed/OtuCounts.tsv",
	header=TRUE,
	sep="\t")

# Import the OTU taxonomic assignments (down to genus level)
otutaxonomy <- read.delim(
	file="./mAbx-Vax-002/data/16S/renamed/OtuTaxonomy.tsv",
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

antit$SampleID <- gsub("c_10", "c_10_Minus", antit$SampleID)

imet <- left_join(imet, antit, by = c("SampleID" = "SampleID"))


############
# TAXONOMY #
############

# Create a rarefaction curve to visualize depth
formotu <- otucounts
row.names(formotu) <- otucounts[,2]
formotu <- formotu[,-c(1:3)]

# Sequence counts across the board
seqcounts <- as.data.frame(rowSums(formotu))
seqcounts$ID <- rownames(seqcounts)
colnames(seqcounts) <- c("Count", "SampleID")

# Rarefy the samples to even depth
rotu <- rrarefy(formotu, sample = 10000)

countdf <- as.data.frame(rowSums(rotu))
colnames(countdf) <- "Counts"
countdf$SampleID <- rownames(countdf)
filternames <- countdf %>%
	filter(Counts >= 10000)

rotu <- subset(rotu, rownames(rotu) %in% filternames$SampleID)


# Format for tidy analysis
gcount <- as.data.frame(rotu)
gcount$SampleID <- rownames(gcount)
gcount <- gcount %>%
	gather(value = Count, key = OTU, -SampleID) %>%
	as.data.frame()


# Add in the taxonomic IDs
tcount <- left_join(gcount, oparse, by = c("OTU" = "OTU"))

# Extract abundance at the ORDER level
gpec <- tcount %>%
	dplyr::select("SampleID", "Order", "Count") %>%
	group_by(SampleID, Order) %>%
	dplyr::summarize(scount = sum(Count)) %>%
	as.data.frame()

# Get the total sequence count per sample
tcount <- gpec %>%
	group_by(SampleID) %>%
	dplyr::summarize(total = sum(scount)) %>%
	as.data.frame()

# Merge in the total counts
mcount <- left_join(gpec, tcount, by = c("SampleID" = "SampleID"))
mcount$relabund <- 100 * mcount$scount / mcount$total

# Confirm that they add to 100
mcount %>%
	group_by(SampleID) %>%
	dplyr::summarize(totalpercent = sum(relabund)) %>%
	as.data.frame()

# Add in the metadata
mcount <- left_join(mcount, imet, by = c("SampleID" = "SampleID")) %>%
	na.omit()

# Filter out the low counts based on mock data accuracy
mcount <- mcount %>%
	filter(scount > 10)

# Order the antibiotics using the of richness ordering below (sorry but this only works
# after running the code below)

mcount$Antibiotic.x <- factor(mcount$Antibiotic.x, levels = of$Antibiotic.x)
mcount <- na.omit(mcount)


##################
# GENUS TAXONOMY #
##################



# Rarefy the samples to even depth
rotu <- rrarefy(formotu, sample = 10000)

countdf <- as.data.frame(rowSums(rotu))
colnames(countdf) <- "Counts"
countdf$SampleID <- rownames(countdf)
filternames <- countdf %>%
	filter(Counts >= 10000)

rotu <- subset(rotu, rownames(rotu) %in% filternames$SampleID)


# Format for tidy analysis
gcount <- as.data.frame(rotu)
gcount$SampleID <- rownames(gcount)
gcount <- gcount %>%
	gather(value = Count, key = OTU, -SampleID) %>%
	as.data.frame()


# Add in the taxonomic IDs
tcount <- left_join(gcount, oparse, by = c("OTU" = "OTU"))

# Extract abundance at the Genus level
gpec <- tcount %>%
	dplyr::select("SampleID", "Genus", "Count") %>%
	group_by(SampleID, Genus) %>%
	dplyr::summarize(scount = sum(Count)) %>%
	as.data.frame()

# Get the total sequence count per sample
tcount <- gpec %>%
	group_by(SampleID) %>%
	dplyr::summarize(total = sum(scount)) %>%
	as.data.frame()

# Merge in the total counts
mcount <- left_join(gpec, tcount, by = c("SampleID" = "SampleID"))
mcount$relabund <- 100 * mcount$scount / mcount$total

# Get the genus taxa with highest average abundance per sample
top_mean_df <- mcount %>%
	group_by(Genus) %>%
	summarise(mean_abund = mean(relabund)) %>%
	arrange(desc(mean_abund)) %>%
	as.data.frame() %>%
	head(n = 10)

mcother <- mcount %>%
	filter(!Genus %in% top_mean_df$Genus) %>%
	group_by(SampleID) %>%
	summarise(relabund = sum(relabund)) %>%
	as.data.frame()

mcother$Genus <- "Other"
mcother <- select(mcother, SampleID, Genus, relabund)

mcountf <- mcount %>%
	filter(Genus %in% top_mean_df$Genus) %>%
	select(SampleID, Genus, relabund)

# Add in the other group to top 10
mcountm <- rbind(mcountf, mcother)

# Confirm that they add to 100
mcountm %>%
	group_by(SampleID) %>%
	dplyr::summarize(totalpercent = sum(relabund)) %>%
	as.data.frame()

# Add in the metadata
mcountm <- left_join(mcountm, imet, by = c("SampleID" = "SampleID")) %>%
	na.omit()

mcountm <- na.omit(mcountm)

# Order by the top taxa averages
taxav <- c(as.character(top_mean_df$Genus), "Other")

mcountm$Genus <- factor(mcountm$Genus, levels=rev(taxav))

getPalette = colorRampPalette(brewer.pal(9, "Set1"))

mcountm$Antibiotic.x <- gsub("Baytril_Enrofloxacin", "Baytril", mcountm$Antibiotic.x)





#############
# DIVERSITY #
#############

# ALPHA
rotu <- rrarefy(formotu, sample = 10000)
rotu[ rotu < 2 ] <- 0
dotu <- as.data.frame(specnumber(rotu))
colnames(dotu) <- "Richness"


dotu$Shannon <- diversity(rotu, index = "shannon")
dotu$sampleid <- row.names(dotu)

# Fix up the metadata file names
imetf <- imet
imetf$SampleID <- gsub("MSA.2002", "MSA_2002", imet$SampleID)
imetf$SampleID <- gsub("MSA.1003", "Mock", imet$SampleID)

# Merge
# Remove na
# Order by mean
mdotu <- left_join(dotu, imetf, by = c("sampleid" = "SampleID")) %>%
	na.omit()

of <- mdotu %>%
	group_by(Antibiotic.x) %>%
	dplyr::summarize(mean = mean(Richness)) %>%
	arrange(desc(mean)) %>%
	as.data.frame()

mdotu$Antibiotic.x <- factor(mdotu$Antibiotic.x, levels = of$Antibiotic.x)

smdotu <- left_join(dotu, imetf, by = c("sampleid" = "SampleID")) %>%
	na.omit()

sof <- smdotu %>%
	group_by(Antibiotic.x) %>%
	dplyr::summarize(mean = mean(Shannon)) %>%
	arrange(desc(mean)) %>%
	as.data.frame()

smdotu$Antibiotic.x <- factor(smdotu$Antibiotic.x, levels = sof$Antibiotic.x)


# Without Mocks

mdotu <- left_join(dotu, imetf, by = c("sampleid" = "SampleID")) %>%
	na.omit()


odist <- vegdist(rotu, method = "bray")

ORD_NMDS <- metaMDS(comm = odist, k=2, trymax = 100)
ORD_FIT = data.frame(MDS1 = ORD_NMDS$points[,1], MDS2 = ORD_NMDS$points[,2])
ORD_FIT$SampleID <- rownames(ORD_FIT)
# Get metadata
ORD_FIT <- left_join(ORD_FIT, imet, by = c("SampleID" = "SampleID"))
# Remove mocks from this part, they will distort ordination
ORD_FIT <- filter(ORD_FIT, SampleID != "Mock") %>%
	filter(SampleID != "MSA_2002") %>%
	filter(Target.x != "Fungi") %>%
	as.data.frame()

ORD_FIT$Antibiotic.x <- gsub("Baytril_Enrofloxacin", "Baytril", ORD_FIT$Antibiotic.x)

theorder <- c("Water", "Baytril", "Tobramycin", "Metronidazole", "Ampicillin", "Neomycin", "Clindamycin", "Vancomycin")

ORD_FIT$Antibiotic.x <- factor(ORD_FIT$Antibiotic.x, levels=theorder)


# Distance correlation

# Quantify distance between water and other groups
sdist <- as.data.frame(as.matrix(vegdist(rotu, method = "bray")))
sdist$samplecompare <- row.names(sdist)
row.names(sdist) <- NULL

sdistg <- gather(sdist, sampleid, dist, -samplecompare)
sdistg <- left_join(sdistg, imet, by = c("samplecompare" = "SampleID"))
sdistg <- left_join(sdistg, imet, by = c("sampleid" = "SampleID"))

sdistg <- filter(sdistg, sampleid != "Mock") %>%
	filter(sampleid != "MSA_2002") %>%
	as.data.frame()

sf <- sdistg %>%
	dplyr::filter(grepl("Water", Antibiotic.y.y))

sfg <- sf %>%
	group_by(samplecompare) %>%
	summarise(mean = mean(dist)) %>%
	as.data.frame()

sfg <- left_join(sfg, imet, by = c("samplecompare" = "SampleID")) %>%
	na.omit() %>%
	as.data.frame()

sfgb <- sfg %>%
	dplyr::filter(!grepl("Fungi", Target.x))

sfgbo <- sfgb

sfgbo$Antibiotic.x <- gsub("Baytril_Enrofloxacin", "Baytril", sfgbo$Antibiotic.x)

sfgbo$Antibiotic.x <- factor(sfgbo$Antibiotic.x, levels=theorder)



####################################
# ABBY EDITS FOR MANUSCRIPT FIGURE #
####################################

## 002 BALB/c titer and disruption correlation figure
distordcorplain2 <- ggplot(sfgbo, aes(x = mean, y = OVA_IgG_titer_ng_ml)) +
   theme_classic() +
   geom_smooth(method = "lm", size=0.2, colour="red") +
   geom_point(shape=21, size=1.5, aes(fill = Antibiotic.x), stroke=.3) +
   xlab("Distance From Water") +
   ylab(expression(paste("anti-OVA IgG titer (",mu,"g/ml)"))) +
   theme(text = element_text(size = 8)) +
   #scale_fill_manual(values = c("#000000", "#FF8000", "#9806C8", "#FFFF00", "#FF80C0", "#80FF00", "#0080FF", "#FF0000"), name = "Antibiotic") +
   scale_fill_manual(values = c("#999999","#FF7F00","#AC5782","#FFE528","#E485B7","#449B75","#596A98","#E41A1C")) +
   theme(text = element_text(family="sans")) +
   theme(axis.text = element_text(family="sans", size = 8)) +
   theme(legend.position = "none")
plot(distordcorplain2)
pdf(
  file = "./mAbx-Vax-002/figures/16S_otu_dist_corr_manuscript.pdf",
  width = 3,
  height = 1.7)
distordcorplain2
dev.off()

##002 BALB/c community barplots
mcountm1 <- mcountm
mcountm1$Antibiotic.x <- factor(mcountm1$Antibiotic.x,
    levels = c("Water","Vancomycin","Clindamycin","Neomycin","Ampicillin",
          "Metronidazole","Tobramycin","Baytril"))
genusplotorder2 <- ggplot(mcountm1, aes(x = SampleID, y = relabund, fill = Genus)) +
  theme_classic() +
  theme(
    #legend.position="bottom",
    legend.position="none",
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    strip.text.x = element_blank(),
    axis.line = element_blank()
  ) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = rev(getPalette(11))) +
  ylab("Percent Relative Abundance") +
  xlab("Sample ID") +
  theme(axis.text = element_text(size = 8)) +
  theme(text = element_text(size = 8)) +
  facet_grid(. ~ Antibiotic.x, scales = "free_x")
plot(genusplotorder2)
pdf(
  file = "./mAbx-Vax-002/figures/16S_taxonomy_genus_manuscript.pdf",
  height = 2,
  width = 3
)
genusplotorder2
dev.off()

##bray curtis
otubray2 <- ggplot(ORD_FIT, aes(x = MDS1, y = MDS2)) +
  theme_classic() +
  #geom_point(size = 2.5) +
  geom_point(shape=21, size=1.5, aes(fill = Antibiotic.x), stroke=.3) +
  #geom_point(shape=21, size=5, aes(fill = Antibiotic.x), stroke=.3) +
  theme(axis.text = element_text(family="sans", size = 8),
    text = element_text(size = 8),
    legend.position="none"
  )+ scale_fill_manual(values = c("#999999","#FF7F00","#AC5782","#FFE528","#E485B7","#449B75","#596A98","#E41A1C"))
plot(otubray2)

pdf(
  file = "./mAbx-Vax-002/figures/16S_otu_braycurtis_manuscript.pdf",
  width = 2.8,
  height = 2.2)
otubray2
dev.off()
