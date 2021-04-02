#!/usr/bin/env Rscript
# Copyright © 2021 Merck Sharp & Dohme Corp. a subsidiary of Merck & Co., Inc., Kenilworth, NJ, USA.
# Geoffrey Hannigan
# Systems Biology Research Group
# Merck Exploratory Science Center, Cambridge

##################
# Load Libraries #
##################

library(ggplot2)
library(tidyr)
library(dplyr)
library(vegan)
library(stringr)
library(RColorBrewer)
#library(Rmisc)
library(cowplot)

################
# Load In Data #
################

# Load in the sequencing depth data
input_depth <- read.delim(
	file = "./NHP-Ab-Vax-001/data/cosmosid_run/data_from_online_portal/read_depth.tsv",
	sep = "\t",
	header = TRUE
	)

input_counts <- read.delim(
	file = "./NHP-Ab-Vax-001/data/cosmosid_run/data_from_online_portal/count_table.tsv",
	sep = "\t",
	header = TRUE
	) %>%
	select(-Taxonomy.Lineage, -Taxonomy.Id)

# Let's make a taxonomy table
taxonomytable <- input_counts %>%
	select(Kingdom, Phylum, Class, Order, Family, Genus, Species, Name)

# Remove from the count table
trimmed_counts <- input_counts %>%
	select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species)

# Transpose the count table for rarefaction
# Remember the names
remember_names <- trimmed_counts$Name
ttrim <- as.data.frame(t(trimmed_counts[,-1]))
colnames(ttrim) <- remember_names
# Add a non-specific filter
ttrim <- ttrim[,colSums(ttrim) > 5000]

metadata <- read.delim(
	file = "./NHP-Ab-Vax-001/data/metadata/metadata.tsv",
	sep = "\t",
	header = TRUE
	)

##############################
# Analysis - Alpha Diversity #
##############################

# Get richness
# Make sure to rarefy dataset
rotu <- rrarefy(ttrim, sample = 30000)
length(rownames(rotu))
rotu <- rotu[c(rowSums(rotu) >= 30000),]
rotu[ rotu < 10 ] <- 0
length(rownames(rotu))

# Lets look at richness here, and compare to controls
dotu <- as.data.frame(specnumber(rotu))
colnames(dotu) <- "Richness"
dotu$sampleid <- row.names(dotu)
# Remove the merck label
dotu$sampleid <- gsub("_MER1264", "", dotu$sampleid)

dmmerge <- left_join(dotu, metadata, by = c("sampleid" = "SampleID"))

dmmerge <- dmmerge

# Organize the days
dmmerge$Day <- gsub("Day_", "", dmmerge$Day)
dmmerge$Day <- as.numeric(dmmerge$Day)
dmmerge <- dmmerge %>%
	arrange(Day)
dmmerge$Day <- as.factor(dmmerge$Day)


#############################
# Analysis - Beta Diversity #
#############################


odist <- vegdist(rotu, method = "bray")

ORD_NMDS <- metaMDS(comm = odist, k=2, trymax = 100)
ORD_FIT = data.frame(MDS1 = ORD_NMDS$points[,1], MDS2 = ORD_NMDS$points[,2])
ORD_FIT$SampleID <- rownames(ORD_FIT)
# Remove the merck label
ORD_FIT$SampleID <- gsub("_MER1264", "", ORD_FIT$SampleID)
# Get metadata
ORD_FIT_M <- left_join(ORD_FIT, metadata, by = c("SampleID" = "SampleID")) %>%
	na.omit()

ORD_FIT_M$Day <- gsub("Day_", "", ORD_FIT_M$Day)
ORD_FIT_M$Day <- as.numeric(ORD_FIT_M$Day)
ORD_FIT_M <- ORD_FIT_M %>%
	arrange(Day)
ORD_FIT_M$Day <- as.factor(ORD_FIT_M$Day)

# Mark groups for coloring
ORD_FIT_M$Class <- "No Antibiotics"
ORD_FIT_M$Class[which(ORD_FIT_M$Group == 4 & ORD_FIT_M$Day == 0)] = "Day 0 Post Vancomycin"
ORD_FIT_M$Class[which(ORD_FIT_M$Group == 4 & ORD_FIT_M$Day == 28)] = "Day 28 Post Vancomycin"
ORD_FIT_M$Class[which(ORD_FIT_M$Group == 4 & ORD_FIT_M$Day == 84)] = "Day 83 Post Vancomycin"
ORD_FIT_M$Class[which(ORD_FIT_M$Group == 4 & ORD_FIT_M$Day == 173)] = "Day 173 Without Vancomycin"

# Order the values
ORD_FIT_M$Class <- factor(
	ORD_FIT_M$Class,
	levels = c(
		"No Antibiotics",
		"Day 0 Post Vancomycin",
		"Day 28 Post Vancomycin",
		"Day 83 Post Vancomycin",
		"Day 173 Without Vancomycin"
	)
)


ord_plot2 <- ORD_FIT_M %>%
  filter(Group == 1 | Group == 4) %>%
  ggplot(aes(x = MDS1, y = MDS2)) +
  theme_classic() +
  #geom_point(size = 2.5) +
  geom_point(shape=21, size=1.5, aes(fill = Class), stroke=.3) +
  scale_fill_manual(values = c(
    "#596A98",
    "#E41A1C",
    "#E485B7",
    "#FF7F00",
    "#449B75")
  )+ #name = "Treatment") +
  theme(legend.position="none") +
  scale_shape_discrete(name = "Group") +
  theme(axis.text = element_text(family="sans", size = 8),
      text = element_text(size = 8))
plot(ord_plot2)


pdf(
  file = "./NHP-Ab-Vax-001/figures/bray_curtis_manuscript.pdf",
  height = 2.2,
  width = 2.8
)
ord_plot2
dev.off()



#################################
# Plotting - Combine With Titer #
#################################

# Titer
tplot <- read.delim(
	file = "./NHP-Ab-Vax-001/data/metadata/nhp_full_titer.tsv",
	sep = "\t",
	header = TRUE
	)

# Make titer numeric
tplot$Titer <- as.numeric(as.character(tplot$Titer))

# Remove the groups that we don't want
tplot <- tplot %>%
	filter(Treatment != "Combo") %>%
	as.data.frame()

# Make the days numeric
tplot$Day <- as.numeric(as.character(tplot$Day))
tplot$Treatment <- gsub("Vanco", "Vancomycin", tplot$Treatment)

# Plot out the titers over time
splot <- summarySE(tplot, measurevar="Titer", groupvars=c("Day","Treatment"))
pd <- position_dodge(1)

xxis <- scale_x_continuous(limits = c(-20, 190), breaks = round(seq(min(lsplot$Day), max(lsplot$Day), by = 20),1))
vvline <- geom_vline(xintercept = c(0, 28, 83, 173), linetype = "dashed", alpha = 0.75, colour = "grey")

forlegend <- ggplot(splot, aes(x=Day, y=Titer, colour=Treatment)) +
	theme_classic() +
	geom_errorbar(aes(ymin=Titer-se, ymax=Titer+se), width=1, position=pd) +
	geom_line(position=pd) +
	geom_point(position=pd) +
	xxis +
	scale_y_log10() +
	ylab("Flu Endpoint Titer") +
	scale_colour_manual(values = c(brewer.pal(3,"Set1")[2], brewer.pal(3,"Set1")[1])) +
	vvline

titer_plot <- ggplot(splot, aes(x=Day, y=Titer, colour=Treatment)) +
	theme_classic() +
	geom_errorbar(aes(ymin=Titer-se, ymax=Titer+se), width=1, position=pd) +
	geom_line(position=pd) +
	geom_point(position=pd) +
	xxis +
	scale_y_log10() +
	ylab("Flu Endpoint Titer") +
	scale_colour_manual(values = c(brewer.pal(3,"Set1")[2], brewer.pal(3,"Set1")[1])) +
	theme(legend.position = "none",
		axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
	vvline

# Richness
line_cow <- ggplot(lsplot, aes(x=Day, y=Richness, colour=Group)) +
	theme_classic() +
	geom_errorbar(aes(ymin=Richness-se, ymax=Richness+se), width=1, position=lpd) +
	geom_line(position=lpd) +
	geom_point(position=lpd) +
	xxis +
	scale_y_continuous(limits = c(0,90)) +
	ylab("Richness") +
	scale_colour_manual(values = c(brewer.pal(3,"Set1")[2], brewer.pal(3,"Set1")[1])) +
	theme(
		legend.position = "none",
		axis.title.x=element_blank()) +
	vvline

legend <- get_legend(
	forlegend
)

firstplot <- plot_grid(titer_plot, line_cow, align = "v", ncol = 1, rel_heights = c(0.9, 1))
legend_plus <- plot_grid(firstplot, legend, rel_widths = c(3, .45))

pdf(
	file = "./NHP-Ab-Vax-001/figures/stacked_line_div_titer.pdf",
	height = 5,
	width = 8
	)

	legend_plus

dev.off()
