#!/usr/bin/env Rscript
# Copyright © 2021 Merck Sharp & Dohme Corp. a subsidiary of Merck & Co., Inc., Kenilworth, NJ, USA.
# Geoffrey Hannigan edited by Begüm Topcuoglu
# Systems Biology Research Group
# Merck Exploratory Science Center, Cambridge

#################
# Set Libraries #
#################

deps = c("vegan", "pheatmap", "dplyr", "cowplot", "viridis" ,"ggplot2", "stringr", "data.table", "tibble","parallel", "tidyverse", "tibble");
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE, repos = 'http://cran.rstudio.com/', dependencies=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

#######################################
# Load In and Preprocess RNAseq Data #
########################################


#----------------------------------------------------------------------------#
load("./mAbx-Vax-005/data/host_rna_seq/VST_Data.RData")
# You can see the transformed count table here, for correlations
vData[1:5,1:5]

VD <- as.data.frame(t(vData))[-1,]

colnames(VD) <- c(t(VD[1,]))

# Remove the columns without names (NA)
VD <- VD[!is.na(names(VD))]


VD <- VD[-1,] %>%
	rownames_to_column("names")

colnames(VD) <- gsub("^", "", names(VD))

VDD <- data.frame(sapply(VD, function(x) as.numeric(as.character(x))))
VDD$names <- VD$names
#----------------------------------------------------------------------------#

#######################################
# Load In RNAseq MetaData #
########################################

# Bring in the metadata for this study

rna_metadata <- as.data.frame(fread(
	file="./mAbx-Vax-005/data/host_rna_seq/ED_mAbx_005.csv",
	header=TRUE,
	sep=","))
# Fix redundant column IDs
colnames(rna_metadata) <- make.unique(colnames(rna_metadata))



#######################################
# Load In 16S rRNA data #
########################################

### 16S RRNA SEQ
# Import the OTU abundance table
# Pull in the large table quickly with fread
otucounts <- as.data.frame(fread(
	file="./mAbx-Vax-005/data/raw/OtuCounts.tsv",
	header=TRUE,
	sep="\t"))

# Import the OTU taxonomic assignments (down to genus level)
otutaxonomy <- read.delim(
	file="./mAbx-Vax-005/data/raw/OtuTaxonomy.tsv",
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


###############################################
# Load In and preprocess 16S rRNA MetaData #
################################################


micro_metadata <- read.delim(
	file="./mAbx-Vax-005/data/metadata/mabxvax005_metadata.tsv",
	header=TRUE,
	sep="\t")

# Make the animal number format same as RNA-seq
micro_metadata$CDonor <- gsub(".*D0\\.CC", "M", micro_metadata$SampleID, perl = TRUE)


#################################
# OTU Count Table Normalization #
#################################


formotu <- otucounts
row.names(formotu) <- otucounts[,2]
formotu <- formotu[,-c(1:3)]

# There is just way too much excess OTU information here
# Cut those out with very few counts overall
formotu <- formotu[,colSums(formotu) > 10]

# What does the depth distribution look like?
qplot(rowSums(formotu))

# Get richness and Shannon diversity
# Make sure to rarefy dataset
rotu <- rrarefy(formotu, sample = 10000)
length(rownames(rotu))
rotu <- rotu[c(rowSums(rotu) >= 10000),]
rotu[ rotu < 2 ] <- 0
length(rownames(rotu))
# Make it a data frame
names_rotu <- rotu %>%
	as.data.frame() %>%
	rownames_to_column("names") %>%
	as.data.frame()

otu_m <- left_join(names_rotu, micro_metadata, by = c("names" = "SampleID"))
# Sanity check (a couple were lost in rarefaction)
nrow(names_rotu)
nrow(micro_metadata)
nrow(otu_m)

otu_f <- otu_m %>%
	filter(SampleType == "Cecal_Content") %>%
	filter(Antibiotics %in% c("Water", "Vancomycin"))
str(otu_f)
otu_final <- otu_f[, -(2975:2987)]


tidy_data <- gather(otu_final, "OTU", "Count", 2:2974)

source("mAbx-Vax-005/bin/relative_abundances.R")

otu_final <- get_relative_count(tidy_data)

spread <- otu_final %>%
  select(names, OTU, relabund)%>%
	spread(key = OTU, value = relabund, fill = 0)


write.csv(spread , "mAbx-Vax-005/data/otu_relative_abundances.csv")

########################################################
# Prepare OTU counts for correlation matrix use  #
########################################################
# Set variance lower limit
variancevar <- 10


# Get list of higher variance OTUs
otus_only <- otu_f[,grep("Otu", colnames(otu_f))]

# Extract the list of OTU names to build corr matrix
OTU_names <- colnames(otus_only[ , which(sapply(otus_only, var) > variancevar)])

otu_labels <- oparse %>%
  filter(OTU %in% OTU_names) %>%
  select(OTU, Genus)

write.csv(otu_f, "mAbx-Vax-005/data/otu_abundances.csv")

########################################################
# Prepare RNAseq counts for correlation matrix use  #
########################################################

# Get the gene names Corey identified as high logFC and important for innate immune system
# Abby picked them accoruding to functional importance


genes <- read.delim(
	file="./mAbx-Vax-005/data/cor_data/gene_list.csv",
	header=TRUE,
	sep=",")

gene_names_intestine_abby <- c("Klf15", "Tle1")
gene_names_spleen_abby <- c("Cxcl5", "Zbtb16")
gene_names_wb_abby <- c("Cxcr6", "Ptpn22", "Crispld2")
gene_names_lymphnode_abby <- as.character(genes$Mesenteric.Lymph.Nodes)
length(gene_names_lymphnode_abby)
# Make a merged RNA seq file as well
VDM <- left_join(VDD, rna_metadata, by = c("names" = "ShortName"))

########################################################
# Build correlation matrix for each site #
########################################################


build_corr <- function(site_name, gene_names, treatment_group1, treatment_group2){

  #--------------- Prepare Datasets ----------------------##
  VDMF <- VDM %>%
	 filter(Tissue == site_name) %>%
   filter(Treatment %in% c(treatment_group1, treatment_group2))


 #--------------- Build correlation matrix ----------------------##

   # Pearson correlation with FDR correction

   genecor <- mclapply(OTU_names, mc.cores = 8, function(i) {
   	xx <- otu_f %>% select_(i) %>% as.matrix() %>% as.numeric()
   	write(i, stderr())
   		midout <- lapply(gene_names, function(j) {
   			write(j, stderr())
   			yy <- VDMF %>% select(j) %>% as.matrix() %>% as.numeric()
   			oc <- cor.test(xx, yy, method = "spearman", adjust = "BH") # FDR corrected P value
   			pval <- oc$p.value
   			cest <- oc$estimate
   			xout <- data.frame(
   				otuname = i,
   				genename = j,
   				pvalue = pval,
   				correlation = cest
   				)
   			return(xout)
   		})
   		fout <- do.call(rbind, midout)
   		return(fout)
   })

   return(genecor)
}


#----------------#
# SMALL INTESTINE
#----------------#
seed=0
# Create a data.matrix with correlation results (coef and adj.pvalue)
gc_bind_small_int <- do.call(rbind, build_corr("Small_Intestine", gene_names_intestine_abby, "WA", "VA")) %>%
  inner_join(otu_labels, by=c("otuname" = "OTU"))%>%
  mutate(site = "Small_intestine") %>%
  unite(otuname, Genus, otuname, remove=FALSE, sep=" (")  %>%
  arrange(desc(otuname)) %>%
  mutate(otuname = paste(otuname, ")", sep="")) %>%
  mutate(otuname=paste0(gsub('TU0*', 'TU ', otuname)))

#----------------#
# LYMPH NODE
#----------------#

# Create a data.matrix with correlation results (coef and adj.pvalue)
gc_bind_lymphnode <- do.call(rbind, build_corr("Mesenteric_Lymph_Node", gene_names_lymphnode_abby, "WA", "VA")) %>%
  inner_join(otu_labels, by=c("otuname" = "OTU"))%>%
  mutate(site = "Lymph_Node") %>%
  unite(otuname, Genus, otuname, remove=FALSE, sep=" (") %>%
  arrange(desc(otuname)) %>%
  mutate(otuname = paste(otuname, ")", sep="")) %>%
  mutate(otuname=paste0(gsub('TU0*', 'TU ', otuname)))


#----------------#
# SPLEEN
#----------------#

# Create a data.matrix with correlation results (coef and adj.pvalue)
gc_bind_spleen <- do.call(rbind, build_corr("Spleen", gene_names_spleen_abby, "WA", "VA")) %>%
  inner_join(otu_labels, by=c("otuname" = "OTU"))%>%
  mutate(site = "Spleen") %>%
  unite(otuname, Genus, otuname, remove=FALSE, sep=" (")  %>%
  arrange(desc(otuname)) %>%
  mutate(otuname = paste(otuname, ")", sep="")) %>%
  mutate(otuname=paste0(gsub('TU0*', 'TU ', otuname)))


#----------------#
# WHOLE BLOOD
#----------------#

# Create a data.matrix with correlation results (coef and adj.pvalue)
gc_bind_wb <- do.call(rbind, build_corr("Whole_Blood", gene_names_wb_abby, "WA", "VA")) %>%
  inner_join(otu_labels, by=c("otuname" = "OTU"))%>%
  mutate(site = "Whole_Blood") %>%
  unite(otuname, Genus, otuname, remove=FALSE, sep=" (")  %>%
  arrange(desc(otuname)) %>%
  mutate(otuname = paste(otuname, ")", sep="")) %>%
  mutate(otuname=paste0(gsub('TU0*', 'TU ', otuname)))

all <-  rbind(gc_bind_small_int, gc_bind_lymphnode, gc_bind_spleen, gc_bind_wb)

#----------------#
# PLOT FOR MANUSCRIPT
#----------------#

all$Genus<- str_remove(all$Genus, "_unclassified")

final <- all %>%
  unite(gene_name_location, genename, site, sep=" -- ", remove=FALSE)%>%
  unite(gene_name_location_otu, gene_name_location, otuname, sep=" -- ", remove=FALSE)%>%
  group_by(Genus)

filtered <- final %>%
filter(pvalue < 0.005)%>%
filter(correlation > 0.8 | correlation < -0.8) %>%
group_by(otuname)%>%
tally() %>%
filter(n>51)

data <- final %>%
filter(pvalue < 0.005) %>%
filter(correlation > 0.8 | correlation < -0.8) %>%
filter(otuname %in% filtered$otuname)

more_filtered <- data %>%
group_by(gene_name_location)%>%
tally() %>%
filter(n>5)

done_filtered <- data %>%
filter(gene_name_location %in% more_filtered$gene_name_location)


# Make a test matrix and format it appropriately
testmatrix <- done_filtered %>%
  select(gene_name_location, correlation, otuname) %>%
  spread(key = gene_name_location, value = correlation) %>%
  ungroup()%>%
  select(-Genus)

# Make NA values 0
testmatrix[is.na(testmatrix)] <- 0

# Keep the otunames as row names
row.names(testmatrix) <- testmatrix$otuname

# Keep the site as an annotation
t <- data.frame(t(testmatrix))
t <- rownames_to_column(t)
t$site_gene <- t$rowname

b <- separate(t, site_gene, c("gene", "site"), sep = " -- ", remove = TRUE)
b <- b[-1,]
annotation <- data.frame(Site = b$site)

# Remove all columns except correlation coefficients
testmatrix <- testmatrix %>%
  select(-otuname)

# Add the rownames to the annotation to recognize the gene names and sites
row.names(annotation) <- colnames(testmatrix)

# Get the gene names only
labels_row <- data.frame(Gene = b$gene)
row.names(labels_row) <- colnames(testmatrix)


outplot <- pheatmap(t(testmatrix), cluster_rows=T, cluster_cols=F, show_colnames=T, fontsize=4, clustering_method = "average", col = colorRampPalette(c("navy", "white", "firebrick3"))(20001), cellwidth=4, cellheight=4, border_color="gray", annotation_row = annotation, labels_row = labels_row$Gene)

pdf("./mAbx-Vax-005/figures/clustered_heatmap.pdf",
  width = 6,
  height = 10)

  outplot

dev.off()
