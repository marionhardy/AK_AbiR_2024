
### This script was written by Marion Hardy for a collaboration with 
## Akshaya Karthikeyan from the Lombard lab

## This sets up your DESeq2 object for further data analysis
## This is the AbiR cell line

library(tidyverse)
library(readxl)
library(DESeq2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)

## Load the raw count matrices, two in this case

counts = read.table("./data/rawcountdata_AbiR.txt", header = T)

## Load metadata (sample annotation file)

coldata = read_xlsx("./data/metadata.xlsx", sheet = "Recoded")
rownames(coldata) = coldata$ID

## Clean up counts matrix and ensembl id

strrep = sub(pattern = "\\.(.*)","",counts$X)
counts$X = strrep
rownames(counts) = counts$X
counts = counts %>% select(!c(X,gene_name))

# Create the full model for comparison of samples
# AK said compare DMSO to Ola day 5 and Ola day 9

dds = DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                              design = ~Condition) 
# checked experimental design, looks good, no recorded confounding variable

# Generate a linear model

dds = DESeq(dds)
resultsNames(dds)

## Checking distribution of counts per sample

as_tibble(assay(dds)) %>%
  gather(sample, value = counts) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ sample) # looks ok

# Checking size factors and dispersion

sizeFactors(dds) # only takes into account the sequencing depth, looks ok
plotDispEsts(dds) # verifies normalization, graph looks a-ok

# Checking PCA

rld = vst(dds)

plotPCA(rld,intgroup="Condition") + 
  theme_bw()+
  labs(title = 'PCA per condition')

# Checking sample similarity

sampleDists <- dist(t(assay(dds)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(dds$Condition, sep="-")
colnames(sampleDistMatrix) <- paste(dds$Condition, sep="-")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# you can see one of the day9 is clustering with day5
# but the pca shows that very little of the variance is explained
# by the length of Ola treatment so it makes sense that they are similar
# enough to cluster together when computing euclidian distances


## Saving the DESeq object

saveRDS(dds, "./data_output/AbiR_DESeq.Rds")
