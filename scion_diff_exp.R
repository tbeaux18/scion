#!/usr/bin/env Rscript

library(tidyverse)
library(magrittr)
library(DESeq2)

run_scRNAseq_ZINB_DESeq <- function(count_matrix, expdesign, colCondition, alpha){
  
  library(DESeq2)
  library(zinbwave)
  
  dds_list = list()
  
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = expdesign, design=as.formula(colCondition))
  
  
  # filtering criteria: at least 5 samples have at least 5 counts; can be changed
  # filter <- rowSums(counts(dds)) >= 10
  filter <- rowSums(assay(dds) > 5) > 5 
  
  dds <- dds[filter,]
  dds_list[["dds"]] <- dds
  
  zinb <- zinbwave(dds, K=0, BPPARAM=SerialParam(), epsilon=1e12)
  dds_list[["zinb"]] <- zinb
  
  zinb.dds <- DESeqDataSet(zinb, design=as.formula(colCondition))
  
  
  # these settings are recommended by latest literature for single cell testing especially with UMI Counts
  # using LRT is suggested;
  # LRT works by determining p-values by finding the difference in deviance between full and reduced models
  # these are NOT log2fold changes
  zinb.dds <- DESeq(zinb.dds, test="LRT", reduced = ~1, sfType="poscounts", minmu=1e-6, minReplicatesForReplace=Inf)
  dds_list[["zinb.dds"]] <- zinb.dds
  
  # grabs DGE results from dds object
  # zinb.res <- results(zinb.dds, alpha=alpha)
  
  # shrinks log fold change use apeglm; apeglm must be installed
  zinb.res <- lfcShrink(zinb.dds, coef=resultsNames(zinb.dds)[2], type="apeglm")
  # reordering results by adjusted pvalue FDR 
  zinb.res <- zinb.res[order(zinb.res$padj),]
  dds_list[["zinb.res"]] <- zinb.res
  
  return(dds_list)
  
}


# need to set working dir
# setwd("~/Documents/single_cell_pipeline/scion")

# arg1 = count data
# arg2 = cell data
args = commandArgs(trailingOnly=TRUE)

# experimental design must be in the same order as the columns of the umi count matrix!
# umi count matrix is in lexographic order based on cell barcode sequence
# replace 'ExperimentDesign.csv' with the absolute path to your Experiment Design sheet
exp_design <- read.csv('ExperimentDesign.csv', stringsAsFactors = F, header = T)
exp_design %<>% arrange(barcode_sequence)

# All count data from zUMIs pipeline
# all.counts <- readRDS(args[1])
# replace this path with the path to the .dgecounts.rds file outputted from the pipeline
all.counts <- readRDS('~/Documents/single_cell_pipeline/scion/zUMIs_output/expression/dnpf_dev.dgecounts.rds')

# create all-sample UMI count matrix using cell_data sample_names as column names           // TODO: warn if column names repeat - create a table?
# umi counts include intron counts (best to use per zUMIs due to accurate UMI collapsing)
umi.count <- as.matrix(all.counts$umicount$inex$all)

# remove initial matrix to reduce RAM use
rm(all.counts)

# THE BELOW LINE MUST BE TRUE IN ORDER TO CONTINUE
# tests if experimental design order is same as umi count matrix
all(exp_design$barcode_sequence == colnames(umi.count))

# converting umi.count to dataframe to convert colnames to human readable sample names
umi.count <- as.data.frame(umi.count)
names(umi.count) <- exp_design$sample_name[match(names(umi.count), exp_design$barcode_sequence)]
umi.count <- as.matrix(umi.count)




# plotting library sizes for each sample
librarySizes <- as.data.frame(colSums(umi.count))
colnames(librarySizes) <- c("lib_size")
ggplot(librarySizes, aes(x=rownames(librarySizes), y=lib_size)) + 
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell ID", y = "Library Size", title="Library Size (# of Total Counts) per Cell")


# plotting PCA of the original data
library(ggfortify)
rlogcounts <- rlog(as.matrix(umi.count))
# run PCA
pcDat <- prcomp(t(rlogcounts))

# plot PCA
# change x and y to principal component numbers you would like to compare
autoplot(pcDat, x=1, y=2, data = exp_design, colour="condition", shape="gfp_status",size=3)



# filter experiment design sheet to fed to compare GFP pos vs. neg DGE
fed.design <- exp_design %>% filter(condition == 'fed')
fed.cts <- as.data.frame(umi.count) %>% select(fed.design$sample_name) %>% as.matrix()
fed_dds <- run_scRNAseq_ZINB_DESeq(count_matrix=fed.cts, expdesign=fed.design, colCondition=~gfp_status, alpha=0.1)
write.csv(fed_dds$zinb.res, file='fed_dge_results.csv')

# filter experimental design sheet to starved to compare GFP pos vs. neg DGE
starve.design <- exp_design %>% filter(condition == 'starved')
starve.cts <- as.data.frame(umi.count) %>% select(starve.design$sample_name) %>% as.matrix()
starve_dds <- run_scRNAseq_ZINB_DESeq(count_matrix=starve.cts, expdesign=starve.design, colCondition=~gfp_status, alpha=0.1)
write.csv(starve_dds$zinb.res, file='starved_dge_results.csv')

# filter experimental design to just GFP positive cells to compare fed vs. starved conditions
gfp.design <- exp_design %>% filter(gfp_status == 'positive')
gfp.cts <- as.data.frame(umi.count) %>% select(gfp.design$sample_name) %>% as.matrix()
gfp_dds <- run_scRNAseq_ZINB_DESeq(count_matrix=gfp.cts, expdesign=gfp.design, colCondition=~condition, alpha=0.1)
write.csv(gfp_dds$zinb.res, file='gfp_dge_results.csv')


# heat map of top 35 genes plotted for each cell. 
# the amount by which each gene deviates in a specific sample from the gene’s average across all samples.
library(genefilter)
library(gplots)
rld <- rlog(gfp_dds$zinb.dds)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35)
heatmap.2( assay(rld)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="column", cexCol = 0.5,
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
