#!/usr/bin/env Rscript

library(tidyverse)
library(magrittr)
library(DESeq2)
library(DEGreport)

# need to set working dir
# setwd("~/Documents/single_cell_pipeline/scion")

# arg1 = count data
# arg2 = cell data
args = commandArgs(trailingOnly=TRUE)

# experimental design must be in the same order as the columns of the umi count matrix!
# umi count matrix is in lexographic order based on cell barcode sequence
exp_design <- read.csv('ExperimentDesign.csv', stringsAsFactors = F, header = T)
exp_design %<>% arrange(barcode_sequence)

# All count data from zUMIs pipeline
# all.counts <- readRDS(args[1])
all.counts <- readRDS('~/Documents/single_cell_pipeline/scion/zUMIs_output/expression/dnpf_dev.dgecounts.rds')

# create all-sample UMI count matrix using cell_data sample_names as column names           // TODO: warn if column names repeat - create a table?
# umi counts include intron counts (best to use per zUMIs due to accurate UMI collapsing)
umi.count <- as.matrix(all.counts$umicount$inex$all)

# remove initial matrix
rm(all.counts)

# converting umi.count to dataframe to convert colnames to human readable sample names
umi.count <- as.data.frame(umi.count)
names(umi.count) <- exp_design$sample_name[match(names(umi.count), exp_design$barcode_sequence)]
umi.count <- as.matrix(umi.count)

# THE BELOW LINE MUST BE TRUE IN ORDER TO CONTINUE
# tests if experimental design order is same as umi count matrix
all(exp_design$sample_name == colnames(umi.count))


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
autoplot(pcDat,
         data = exp_design, 
         colour="condition", 
         shape="gfp_status",
         size=3)


# Get log2 counts per million
logcounts <- log2(umi.count + 1)
# make a colour vector
statusCol <- as.numeric(factor(exp_design$condition)) + 1
# Check distributions of samples using boxplots
logstack <- stack(as.data.frame(logcounts))
ggplot(logstack) + 
  geom_boxplot(aes(x=ind, y=values)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell ID", y = "Log2(Counts)")

boxplot(logcounts, 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,
        col=statusCol)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(as.matrix(logcounts)), col="blue")





run_scRNAseq_ZINB_DESeq <- function(count_matrix, expdesign, colCondition){
  
  library(DESeq2)
  library(zinbwave)
  
  dds_list = list()
  
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = expdesign, design=as.formula(colCondition))
  

  # filtering criteria: at least 5 samples have at least 5 counts
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
  # zinb.res <- results(zinb.dds)
  
  # shrinks log fold change use apeglm; apeglm must be installed
  zinb.res <- lfcShrink(zinb.dds, coef=resultsNames(zinb.dds)[2], type="apeglm")
  zinb.res <- zinb.res[order(zinb.res$padj),]
  dds_list[["zinb.res"]] <- zinb.res
  return(dds_list)
  
}


# filtering to just GFP positive cells
fed.design <- exp_design %>% filter(condition == 'fed')
fed.cts <- as.data.frame(umi.count) %>% select(fed.design$sample_name) %>% as.matrix()

fed_dds <- scRNAseq_ZINB_DESeq(count_matrix=fed.cts, expdesign=fed.design, colCondition=~gfp_status)

# filtering to just GFP positive cells
gfp.design <- exp_design %>% filter(gfp_status == 'positive')
gfp.umi.count <- as.data.frame(umi.count) %>% select(gfp.design$sample_name) %>% as.matrix()

gfp_dds <- scRNAseq_ZINB_DESeq(count_matrix=fed.cts, expdesign=gfp.design, colCondition=~condition)


plotMA(zinb.res, ylim=c(-3,3), cex=.8)
abline(h=c(-1,1), col="dodgerblue", lwd=2)

mcols(dds)

res <- results(zinb.dds, independentFiltering=FALSE)
plot(mcols(zinb.dds)$condition_starved_vs_fed, res$log2FoldChange, ylim=c(-4,4))
abline(0,1,col="red")


