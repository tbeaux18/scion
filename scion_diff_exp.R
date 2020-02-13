#!/usr/bin/env Rscript

library(tidyverse)
library(magrittr)
library(DESeq2)

run_scRNAseq_ZINB_DESeq <- function(count_matrix, expdesign, colCondition, alpha){
  
  library(DESeq2)
  library(zinbwave)
  
  dds_list = list()
  
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = expdesign, design=as.formula(colCondition))
  
  
  # filtering criteria: the sum of all counts down row must be greater than or equal to 5, regardless of number of samples
  filter <- rowSums(counts(dds)) >= 5
  # filter <- rowSums(assay(dds) > 1) > 5
  
  dds <- dds[filter,]
  dds_list[["dds"]] <- dds
  
  # params are based on current literature
  zinb <- zinbwave(dds, K=0, BPPARAM=SerialParam(), epsilon=1e12)
  dds_list[["zinb"]] <- zinb
  
  zinb.dds <- DESeqDataSet(zinb, design=as.formula(colCondition))
  
  
  # these settings are recommended by latest literature for single cell testing especially with UMI Counts
  # using LRT is suggested;
  # LRT works by determining p-values by finding the difference in deviance between full and reduced models
  # these are NOT testing log2fold ratios
  # The LRT examines two
  # models for the counts, a full model with a certain number of terms and a reduced model, in which some of
  # the terms of the full model are removed. The test determines if the increased likelihood of the data using the
  # extra terms in the full model is more than expected if those extra terms are truly zero.
  zinb.dds <- DESeq(zinb.dds, test="LRT", reduced = ~1, sfType="poscounts", minmu=1e-6, minReplicatesForReplace=Inf)
  p <- plotDispEsts(zinb.dds)
  dds_list[["zinb.dds"]] <- zinb.dds
  dds_list[["dispersion.plot"]] <- p
  
  rld <- rlog(zinb.dds)
  dds_list[["rld"]] <- rld
  
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



library(ggfortify)
# plotting PCA of the original data
rlogcounts <- rlog(as.matrix(umi.count))
# run PCA
pcDat <- prcomp(t(rlogcounts))
# plot PCA
# change x and y to principal component numbers you would like to compare
autoplot(pcDat, x=1, y=2, data = exp_design, colour="condition", shape="gfp_status",size=3)



# filter experimental design to just GFP positive cells to compare fed vs. starved conditions
gfp.design <- exp_design %>% filter(gfp_status == 'positive')
gfp.cts <- as.data.frame(umi.count) %>% select(gfp.design$sample_name) %>% as.matrix()
gfp_dds <- run_scRNAseq_ZINB_DESeq(count_matrix=gfp.cts, expdesign=gfp.design, colCondition=~condition, alpha=0.1)
#gfp_dds_5 <- run_scRNAseq_ZINB_DESeq(count_matrix=gfp.cts, expdesign=gfp.design, colCondition=~condition, alpha=0.05)
write.csv(gfp_dds$zinb.res, file='gfp_dge_results.csv')


### reduced based on one starved cell showing more counts than most other cells
gfp.design.reduced <- gfp.design %>% filter(sample_name != 'sta-gfp-2L(S2.2)')
gfp.cts.reduced <- as.data.frame(umi.count) %>% select(gfp.design.reduced$sample_name) %>% as.matrix()
gfp_dds.reduced <- run_scRNAseq_ZINB_DESeq(count_matrix=gfp.cts.reduced, expdesign=gfp.design.reduced, colCondition=~condition, alpha=0.1)



library(genefilter)
library(gplots)
library(RColorBrewer)


# shows heat maps of top 35 highly expressed genes for GFP positive between starve vs. fed states
select <- order(rowMeans(counts(gfp_dds$zinb.dds, normalized=TRUE)),decreasing=TRUE)[1:35]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

# heat map of raw counts
heatmap.2(counts(gfp_dds$zinb.dds,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10,6))
# heat map of regularized log transformation counts
heatmap.2(assay(gfp_dds$rld)[select,], col = hmcol,
          Rowv = TRUE, Colv = TRUE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))



# no Ilp2 or Ilp1
# make a character vector such as this with the names of the genes
# header for column 1 must be gene_name
# you can replace the neuro_peptide.csv with whatever the csv filename of the new gene name set
neuro_df <- read.csv('neuro_peptide.csv', stringsAsFactors = F, header = T)
neuro_df$gene_name <- toupper(neuro_df$gene_name)
# creates a blue to red color palette
neuro.hmcol <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)

# finds the indices of the desired genes in the character vector you make
# make sure to change the neuro_list to whatever character vector variable name
idx <- which(toupper(rownames(assay(gfp_dds$rld))) %in% neuro_df$gene_name, arr.ind=TRUE)
# generates the heat map in order
# need to fix the number of genes showing up on the y axis

neuro_gene_rld <- assay(gfp_dds$rld)[idx,]
sort_col_names <- sort(colnames(neuro_gene_rld))
neuro_gene_rld_sorted <- neuro_gene_rld[,sort_col_names]
neuro_ordered_rowmean <- neuro_gene_rld_sorted[order(rowMeans(neuro_gene_rld_sorted), decreasing = TRUE), ]

neuro.hm <- heatmap.2(neuro_ordered_rowmean, col = neuro.hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6), symbreaks = FALSE)

# regularized log transformation
# remove the dependence of the variance on the mean, 
# particularly the high variance of the logarithm of count data when the mean is low.
# rowvars fxn grabs variance of numeric array
topVarGenes <- head(order(rowVars(assay(gfp_dds$rld)), decreasing=TRUE ), 35)

# generates heat map based on variance
# heat map of top 35 genes plotted for each cell. 
# the amount by which each gene deviates in a specific sample from the gene’s average across all samples
heatmap.2( assay(gfp_dds$rld)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="column", cexCol = 0.5,
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))


# heat map of sample to sample distance
distsRL <- dist(t(assay(gfp_dds$rld)))
mat <- as.matrix(distsRL)
hc <- hclust(distsRL)
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none",
          col = rev(hmcol), margin=c(13, 13))


# plotting PCA of GFP positive fed vs. starved
# can only plot PC1 vs. PC2; function from DESeq2 library
plotPCA(gfp_dds$rld, intgroup=c("condition"))



##############################################################################

# BELOW IS JUST QC RUNS COMPARING GFP STATUS

##############################################################################

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