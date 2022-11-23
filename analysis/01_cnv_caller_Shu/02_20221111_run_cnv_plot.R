# modified from Shu's script modifed from 2019 NBT Massively parallel single-cell chromatin landscapes...
# https://github.com/GreenleafLab/10x-scATAC-2019/tree/master/code 01_Filter_Cells_v2.R and 08_Run_scCNV_v2.R
library(ArchR)
library(ComplexHeatmap)
set.seed(1)
addArchRGenome("hg38")
library(magrittr)
library(ggplot2)
library(Rcpp)
library(viridis)
library(Matrix)
library(SummarizedExperiment)
library(matrixStats)
library(readr)
library(GenomicRanges)
library(edgeR)
library(Seurat)
packageVersion('Seurat')
library(BSgenome.Hsapiens.UCSC.hg38)
source('/Users/xwyan/Dropbox/LAB/github/2022_cnv/cnv/analysis/01_cnv_caller_Shu/03_20221122_functions.r')

# ----------------------------
# Get Inputs
# ----------------------------
input_dir <- c("/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_scMultiome_ColoDM-ColoHSR/COLO320HSR_5k/03_analysis/")
save_dir <- c("/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_scMultiome_ColoDM-ColoHSR/COLO320HSR_5k/03_analysis/")
sample <- c("COLO320HSR_5k_rep1")
chr = 8
used.chr<- paste0(c('chr'), chr)

# ----------------------------
# Window generation
# ----------------------------
# blacklist was downloaded from:https://github.com/Boyle-Lab/Blacklist/tree/master/lists
# blacklist was originally described in paper The ENCODE Blacklist: Identification of Problematic Regions of the Genome (June 2019)
# blacklist <- import.bed("/Users/xwyan/Dropbox/LAB/Seqbasic/resources/blacklist/hg38-blacklist.v2.bed")
# windows <- makeWindows(genome = BSgenome.Hsapiens.UCSC.hg38, blacklist = blacklist)
# saveRDS(windows, paste0(save_dir, 'hg38_blacklist_window1e6_sliding2e5.rds'))

# ----------------------------
# Window and ATAC reads loading
# ----------------------------
windows1 <- readRDS(paste0(input_dir, 'hg38_blacklist_window1e6_sliding2e5.rds'))
windows1_chr <- windows1[windows1@seqnames %in% used.chr,]

fragments <- readRDS(paste0(input_dir, sample, "_chr", chr, ".gr.rds"))

# only run cnv with chr in ecDNA
cnaObj_dm2 <- scCNA_dm(windows1_chr, fragments, neighbors = 1000, LFC = 1.25, FDR = 0.1, force = TRUE, used.chr = used.chr) 
cnaObj <- cnaObj_dm2
data_log2FC <- t(cnaObj@assays@data@listData$log2FC)
data_counts <- t(cnaObj@assays@data@listData$counts)
rownames(data_log2FC) <- rownames(cnaObj@colData)
colnames(data_log2FC) <- cnaObj@rowRanges@ranges@start
rownames(data_counts) <- rownames(cnaObj@colData)
colnames(data_counts) <- cnaObj@rowRanges@ranges@start
df_log2FC <- data.frame(data_log2FC)
df_counts <- data.frame(data_counts)
write.table(df_log2FC, paste0(save_dir, sample, "_window1e6_sliding2e5_chr8_log2FC.txt"), sep='\t')
write.table(df_counts, paste0(save_dir, sample, "_window1e6_sliding2e5_chr8_counts.txt"), sep='\t')



