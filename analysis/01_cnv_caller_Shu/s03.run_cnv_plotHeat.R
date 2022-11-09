#20210721 zs#
#modified from 2019 NBT Massively parallel single-cell chromatin landscapes...
#https://github.com/GreenleafLab/10x-scATAC-2019/tree/master/code 01_Filter_Cells_v2.R and 08_Run_scCNV_v2.R
library(ArchR)
library(ComplexHeatmap)
set.seed(1)
addArchRGenome("hg19")
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
library(BSgenome.Hsapiens.UCSC.hg19)

#-----------------
# Reading Fragment Files
#-----------------
#message("Reading in fragment files...")
#file_fragments <- read.table('input_fragments.txt')
#libname <- read.table('input_libName.txt')
#listOfDataFrames<-rep(list(NA),nrow(file_fragments))
#for(i in 1:nrow(file_fragments)){
#listOfDataFrames[[i]] <- data.frame(readr::read_tsv(file_fragments[i,1], col_names=FALSE))
#name<-libname[i,1]
#listOfDataFrames[[i]]$V4 <- paste0(name,"#",listOfDataFrames[[i]]$V4)
#}
#fragments <- do.call("rbind", listOfDataFrames) #need too much space
#saveRDS(fragments,'all_fragments20210721.rds')

#fragments <- data.frame(readr::read_tsv('/gpfs1/tangfuchou_pkuhpc/zhangshu/project/gc/gc_scATAC/gc_patient/s05.analysis/archr/fragments_lib/all_lib_F4q30_CB_merged.tsv', col_names=FALSE))
#fragments <- GRanges(
  #seqnames = fragments[,1], 
  #IRanges(fragments[,2]+1, fragments[,3]), 
  #RG = fragments[,4], 
  #N = fragments[,5]
  #)
#fragments$RG <- paste0(name,"#",fragments$RG)
#saveRDS(fragments,'all_fragments20210721_raw.gr.rds')

#select cells
#proj <- loadArchRProject(path = "./Save-Proj2/", force = TRUE, showLogo = TRUE)
#used.cells<-proj$cellNames

#message("Filtering Lowly Represented Cells...")
##tabRG <- table(fragments$RG)
##keep <- names(tabRG)[which(tabRG >= minFrags)]
##fragments <- fragments[fragments$RG %in% keep,]
#fragments <- fragments[fragments$RG %in% used.cells,]
#used.chr<-paste0(c('chr'),c(1:22,'X'))
#fragments <- fragments[fragments@seqnames %in% used.chr,]
#fragments <- sort(sortSeqlevels(fragments))

#Save
#saveRDS(fragments, 'all_fragments_20210721_used.gr.rds')


#----------------------------
# Get Inputs
#----------------------------
#blacklist <- import.bed("./hg19-blacklist.v2.bed.gz")
#windows <- makeWindows(genome = BSgenome.Hsapiens.UCSC.hg19, blacklist = blacklist)
windows1 <- readRDS('hg19_blacklist_window10e6_sliding2e6.rds')
windows2 <- readRDS('hg19_blacklist_window1e6_sliding2e5.rds')
dir <- c("/storage/shuzhang/project/COLODM_COLOHSR_scRNA_ATAC/scATAC/COLODM_HSR")
proj <- loadArchRProject(path = dir, force = TRUE, showLogo = TRUE)
used.cells <- proj$cellNames #[which(proj$patient %in% c("COLO320DM")),]
samples <- as.character(read.table('samples.txt')$V1) 
#i=1
for(i in 1:length(samples)) {
#getFragmentsFromArrow() don't collect count number
print(samples[i])
fragments <- readRDS(paste0(samples[i],"_used.gr.rds"))
fragments <- fragments[fragments$RG %in% used.cells,]
print(head(fragments))
#----------------------------
# Run CNV
#----------------------------
source('./s02.Functions_cnv_greenleaf_2019nbt.r')
#rm.chrs <- c("chrM","chrX","chrY")
#cnaObj_all <- scCNA(windows1, fragments, neighbors = 100, LFC = 1.25, FDR = 0.1, force = FALSE, remove = rm.chrs)
#saveRDS(cnaObj_all, paste0(samples[i],"_scATAC_CNV_rmchr_win10e6_fragments.rds"))
used.chrs <- c("chr6","chr8","chr13","chr16")
#cnaObj_dm1 <- scCNA_dm(windows1, fragments, neighbors = 100, LFC = 1.25, FDR = 0.1, force = FALSE, used.chr = used.chrs) #only run cnv with chr in ecDNA
#saveRDS(cnaObj_dm1, paste0(samples[i],"_scATAC_CNV_usedchr_win10e6_fragments.rds"))
cnaObj_dm2 <- scCNA_dm(windows2, fragments, neighbors = 1000, LFC = 1.25, FDR = 0.1, force = TRUE, used.chr = used.chrs) #only run cnv with chr in ecDNA
saveRDS(cnaObj_dm2, paste0(samples[i],"_scATAC_CNV_usedchr_win1e6_fragments.rds"))

#----------------------------
# Plot
#----------------------------
library(RColorBrewer)
library(plotly)
library(pheatmap)


pdf(paste0(samples[i],'_scATAC_CNV_usedchr_win1e6_fragments_heat.pdf'),width=16,height=10)
#p <- cnv_heat(cnaObj=cnaObj_all,proj=proj)
#print(p)
#p1 <- cnv_heat(cnaObj=cnaObj_dm1,proj=proj)
#print(p1) #scATAC_CNV_usedchr_win10e6
p2 <- cnv_heat(cnaObj=cnaObj_dm2,proj=proj)
print(p2) #scATAC_CNV_usedchr_win1e6
dev.off()
}
