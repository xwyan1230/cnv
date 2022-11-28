set.seed(1)
addArchRGenome("hg38")
library(GenomicRanges)
library(trackViewer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringr)


# ----------------------------
# Get Inputs
# ----------------------------
input_dir <- c("/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/")
save_dir <- c("/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/")
sample <- c("COLO320DM_5k_rep1")
chr = 8
used.chr<- paste0(c('chr'), chr)

G1 = read.table(paste0(input_dir, sample, '_sorted_245_G1.txt'), sep='\t')
S = read.table(paste0(input_dir, sample, '_sorted_245_S.txt'), sep='\t')
G2M = read.table(paste0(input_dir, sample, '_sorted_245_G2M.txt'), sep='\t')

Colo320DM <- readRDS(paste0(input_dir, sample, "_chr", chr, ".gr.rds"))
Colo320DM$cell <- str_split_fixed(Colo320DM$RG, '#', 2)[,2]

Colo320DM_G1 <- Colo320DM
Colo320DM_G1 <- Colo320DM_G1[Colo320DM_G1$cell %in% G1$V1]
Colo320DM_S <- Colo320DM
Colo320DM_S <- Colo320DM_S[Colo320DM_S$cell %in% S$V1]
Colo320DM_G2M <- Colo320DM
Colo320DM_G2M <- Colo320DM_G2M[Colo320DM_G2M$cell %in% G2M$V1]

pdf(paste0(save_dir, sample,'_track_127600001-127650001.pdf'),width=16,height=10)
p1 <- plotGRanges(Colo320DM_G1, Colo320DM_S, Colo320DM_G2M, range=GRanges("chr8", IRanges(127600001,127650001)))
print(p1)
dev.off()
# addGuideLine(guideLine=c(127610001, 127620001, 127630001), vp=p1)
# addArrowMark(list(x=unit(.5, "npc"), y=unit(.39, "npc")), col="blue")