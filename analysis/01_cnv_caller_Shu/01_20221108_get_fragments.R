library(ArchR)
set.seed(1)
addArchRGenome("hg38")
library(readr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

##-----------------
## Reading Fragment Files
##-----------------

message("Reading in fragment files...")
input_dir <- c("/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_scMultiome_ColoDM-ColoHSR/COLO320HSR_5k/02_cellranger_output_King/COLO320HSR_5k_rep1_hg38/")
save_dir <- c("/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_scMultiome_ColoDM-ColoHSR/COLO320HSR_5k/03_analysis/")
sample <- c("COLO320HSR_5k_rep1")
# for whole genome: chr = c(1:22,'X','Y','M')
chr = c(8)
chr_save = c("_chr8")

fragments <- data.frame(readr::read_tsv(paste0(input_dir,'atac_fragments.tsv.gz'), col_names=FALSE))
  
# transfer to gr file
fragments <- GRanges(
  # check if seqnames with chr in your fragments.tsv.gz
  # seqnames = paste0("chr",fragments[,1]), 
  seqnames = (fragments[,1]), 
  # start and end
  IRanges(fragments[,2]+1, fragments[,3]), 
  # cell name
  RG = fragments[,4], 
  # count number
  N = fragments[,5] 
)

# change cell name to ArchR name pattern
fragments$RG <- paste0(sample,"#",fragments$RG) 

# select cells
message("Filtering Lowly Represented Cells...")
minFrags=1e4
# how many fragments in each cell barcode
tabRG <- table(fragments$RG) 
keep <- names(tabRG)[which(tabRG >= minFrags)]
fragments <- fragments[fragments$RG %in% keep,]
# make sure if seqnames with 'chr'
used.chr<- paste0(c('chr'), chr) 
fragments <- fragments[fragments@seqnames %in% used.chr,]

fragments <- sort(sortSeqlevels(fragments))

##Save
saveRDS(fragments, paste0(save_dir, sample, chr_save, '.gr.rds'))

