library(ArchR)
set.seed(1)
addArchRGenome("hg19")
library(readr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
##-----------------
## Reading Fragment Files
##-----------------
message("Reading in fragment files...")
dir <- c("/storage/shuzhang/project/COLODM_COLOHSR_scRNA_ATAC/scATAC/COLODM_HSR")
proj <- loadArchRProject(path = dir, force = TRUE, showLogo = TRUE)
used.cells<-proj$cellNames

input_dir <- c("/storage/shuzhang/project/COLODM_COLOHSR_scRNA_ATAC/data/")
samples <- as.character(read.table('samples.txt')$V1)
#i=1
for(i in 1:length(samples)) {
fragments <- data.frame(readr::read_tsv(paste0(input_dir,samples[i],'_atac_fragments.tsv.gz'), col_names=FALSE))

##tranfer to gr file
fragments <- GRanges(
  seqnames = paste0("chr",fragments[,1]), #check if seqnames with chr in your fragments.tsv.gz
  IRanges(fragments[,2]+1, fragments[,3]), #start and end
  RG = fragments[,4], #cell name
  N = fragments[,5] #count number
  )
fragments$RG <- paste0(samples[i],"#",fragments$RG) #change cell name to ArchR name pattern

##select cells
message("Filtering Lowly Represented Cells...")
minFrags=1e4
tabRG <- table(fragments$RG) #how many fragments in each cell barcode
keep <- names(tabRG)[which(tabRG >= minFrags)]
fragments <- fragments[fragments$RG %in% keep,]
fragments <- fragments[fragments$RG %in% used.cells,]
used.chr<- paste0(c('chr'),c(1:22,'X','Y','M')) # make sure if seqnames with 'chr'
fragments <- fragments[fragments@seqnames %in% used.chr,]

fragments <- sort(sortSeqlevels(fragments))

##Save
saveRDS(fragments, paste0(samples[i],'_used.gr.rds'))
}
