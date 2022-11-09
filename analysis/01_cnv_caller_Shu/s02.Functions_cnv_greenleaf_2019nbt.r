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
# scATAc CNV functions
#----------------------------
"%ni%" <- Negate("%in%") #not in

countInsertions <- function(query, fragments, by = "RG"){
  #Count By Fragments Insertions
  inserts <- c(
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(start(fragments), start(fragments)), RG = mcols(fragments)[,by]),
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(end(fragments), end(fragments)), RG = mcols(fragments)[,by])
  )
  by <- "RG"
  overlapDF <- DataFrame(findOverlaps(query, inserts, ignore.strand = TRUE, maxgap=-1L, minoverlap=0L, type = "any"))
  overlapDF$name <- mcols(inserts)[overlapDF[, 2], by]
  overlapTDF <- transform(overlapDF, id = match(name, unique(name)))
  #Calculate Overlap Stats
  inPeaks <- table(overlapDF$name)
  total <- table(mcols(inserts)[, by])
  total <- total[names(inPeaks)]
  frip <- inPeaks / total
  #Summarize
  sparseM <- Matrix::sparseMatrix(
    i = overlapTDF[, 1], 
    j = overlapTDF[, 4],
    x = rep(1, nrow(overlapTDF)), 
    dims = c(length(query), length(unique(overlapDF$name))))
  colnames(sparseM) <- unique(overlapDF$name)
  total <- total[colnames(sparseM)]
  frip <- frip[colnames(sparseM)]
  out <- list(counts = sparseM, frip = frip, total = total)
  return(out)
}

##For CNV, default windowSize = 10e6, slidingSize = 2e6
##For ecDNA, windowSize = 1e6, slidingSize = 2e5
window.Size = 10e6; sliding.Size = 2e6
makeWindows <- function(genome, blacklist, windowSize = window.Size, slidingSize = sliding.Size){
	chromSizes <- GRanges(names(seqlengths(genome)), IRanges(1, seqlengths(genome)))
	chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
	windows <- slidingWindows(x = chromSizes, width = windowSize, step = slidingSize) %>% unlist %>% .[which(width(.)==windowSize),]
	mcols(windows)$wSeq <- as.character(seqnames(windows))
  	mcols(windows)$wStart <- start(windows)
  	mcols(windows)$wEnd <- end(windows)
	message("Subtracting Blacklist...")
	windowsBL <- lapply(seq_along(windows), function(x){
			if(x %% 100 == 0){
				message(sprintf("%s of %s", x, length(windows)))
			}
			gr <- GenomicRanges::setdiff(windows[x,], blacklist)
			mcols(gr) <- mcols(windows[x,])
			return(gr)
		})
	names(windowsBL) <- paste0("w",seq_along(windowsBL))
	windowsBL <- unlist(GRangesList(windowsBL), use.names = TRUE)
	mcols(windowsBL)$name <- names(windowsBL)
	message("Adding Nucleotide Information...")
	windowSplit <- split(windowsBL, as.character(seqnames(windowsBL)))
	windowNuc <- lapply(seq_along(windowSplit), function(x){
		message(sprintf("%s of %s", x, length(windowSplit)))
	    chrSeq <- Biostrings::getSeq(genome,chromSizes[which(seqnames(chromSizes)==names(windowSplit)[x])])
	    grx <- windowSplit[[x]]
	    aFreq <- alphabetFrequency(Biostrings::Views(chrSeq[[1]], ranges(grx)))
	    mcols(grx)$GC <- rowSums(aFreq[, c("G","C")]) / rowSums(aFreq)
	    mcols(grx)$AT <- rowSums(aFreq[, c("A","T")]) / rowSums(aFreq)
	    return(grx)
	  }) %>% GRangesList %>% unlist %>% sortSeqlevels %>% sort
	windowNuc$N <- 1 - (windowNuc$GC + windowNuc$AT)
	windowNuc
}

scCNA <- function(windows, fragments, neighbors = 100, LFC = 1.5, FDR = 0.1, force = FALSE, remove = c("chrM","chrX","chrY")){
	#Keep only regions in filtered chromosomes
	windows   <- GenomeInfoDb::keepStandardChromosomes(windows, pruning.mode = "coarse")
	fragments <- GenomeInfoDb::keepStandardChromosomes(fragments, pruning.mode = "coarse")
	windows <- windows[seqnames(windows) %ni% remove]
	fragments <- fragments[seqnames(fragments) %ni% remove]

	#Count Insertions in windows
	message("Getting Counts...")
	counts <- countInsertions(windows, fragments, by = "RG")[[1]]
	message("Summarizing...")
	windowSummary <- GRangesList()
	countSummary <- matrix(nrow=length(unique(windows$name)), ncol = ncol(counts))
	for(x in seq_along(unique(mcols(windows)$name))){
		if(x %% 100 == 0){
			message(sprintf("%s of %s", x, length(unique(mcols(windows)$name))))
		}
		idx <- which(mcols(windows)$name == unique(mcols(windows)$name)[x])
		wx <- windows[idx,]
		wo <- GRanges(mcols(wx)$wSeq , ranges = IRanges(mcols(wx)$wStart, mcols(wx)$wEnd))[1,]
		mcols(wo)$name <- mcols(wx)$name[1]
		mcols(wo)$effectiveLength <- sum(width(wx))
		mcols(wo)$percentEffectiveLength <- 100*sum(width(wx))/width(wo)
		mcols(wo)$GC <- sum(mcols(wx)$GC * width(wx))/width(wo)
		mcols(wo)$AT <- sum(mcols(wx)$AT * width(wx))/width(wo)
		mcols(wo)$N <- sum(mcols(wx)$N * width(wx))/width(wo)
		countSummary[x,] <- Matrix::colSums(counts[idx,,drop=FALSE])
		windowSummary[[x]] <- wo
	}
	windowSummary <- unlist(windowSummary)
	
	#Keep only regions with less than 0.1% N
	keep <- which(windowSummary$N < 0.001) 
	windowSummary <- windowSummary[keep,]
	countSummary <- countSummary[keep,]
	
	#Now determine the nearest neighbors by GC content
	message("Computing Background...")
	bdgMean <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
	bdgSd <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
	log2FC <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
	z <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
	pval <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))

	for(x in seq_len(nrow(countSummary))){
		if(x %% 100 == 0){
			message(sprintf("%s of %s", x, nrow(countSummary)))
		}
		#Get Nearest Indices
		idxNN <- head(order(abs(windowSummary$GC[x] - windowSummary$GC)), neighbors + 1)
		idxNN <- idxNN[idxNN %ni% x]
		#Background
		if(any(colMeans(countSummary[idxNN, ])==0)){
			if(force){
				message("Warning! Background Mean = 0 Try a higher neighbor count or remove cells with 0 in colMins")
			}else{
				stop("Background Mean = 0!")
			}
		}
		bdgMean[x, ] <- colMeans(countSummary[idxNN, ])
		bdgSd[x, ] <- matrixStats::colSds(countSummary[idxNN, ])
		log2FC[x, ] <- log2((countSummary[x, ]+1e-5) / (bdgMean[x, ]+1e-5))
		z[x, ] <- (countSummary[x,] - bdgMean[x, ]) / bdgSd[x, ]
		pval[x, ] <- 2*pnorm(-abs(z[x, ]))
	}
	padj <- apply(pval, 2, function(x) p.adjust(x, method = "fdr"))
	CNA <- matrix(0, nrow=nrow(countSummary), ncol=ncol(countSummary))
	CNA[which(log2FC >= LFC & padj <= FDR)] <- 1

	se <- SummarizedExperiment(
		assays = SimpleList(
				CNA = CNA,
				counts = countSummary,
				log2FC = log2FC,
				padj = padj,
				pval = pval,
				z = z,
				bdgMean = bdgMean,
				bdgSd = bdgSd
			),
		rowRanges = windowSummary
	)
	colnames(se) <- colnames(counts)

	return(se)
}


scCNA_dm <- function(windows, fragments, neighbors = 100, LFC = 1.5, FDR = 0.1, force = TRUE, used.chr = c("chr8")){
        #Keep only regions in filtered chromosomes
        windows   <- GenomeInfoDb::keepStandardChromosomes(windows, pruning.mode = "coarse")
        fragments <- GenomeInfoDb::keepStandardChromosomes(fragments, pruning.mode = "coarse")
        windows <- windows[seqnames(windows) %in% used.chr]
        fragments <- fragments[seqnames(fragments) %in% used.chr]

        #Count Insertions in windows
        message("Getting Counts...")
        counts <- countInsertions(windows, fragments, by = "RG")[[1]]
        message("Summarizing...")
        windowSummary <- GRangesList()
        countSummary <- matrix(nrow=length(unique(windows$name)), ncol = ncol(counts))
        for(x in seq_along(unique(mcols(windows)$name))){
                if(x %% 100 == 0){
                        message(sprintf("%s of %s", x, length(unique(mcols(windows)$name))))
                }
                idx <- which(mcols(windows)$name == unique(mcols(windows)$name)[x])
                wx <- windows[idx,]
                wo <- GRanges(mcols(wx)$wSeq , ranges = IRanges(mcols(wx)$wStart, mcols(wx)$wEnd))[1,]
                mcols(wo)$name <- mcols(wx)$name[1]
                mcols(wo)$effectiveLength <- sum(width(wx))
                mcols(wo)$percentEffectiveLength <- 100*sum(width(wx))/width(wo)
                mcols(wo)$GC <- sum(mcols(wx)$GC * width(wx))/width(wo)
                mcols(wo)$AT <- sum(mcols(wx)$AT * width(wx))/width(wo)
                mcols(wo)$N <- sum(mcols(wx)$N * width(wx))/width(wo)
                countSummary[x,] <- Matrix::colSums(counts[idx,,drop=FALSE])
                windowSummary[[x]] <- wo
        }
        windowSummary <- unlist(windowSummary)

        #Keep only regions with less than 0.1% N
        keep <- which(windowSummary$N < 0.001)
        windowSummary <- windowSummary[keep,]
        countSummary <- countSummary[keep,]
        #Now determine the nearest neighbors by GC content
        message("Computing Background...")
        bdgMean <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
        bdgSd <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
        log2FC <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
        z <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
        pval <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))

        for(x in seq_len(nrow(countSummary))){
                if(x %% 100 == 0){
                        message(sprintf("%s of %s", x, nrow(countSummary)))
                }
                #Get Nearest Indices
                idxNN <- head(order(abs(windowSummary$GC[x] - windowSummary$GC)), neighbors + 1)
                idxNN <- idxNN[idxNN %ni% x]
                #Background
                if(any(colMeans(countSummary[idxNN, ])==0)){
                        if(force){
                                message("Warning! Background Mean = 0 Try a higher neighbor count or remove cells with 0 in colMins")
                        }else{
                                stop("Background Mean = 0!")
                        }
                }
                bdgMean[x, ] <- colMeans(countSummary[idxNN, ])
                bdgSd[x, ] <- matrixStats::colSds(countSummary[idxNN, ])
                log2FC[x, ] <- log2((countSummary[x, ]+1e-5) / (bdgMean[x, ]+1e-5))
                z[x, ] <- (countSummary[x,] - bdgMean[x, ]) / bdgSd[x, ]
                pval[x, ] <- 2*pnorm(-abs(z[x, ]))
        }
        padj <- apply(pval, 2, function(x) p.adjust(x, method = "fdr"))
        CNA <- matrix(0, nrow=nrow(countSummary), ncol=ncol(countSummary))
        CNA[which(log2FC >= LFC & padj <= FDR)] <- 1

        se <- SummarizedExperiment(
                assays = SimpleList(
                                CNA = CNA,
                                counts = countSummary,
                                log2FC = log2FC,
                                padj = padj,
                                pval = pval,
                                z = z,
                                bdgMean = bdgMean,
                                bdgSd = bdgSd
                        ),
                rowRanges = windowSummary
        )
        colnames(se) <- colnames(counts)

        return(se)
}

library(RColorBrewer)
library(plotly)
library(pheatmap)
library(rtracklayer)

cnv_heat <- function(cnaObj=cnaObj_dm1, proj=proj) {
data<-t(cnaObj@assays@data@listData$log2FC)
rownames(data) <- rownames(cnaObj@colData)
chr_window_num <- as.data.frame(t(table(seqnames(cnaObj@rowRanges))))
chr_window_num <- chr_window_num[,2:3]; colnames(chr_window_num) <- c("Chr","Number")
#chr_window_num <- chr_window_num[order(factor(chr_window_num$Chr,levels = paste0(c('chr'),c(1:22,"X")),ordered=T)),]
for(i in 2:nrow(chr_window_num)){
  chr_window_num[i,2] <-chr_window_num[i,2]+chr_window_num[i-1,2]
}

bed_file <- c("COLODM_ecDNA.bed")
acd.gr <- import(bed_file, genome = "hg19")

peaks_Inacd <- findOverlaps(acd.gr, cnaObj@rowRanges)
annotation_col = data.frame(chr=cnaObj@rowRanges@seqnames,ecDNA="N")
annotation_col[unique(subjectHits(peaks_Inacd)),]$ecDNA <- c("Y")

cell.info<-data.frame(cell=proj$cellNames,cluster=proj$Clusters,patient=proj$patient)
rownames(cell.info) <- cell.info$cell
annotation_row = cell.info[intersect(rownames(cell.info),rownames(data)),c('cluster','patient')]
annotation_row <- annotation_row[with(annotation_row,order(cluster)),]

data<-data[rownames(annotation_row),]
#print(dim(data))
#gc.colorp <- c(brewer.pal(9,"Set1"),brewer.pal(8,"Dark2"),brewer.pal(8,"Paired"))
##for epi cells
#gc.colorc <- c(brewer.pal(10,"Set3"),brewer.pal(3,"Set2"))
#gc.patient <- unique(annotation_row$patient)
#gc.cluster <- unique(annotation_row$cluster)
#ann_colors = list(patient=setNames(gc.colorp,gc.patient),cluster=setNames(gc.colorc,gc.cluster))
bk = unique(c(seq(-5,5, length=100)))
p <- pheatmap(data,col=colorRampPalette(c("navy","white","red"))(100),breaks = bk,
             cluster_rows =F, cluster_cols = F,scale = "none", #gaps_col = chr_window_num$Number,
             annotation_row = annotation_row, annotation_col = annotation_col, #annotation_colors = ann_colors,#
            legend = TRUE,show_rownames = F, show_colnames = F,
            border_color=NA,fontsize = 9,fontsize_row = 10,fontsize_col = 10)
}
