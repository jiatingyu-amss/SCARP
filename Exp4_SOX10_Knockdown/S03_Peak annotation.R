library(ChIPseeker)
library(ggplot2)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggupset)
library(ggimage)
require(clusterProfiler)
source('./gene_richment (clusterProfiler).R')
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene



peaks <- readPeakFile('./SOX10 data/Sox10KD_peaks.bed')


#=========================================================================
#             Peaks annotation (bed file + txdb object)
#=========================================================================
# 1. back
# =============================================================
back_peaks <- readPeakFile('./SOX10 data/Sox10KD_peaks.bed')
back_x_annogene <- annotatePeak(back_peaks,
                                tssRegion = c(-1000, 1000),
                                TxDb = txdb,
                                annoDb = 'org.Hs.eg.db')
back_genes <- as.data.frame(back_x_annogene)
back_genes <- back_genes[grepl('Promoter', back_genes$annotation),]
# write.csv(back_genes, './SOX10 data/annotated_genes_1000bp.csv')



# # or use gene transforamtion
# back_genes <- seq2gene(back_peaks,
#                        tssRegion = c(-1000, 1000),
#                        flankDistance = 1000,
#                        TxDb=txdb)

# 2. diff
# =============================================================
diff_peaks <- readPeakFile('./Data/diff_sites200_peaks.bed')
diff_peaks <-  diff_peaks[1:50]
# diff_peaks <- readPeakFile('./Data/promoter/diff_sites200_peaks_promoter (2).bed')
#
# diff_x_annogene <- annotatePeak(diff_peaks,
#                                 tssRegion = c(-1000, 1000),
#                                 TxDb=txdb,
#                                 annoDb = 'org.Hs.eg.db')
# diff_genes <- as.data.frame(diff_x_annogene)
# diff_genes <- diff_genes[grepl('Promoter', 
#                                diff_genes$annotation),]

# # or use gene transforamtion
diff_genes <- seq2gene(diff_peaks,
                       tssRegion = c(-1000, 1000),
                       flankDistance = 1000,
                       TxDb=txdb)


# 3. enrich_result
# =============================================================
# print(paste0('diff_genes: ', length(unique(diff_genes$geneId))))
# print(paste0('back_genes: ',length(unique(back_genes$geneId))))

print(paste0('diff_genes: ', length(unique(diff_genes))))
print(paste0('back_genes: ',length(unique(back_genes))))
enrich_result_ENTREZID <- enrichGO_KEGG(gene_list = diff_genes,
                                        keytype = 'ENTREZID',
                                        back_gene = back_genes)

# enrich_result_ENTREZID <- enrichGO_KEGG(gene_list = diff_genes$geneId,
#                                         keytype = 'ENTREZID',
#                                         back_gene = back_genes$geneId)

# enrich_result_ENSEMBL<- enrichGO_KEGG(gene_list = diff_genes$ENSEMBL,
#                                      keytype ='ENSEMBL',
#                                      back_gene = back_genes$ENSEMB)



#=========================================================================
#             Data Visualization 
#=========================================================================
# covplot(diff_peaks)
# plotAnnoPie(back_x_annogene)
# upsetplot(back_x_annogene, vennpie=TRUE)
# plotDistToTSS(back_x_annogene,
#               title = 'Distribution of ATAC-seq loci relative to TSS',
#               ylab = 'Peaks (%) (5\'-> 3\')')
# 
# # only keep promoter sites
# promoter_anno <- subset(back_x_annogene, abs(distanceToTSS)<1000)
# promoter_anno <- subset(promoter_anno, grepl('Promoter', annotation))
# y_promoter_anno <- as.data.frame(promoter_anno)
# 
# plotDistToTSS(promoter_anno,
#               title = 'Distribution of ATAC-seq loci relative to TSS',
#               ylab = 'Peaks (%) (5\'-> 3\')')
# upsetplot(promoter_anno, vennpie=TRUE)
# plotAnnoPie(promoter_anno)



barplot(enrich_result_ENTREZID$CC, x = "GeneRatio", 
        color = "p.adjust",  
        showCategory = 15,  
        title = "Cellular Components")

barplot(enrich_result_ENTREZID$MF, x = "GeneRatio", 
        color = "p.adjust",  
        showCategory = 15,  
        title = "Molecular Function")

barplot(enrich_result_ENTREZID$BP, x = "GeneRatio", 
        color = "p.adjust",  
        showCategory = 15,  
        title = "Biological Process")

barplot(enrich_result_ENTREZID$KEGG, x = "GeneRatio", 
        color = "p.adjust",  
        showCategory = 15,  
        title = "KEGG")

cnetplot(enrich_result_ENTREZID$KEGG,
         foldChange = diff_genes$SYMBOL,
         showCategory = 15, colorEdge = T,
         node_label="gene",
         layout = "gem")

# heatplot(enrich_result$CC)



