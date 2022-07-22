library(clusterProfiler)
library(org.Hs.eg.db)
library(topGO)
library(DOSE)
# library(doseplot)
library(patchwork)
library(enrichplot)
library(ReactomePA)
library(forcats)
library(ggstance)


# 
# file_name1 <- './Data/Top100sites_3UTR.csv'
# file_name2 <- './Data/Top100sites_promoter.csv'
# genes1 <- read.csv(file_name1)
# genes2 <- read.csv(file_name2)
# genes <- rbind(genes1, genes2)
# # 
# 
# 
# 
# CC <- enrichGO(gene = genes$ENSEMBL,  #Gene list
#                keyType = "ENSEMBL",  #Gene ID type
#                OrgDb=org.Hs.eg.db,  #Species
#                ont = "CC",   
#                pvalueCutoff = 1,  #pvalue threshold
#                pAdjustMethod = "fdr",  #Multiple hypothesis testing correction method
#                minGSSize = 1,   #Minimum set of genes to annotate
#                maxGSSize = 500,  #Maximun set of genes to annotate
#                qvalueCutoff = 1,  
#                readable = TRUE)  #Gene ID to gene name conversion
# 
# MF <- enrichGO(gene = genes$ENSEMBL,  #Gene list
#                keyType = "ENSEMBL",  #Gene ID type
#                OrgDb=org.Hs.eg.db,  #Species
#                ont = "MF",   
#                pvalueCutoff = 1,  #pvalue threshold
#                pAdjustMethod = "fdr",  #Multiple hypothesis testing correction method
#                minGSSize = 1,   #Minimum set of genes to annotate
#                maxGSSize = 500,  #Maximun set of genes to annotate
#                qvalueCutoff = 1,  
#                readable = TRUE)  #Gene ID to gene name conversion
# 
# 
# BP <- enrichGO(gene = genes$ENSEMBL,  #Gene list
#                keyType = "ENSEMBL",  #Gene ID type
#                OrgDb=org.Hs.eg.db,  #Species
#                ont = "BP",   
#                pvalueCutoff = 1,  #pvalue threshold
#                pAdjustMethod = "fdr",  #Multiple hypothesis testing correction method
#                minGSSize = 1,   #Minimum set of genes to annotate
#                maxGSSize = 500,  #Maximun set of genes to annotate
#                qvalueCutoff = 1,  
#                readable = TRUE)  #Gene ID to gene name conversion
# 
# # keytypes(org.Hs.eg.db)
# ENTREZID <- bitr(gene = genes$ENSEMBL,
#                  fromType = "ENSEMBL",
#                  toType = 'ENTREZID',
#                  OrgDb = 'org.Hs.eg.db')
# 
# 
# KEGG <- enrichKEGG(gene = ENTREZID$ENTREZID,
#                  keyType = "kegg", 
#                  organism = 'hsa',
#                  pvalueCutoff = 1,
#                  pAdjustMethod = "fdr",  
#                  qvalueCutoff = 1)
# 
# 
# barplot(CC,  #
#         x = "GeneRatio",  # x stick
#         color = "p.adjust",  #right y stick
#         showCategory = 10,  # show top 20 dots
#         title = "Cellular Components" 
# )
# barplot(MF,  #
#         x = "GeneRatio",  # x stick
#         color = "p.adjust",  #right y stick
#         showCategory = 10,  # show top 20 dots
#         title = "Molecular Function" 
# )
# barplot(BP,  #
#         x = "GeneRatio",  # x stick
#         color = "p.adjust",  #right y stick
#         showCategory = 10,  # show top 20 dots
#         title = "Biological Process"
# )
# barplot(KEGG,  #
#         x = "GeneRatio",  # x stick
#         color = "p.adjust",  #right y stick
#         showCategory = 10,  # show top 20 dots
#         title = "KEGG" 
# )
# 
# # gseaplot2(CC, genes$ENSEMBL)
# 
# # enrichMap(CC, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
# 
# # 
# cnetplot(CC, foldChange=genes$ENSEMBL,
#          showCategory = 3, colorEdge = T,
#          node_label="gene",
#          layout = "gem")
# 
# heatplot(CC)




# enrichGO_KEGG <- function(gene_list, keytype, back_gene){
#   library(clusterProfiler)
#   library(org.Hs.eg.db)
#   library(topGO)
#   library(DOSE)
#   library(patchwork)
#   library(enrichplot)
#   library(ReactomePA)
#   library(forcats)
#   library(ggstance)
#   
#   CC <- enrichGO(gene = gene_list,  
#                  keyType = keytype, 
#                  OrgDb=org.Hs.eg.db,  
#                  ont = "CC",   
#                  pvalueCutoff = 1,  
#                  pAdjustMethod = "fdr",  
#                  universe = back_gene,
#                  minGSSize = 1,   
#                  maxGSSize = 500,  
#                  qvalueCutoff = 1,  
#                  readable = TRUE)  
#   
#   MF <- enrichGO(gene = gene_list,  
#                  keyType = keytype,  
#                  OrgDb=org.Hs.eg.db, 
#                  ont = "MF",   
#                  pvalueCutoff = 1, 
#                  pAdjustMethod = "fdr",  
#                  universe = back_gene,
#                  minGSSize = 1, 
#                  maxGSSize = 500, 
#                  qvalueCutoff = 1,  
#                  readable = TRUE)  
#   
#   BP <- enrichGO(gene = gene_list,  
#                  keyType = keytype,  
#                  OrgDb=org.Hs.eg.db,  
#                  ont = "BP",   
#                  pvalueCutoff = 1,  
#                  pAdjustMethod = "fdr",  
#                  universe = back_gene,
#                  minGSSize = 1,  
#                  maxGSSize = 500, 
#                  qvalueCutoff = 1,  
#                  readable = TRUE)  
#   
#   if (keytype == "ENTREZID"){
#     KEGG <- enrichKEGG(gene = gene_list,
#                        keyType = "kegg", 
#                        organism = 'hsa',
#                        pvalueCutoff = 1,
#                        pAdjustMethod = "fdr",  
#                        universe = back_gene,
#                        qvalueCutoff = 1)
#     return(list(CC=CC, 
#                 MF=MF, 
#                 BP=BP, 
#                 KEGG=KEGG))}
#   else{
#     return(list(CC=CC, 
#                 MF=MF, 
#                 BP=BP))
#   }
#   
# }
# 








#############################################################################
#############################################################################

gene_list <-  read.csv('../cytoscape/scand default node.csv')
gene_list <- gene_list$commonName
gene_list <- bitr(gene = gene_list,
                  fromType = "SYMBOL",
                  toType = 'ENTREZID',
                  OrgDb = 'org.Hs.eg.db')


KEGG <- enrichKEGG(gene = gene_list$ENTREZID,
                   keyType = "kegg", 
                   organism = 'hsa',
                   pvalueCutoff = 1,
                   pAdjustMethod = "fdr",  
                   universe = NULL,
                   qvalueCutoff = 1)

barplot(KEGG, x = "GeneRatio", 
        color = "p.adjust",  
        showCategory = 15,  
        title = "KEGG")

library(stringr)
for (ii in 1:length(KEGG@result[["geneID"]])){
  set_ <- strsplit(KEGG@result[["geneID"]][ii],'/')
  
  for(i in 1:length(set_[[1]])){
    set_[[1]][i] <- gene_list[gene_list$ENTREZID==set_[[1]][i],]$SYMBOL
  }
  
  KEGG@result[["geneID"]][ii] <- str_c(set_[[1]],collapse='/')
}

  
cnetplot(KEGG,
         showCategory = 5, 
         colorEdge = T,
         node_label="all",
         circular = T)




