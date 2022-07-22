library(Signac)
library(Seurat)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
set.seed(1234)



load('./filtered_SNARE.RData')



# =====================================================================
#                     import  my  data and filter data
# =====================================================================
SCARP_ATAC_Cells_df <-
  read.table(
    './SCARP_ATAC_Cells_df.txt',
    sep = ',',
    row.names = 1,
    header = T
  )

SCARP_ATAC_Peaks_df <-
  read.table(
    './SCARP_ATAC_Peaks_df.txt',
    sep = ',',
    row.names = 1,
    header = T
  )

common_cells <- intersect(rownames(SCARP_ATAC_Cells_df),
                          rownames(snare@meta.data))
atac <- atac[rownames(SCARP_ATAC_Peaks_df),
             common_cells]
rna <- rna[snare@assays[["RNA"]]@data@Dimnames[[1]],
           common_cells]


# create a Seurat object and add the assays
snare <- CreateSeuratObject(counts = rna)

snare[['ATAC']] <- CreateChromatinAssay(
  counts = atac,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = fragments
)

# add the gene information to the object
Annotation(snare[["ATAC"]]) <- annotations



# =================================================================
# add cell type annotaion
SNARE_cell_type <- read.csv('./SNARE_Celltype.csv', row.names = 1)
help_id <-
  data.frame(1:length(unique(SNARE_cell_type$SNARE_Celltype)))
rownames(help_id) <- unique(SNARE_cell_type$SNARE_Celltype)
temp_name <- c()
for (i in 1:dim(SNARE_cell_type)[1]) {
  temp_name <-
    c(temp_name, help_id[SNARE_cell_type$SNARE_Celltype[i], ])
}
SNARE_cell_type$id <- temp_name
rm(temp_name, i, help_id)


snare$standard_Seurat_lablled <-
  SNARE_cell_type[rownames(snare@meta.data),]$SNARE_Celltype



load('./filter_SNARE_celltype.RData')




# ===========================================================================
#                                 SCARP
# ============================================================================
snare[["SCARP_feat"]] <- CreateDimReducObject(
  embeddings = as.matrix(SCARP_ATAC_Cells_df[common_cells, ]),
  key = "scarp_",
  assay = 'ATAC'
)

choose_dim <- dim(SCARP_ATAC_Cells_df)[2]

snare <- FindNeighbors(snare,
                       reduction = "SCARP_feat",
                       dims = 1:choose_dim, )

snare <- RunUMAP(
  snare,
  reduction = 'SCARP_feat',
  dims = 1:choose_dim,
  reduction.name = "umap.atac.scarp"
)


DimPlot(snare,
        reduction = 'umap.atac.scarp',
        label = F,
        group.by = 'standard_Seurat_lablled') +
  ggtitle('SCARP ATAC')







# ===========================================================================
#                                 cisTOPIC
# ============================================================================
cisTOPIC_cell_embedding <-
  read.table('./Data/cisTOPIC_cell_embedding.csv', sep = '\t')
cisTOPIC_cell_embedding <- t(cisTOPIC_cell_embedding)
rownames(cisTOPIC_cell_embedding) <-
  snare@assays[["ATAC"]]@data@Dimnames[[2]]


snare[["cisTOPIC_feat"]] <- CreateDimReducObject(
  embeddings = as.matrix(cisTOPIC_cell_embedding),
  key = "cistopic_",
  assay = 'ATAC'
)

choose_dim <- dim(cisTOPIC_cell_embedding)[2]

snare <- FindNeighbors(snare,
                       reduction = "cisTOPIC_feat",
                       dims = 1:choose_dim)

snare <- RunUMAP(
  snare,
  reduction = 'cisTOPIC_feat',
  dims = 1:choose_dim,
  reduction.name = "umap.atac.cistopic"
)


DimPlot(snare,
        label = F,
        group.by = 'standard_Seurat_lablled',
        reduction = 'umap.atac.cistopic',) +
  labs(title = 'cisTOPIC ATAC')
# save: 600*500





# ===========================================================================
#                                 scAND
# ============================================================================
scAND_cell_embedding <-
  read.csv('./Data/scAND_cell_embedding.csv', row.names = 1)

rownames(scAND_cell_embedding) <-
  snare@assays[["ATAC"]]@data@Dimnames[[2]]


snare[["scAND_feat"]] <- CreateDimReducObject(
  embeddings = as.matrix(scAND_cell_embedding),
  key = "scand_",
  assay = 'ATAC'
)

choose_dim <- dim(scAND_cell_embedding)[2]

snare <- FindNeighbors(snare,
                       reduction = "scAND_feat",
                       dims = 1:choose_dim)

snare <- RunUMAP(
  snare,
  reduction = 'scAND_feat',
  dims = 1:choose_dim,
  reduction.name = "umap.atac.scand"
)


DimPlot(snare,
        label = F,
        group.by = 'standard_Seurat_lablled',
        reduction = 'umap.atac.scand',) +
  labs(title = 'scAND ATAC')
# save: 600*500





# ===========================================================================
#                                 RNA
# ============================================================================


DefaultAssay(snare) <- "RNA"

snare <- FindVariableFeatures(snare, nfeatures = 3000)
snare <- NormalizeData(snare)
snare <- ScaleData(snare)
snare <- RunPCA(snare, npcs = 30)
snare <- RunUMAP(snare, dims = 1:30, reduction.name = "umap.rna")
snare <- FindNeighbors(snare, dims = 1:30)
snare <- FindClusters(snare, resolution = 0.5, algorithm = 3)
DimPlot(snare,
        label = T,
        group.by = 'standard_Seurat_lablled',
        reduction = 'umap.rna',) +
  labs(title = 'Seurat RNA')



# ===================================================
#  计算各方法轮廓系
# ===================================================
library(factoextra)
library(cluster)

# rename cell type id
new.cluster.ids <- data.frame(id = 1:13)
rownames(new.cluster.ids) <- c(
  "Astro",
  "L2/3 IT",
  "L4",
  "L5 CT",
  "L5 PT",
  "L6 IT",
  "L6 IT.2",
  "L6b",
  "NP",
  "Oligo",
  "Pvalb",
  "Sst",
  "Vip/Lamp5"
)
SNARE_cell_type$id <-
  new.cluster.ids[SNARE_cell_type$SNARE_Celltype, ]


x <- SNARE_cell_type[common_cells, ]$id
names(x) <- rownames(snare@meta.data)


# use UMAP 2-dimensional
feat <- 'umap.rna'
feat <- 'umap.atac.scarp'
feat <- 'umap.atac.cistopic'
feat <- 'umap.atac.scand'

sil <-
  silhouette(x, dist(snare@reductions[[feat]]@cell.embeddings))
fviz_silhouette(sil)

# 600*400
