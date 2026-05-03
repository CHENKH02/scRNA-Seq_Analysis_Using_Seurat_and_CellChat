#----------------------------------
## 1. Load the HDF5 count matrix:
#----------------------------------

#Load Seurat package for scRNA-seq analysis
library(Seurat)

#Read 10X Genomics HDF5 count matrix (genes × cells)
sc_counts = Read10X_h5("5k_Human_Donor1_PBMC_3p_gem-x_5k_Human_Donor1_PBMC_3p_gem-x_count_sample_filtered_feature_bc_matrix.h5")


#----------------------------------
## 2. Create the Seurat object:
#----------------------------------

#Initialize Seurat object, keep genes expressed in at least 3 cells, keep cells with at least 200 detected genes
seurat_obj = CreateSeuratObject(
  counts = sc_counts,
  min.cells = 3,
  min.features = 200
)

#Inspect object structure and dimensions
seurat_obj
dim(seurat_obj)

#Check metadata that contains per-cell QC metrics
head(seurat_obj@meta.data)


#-------------------------------------
## 3. Load cell type annotation table:
#-------------------------------------

#Load external annotation table
cell_type_annotation = read.csv("cell_types.csv")

#Inspect table structure
head(cell_type_annotation)
colnames(cell_type_annotation)


#---------------------------------------------------------
## 4. Quality control(calculate mitochondrial percentage):
#---------------------------------------------------------

#Calculate percentage of mitochondrial genes (indicator of cell stress)
seurat_obj[["percent.mt"]] = PercentageFeatureSet(
  seurat_obj,
  pattern = "^MT-"
)

#Load plotting libraries
library(ggplot2)
library(tidyr)

#Extract metadata for visualization
meta <- seurat_obj@meta.data

#QC metric 1: number of detected genes per cell
p1 <- ggplot(meta, aes(x = "", y = nFeature_RNA)) +
  geom_violin(fill = "lightblue") +
  theme_classic() +
  labs(title = "nFeature_RNA", x = "", y = "Value")

#QC metric 2: total counts per cell
p2 <- ggplot(meta, aes(x = "", y = nCount_RNA)) +
  geom_violin(fill = "lightgreen") +
  theme_classic() +
  labs(title = "nCount_RNA", x = "", y = "Value")

#QC metric 3: mitochondrial percentage
p3 <- ggplot(meta, aes(x = "", y = percent.mt)) +
  geom_violin(fill = "navy") +
  theme_classic() +
  labs(title = "percent.mt", x = "", y = "Value")

#Combine the plots
library(patchwork)
qc_plot <- p1 + p2 + p3
qc_plot

#Save QC plots
ggsave("qc.png", plot = qc_plot, width = 9, height = 4)

#Explore QC relationships,Check relationships between QC metrics
FeatureScatter(
  seurat_obj,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

FeatureScatter(
  seurat_obj,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

FeatureScatter(
  seurat_obj,
  feature1 = "nFeature_RNA",
  feature2 = "percent.mt"
)


#------------------
## 5. filter cells:
#------------------

#Remove low gene count cells (<200), potential doublets (>2500 genes), high mitochondrial cells (>5%)
seurat_obj = subset(
  seurat_obj,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 2500 &
    percent.mt < 5
)

#Check remaining cells
seurat_obj
dim(seurat_obj)


#------------------------
## 6. Normalize the data:
#------------------------

# Normalize counts per cell using 'LogNormalize'
seurat_obj = NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


#-------------------------------------
## 7. Identify 2,000 highly variable genes:
#-------------------------------------

#Select top 2000 highly variable genes that capture biological signal
seurat_obj = FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000
)

#Inspect top variable genes
top10 <- head(VariableFeatures(seurat_obj), 10)
top10


#----------------------------------------------
## 8. Dimensionality reduction (Scaling + PCA):
#----------------------------------------------

#Scale data (center and standardize) and regress out mitochondrial effect
all_genes = rownames(seurat_obj)

seurat_obj = ScaleData(
  seurat_obj,
  features = all_genes,
  vars.to.regress = "percent.mt"
)

#Run PCA on variable genes
seurat_obj = RunPCA(
  seurat_obj,
  features = VariableFeatures(seurat_obj),
  npcs = 50
)

#Visualize PCA loadings (top genes in PC1 and PC2)
VizDimLoadings(
  seurat_obj,
  dims = 1:2,
  reduction = "pca"
)

#PCA plot of cells
DimPlot(
  seurat_obj,
  reduction = "pca"
)

#Inspect first 5 PCs using Heatmap
DimHeatmap(
  seurat_obj,
  dims = 1:5,
  cells = 500,
  balanced = TRUE
)

#Use ElbowPlot() to determine the number of PCs
ElbowPlot(seurat_obj, ndims = 50)


#-----------------------
## 9. Clustering + UMAP:
#-----------------------

#Construct nearest neighbor graph
seurat_obj = FindNeighbors(
  seurat_obj,
  dims = 1:10
)

#Cluster cells
seurat_obj = FindClusters(
  seurat_obj,
  resolution = 0.5
)

#Check cluster sizes
table(Idents(seurat_obj))

#UMAP visualization
seurat_obj = RunUMAP(
  seurat_obj,
  dims = 1:10
)

p_umap <- DimPlot(
  seurat_obj,
  reduction = "umap",
  label = TRUE
)

#Save UMAP
ggsave("umap.png", plot = p_umap, width = 6, height = 5)

p_umap

#Resolution testing - 0.2(lower resolution, fewer clusters)
seurat_obj = FindClusters(
  seurat_obj,
  resolution = 0.2
)

table(Idents(seurat_obj))

#Resolution testing - 1.0(higher resolution, more clusters)
seurat_obj = FindClusters(
  seurat_obj,
  resolution = 1.0
)

table(Idents(seurat_obj))


#---------------------------------
## 10. Marker gene identification:
#---------------------------------

#Identify cluster-specific marker genes
markers = FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.25
)

#Select top 2 markers per cluster
top2_markers <- markers[order(markers$cluster, -markers$avg_log2FC), ]

top2_markers <- do.call(rbind, lapply(split(top2_markers, top2_markers$cluster), head, 2))


#Heatmap of marker genes
p_heatmap <- DoHeatmap(
  seurat_obj,
  features = top2_markers$gene
) 

#Save Heatmap
ggsave("heatmap.png", plot = p_heatmap, width = 20, height = 6)

p_heatmap


#--------------------------------
## 11. Add cell type annotations:
#--------------------------------

#Check data
head(colnames(seurat_obj))
head(cell_type_annotation$barcode)

#Match barcodes to annotation table
cell_type_vector = cell_type_annotation$coarse_cell_type[
  match(colnames(seurat_obj), cell_type_annotation$barcode)
]

#Add to metadata
seurat_obj = AddMetaData(
  seurat_obj,
  metadata = cell_type_vector,
  col.name = "cell_type"
)

#Check annotation
head(seurat_obj@meta.data)
table(seurat_obj$cell_type)

#Set identity to cell type
Idents(seurat_obj) = "cell_type"

#Annotated UMAP
DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE
)

#Compare clusters vs cell types
seurat_obj$seurat_cluster = seurat_obj$seurat_clusters

cluster_celltype_table = table(
  seurat_obj$seurat_cluster,
  seurat_obj$cell_type
)

cluster_celltype_prop = prop.table(
  cluster_celltype_table,
  margin = 1
)

#Create dataframe
cluster_celltype_df = as.data.frame(cluster_celltype_prop)

colnames(cluster_celltype_df) = c(
  "seurat_cluster",
  "cell_type",
  "proportion"
)

#Barplot
ggplot(
  cluster_celltype_df,
  aes(x = seurat_cluster, y = proportion, fill = cell_type)
) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(
    x = "Seurat cluster",
    y = "Proportion",
    fill = "Cell type"
  )

table(seurat_obj$cell_type)

#--------------------------------------------------------------------
## 12. Pathway enrichment(ORA) for cell type marker genes using KEGG:
#--------------------------------------------------------------------

#Use cell types as identity
Idents(seurat_obj) = "cell_type"

#Marker genes per cell type
celltype_markers = FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.25
)

#Load KEGG pathways
library(msigdbr)

#KEGG gene sets
kegg = msigdbr(
  species = "Homo sapiens",
  collection = "C2",
  subcollection = "CP:KEGG_LEGACY"
)

kegg = kegg[, c("gs_name", "gene_symbol")]

#ORA analysis per cell type
library(clusterProfiler)

ora_results = data.frame()

for (ct in unique(celltype_markers$cluster)) {
  
  sub_markers = subset(
    celltype_markers,
    cluster == ct
  )
  
  gene_list = unique(sub_markers$gene)
  
  ora = enricher(
    gene_list,
    TERM2GENE = kegg,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
  
  if (!is.null(ora)) {
    
    ora_df = as.data.frame(ora)
    
    if (nrow(ora_df) > 0) {
      ora_df$cell_type = ct
      ora_results = rbind(ora_results, ora_df)
    }
  }
}

#Keep significant pathways
ora_results = subset(
  ora_results,
  p.adjust <= 0.1
)

#Compute RichFactor
ora_results$RichFactor = sapply(
  strsplit(ora_results$GeneRatio, "/"),
  function(x) as.numeric(x[1]) / as.numeric(x[2])
)

#Clean pathway names
ora_results$Description = gsub("KEGG_", "", ora_results$Description)
ora_results$Description = gsub("_", " ", ora_results$Description)


#Dot plot for ORA
p_ora <- ggplot(
  ora_results,
  aes(
    x = cell_type,
    y = Description
  )
) +
  geom_point(
    aes(
      size = RichFactor,
      fill = p.adjust
    ),
    shape = 21,
    color = "black"
  ) +
  theme_bw() +
  labs(
    x = "Cell type",
    y = "KEGG pathway",
    size = "RichFactor",
    fill = "Adjusted p-value"
  )

ggsave("ora.png", plot = p_ora, width = 7, height = 5)

p_ora

#--------------------------------------------------
## 13. Cell–cell communication analysis (CellChat):
#--------------------------------------------------

#Create CellChat object
library(CellChat)

table(seurat_obj$cell_type)

cellchat = createCellChat(
  object = seurat_obj,
  group.by = "cell_type"
)

#Use human ligand-receptor database
CellChatDB = CellChatDB.human
cellchat@DB = CellChatDB

#Preprocess CellChat data
cellchat = subsetData(cellchat)
cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat)

#Compute communication probabilities
cellchat = computeCommunProb(cellchat)

cellchat = filterCommunication(
  cellchat,
  min.cells = 5
)

#Pathway-level probabilities analysis
cellchat = computeCommunProbPathway(cellchat)

#Aggregate network
cellchat = aggregateNet(cellchat)

#Circle plot for number of interactions
png("cellchat_count.png", width = 800, height = 800)

groupSize = as.numeric(table(cellchat@idents))

netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of interactions"
)

dev.off()

#Circle plot for interaction strength
png("cellchat_weight.png", width = 800, height = 800)

netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction strength"
)

dev.off()

#Sender vs receiver roles
png("cellchat_roles.png", width = 800, height = 800)

cellchat = netAnalysis_computeCentrality(
  cellchat,
  slot.name = "netP"
)

netAnalysis_signalingRole_scatter(cellchat)

dev.off()




