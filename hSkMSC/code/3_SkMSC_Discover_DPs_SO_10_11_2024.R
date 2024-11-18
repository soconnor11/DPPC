# docker run -it -v '/media:/files' cplaisier/ccafv2_seurat5_presto

#--------------------------------
# Set up section / load packages
#--------------------------------

#remotes::install_version("matrixStats", version="1.1.0")
#library(matrixStats)
library(dplyr)
library(Seurat)
library(SeuratDisk)
library(keras)
library(patchwork)
library(ggplot2)
library(grid)
library(gridExtra)
library(writexl)
library(data.table)
library(readr)
library("org.Hs.eg.db")
#library(aricode)
#library(reticulate)
use_python('/usr/bin/python3')
library(ccAFv2)
#devtools::install_github('immunogenomics/presto')
library(presto)
library(pheatmap)

# Set working directory
setwd("/files")
mgenes = read.csv(system.file('extdata', 'ccAFv2_genes.csv', package='ccAFv2'), header=TRUE, row.names=1)[,paste0('human_ensembl')]

#install.packages("Matrix", type = "source")
#install.packages("irlba", type = "source")

# Some features to investigate
features1 = c('S100B', 'SOX2', 'SOX4', 'MKI67', 'APOE', 'VIM', 'CLU', 'FABP7','OLIG1','OLIG2', 'DLL3', 'HES6')
# convert to ensembl IDs
ensembl_features1 = mapIds(org.Hs.eg.db, keys = features1, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
ensembl_features1 = na.omit(data.frame(ensembl_features1))
ensembl_features1_plot = ensembl_features1$ensembl_features1

features2 = c("MBP", "PLP1", "ETNPPL", "CD14","CX3CR1","PTPRC", "RBFOX3")
# convert to ensembl IDs
ensembl_features2 = mapIds(org.Hs.eg.db, keys = features2, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
ensembl_features2 = na.omit(data.frame(ensembl_features2))
ensembl_features2_plot = ensembl_features2$ensembl_features2

features3 = c("CCND1", "CCNE2", "CCNA2", "CCNB1", "CDK1", "CDK2")
# convert to ensembl IDs
ensembl_features3 = mapIds(org.Hs.eg.db, keys = features3, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
ensembl_features3 = na.omit(data.frame(ensembl_features3))
ensembl_features3_plot = ensembl_features3$ensembl_features3

# Plotting order & colors
ccSeurat_order = c("G1", "S", "G2M")
ccSeurat_colors = c("G1"="#f37f73", "S"="#8571b2", "G2M"="#3db270")
ccAFv2_order = c('qG0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1', 'Unknown')
ccAFv2_colors = c("qG0" = "#d9a428", "G1" = "#f37f73", "Late G1" = "#1fb1a9",  "S" = "#8571b2", "S/G2" = "#db7092", "G2/M" = "#3db270" ,"M/Early G1" = "#6d90ca",  "Unknown" = "#D3D3D3")


#---------------------------------------------------
# Load in data
#--------------------------------------------------

# Set up
tag = 'HSkMSC'
data_dir = 'Qiu/SkMSC/usftp21.novogene.com/01.RawData/HSkMSC_cr/seurat_objects'
resdir2 = 'Qiu/SkMSC/usftp21.novogene.com/01.RawData/HSkMSC_cr/analysis_output'

# Load data
seurat1 = readRDS(file.path(data_dir, 'HSkMSC_cr_normalized_ensembl.rds'))

# Redo scTransform to include all genes for ccAFv2 application and downstream analyses
seurat1 = SCTransform(seurat1, return.only.var.genes = FALSE, verbose = FALSE)
seurat1 = PredictCellCycle(seurat1, do_sctransform=FALSE) # 854 ccAFv2 mgenes

# Order ccAFv2 calls
sub1 = ccAFv2_order %in% factor(seurat1$ccAFv2)
seurat1$ccAFv2 = factor(seurat1$ccAFv2, levels = ccAFv2_order[sub1])

# Load in ccSeurat calls
ccSeurat_calls = read.csv(file.path(resdir2, paste0(tag, '_cr_ccSeurat_calls.csv')), row.names = 'X')
seurat1$Phase = ccSeurat_calls$x

# Order ccSeurat calls
sub2 = ccSeurat_order %in% factor(seurat1$Phase)
seurat1$Phase = factor(seurat1$Phase, levels = ccSeurat_order[sub2])

# Run downstream analysis
seurat1 = RunPCA(seurat1)
#pdf(file.path(resdir2, paste0(tag, '_elbow_plot.pdf')))
#ElbowPlot(seurat1)
#dev.off()
seurat1 = FindNeighbors(seurat1, dims = 1:10)
seurat1 = FindClusters(seurat1)
seurat1 = RunUMAP(seurat1, dims = 1:10)

# Plot
d1 = DimPlot(seurat1, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
d2 = DimPlot(seurat1, reduction = "umap", group.by = "Phase", cols = ccSeurat_colors[sub2])
d3 = DimPlot(seurat1, reduction = "umap", group.by = "ccAFv2", cols = ccAFv2_colors[sub1])
lst = list(d1, d2, d3)
pdf(file.path(resdir2, paste0(tag, '_umap.pdf')), height = 8, width = 8)
grid.arrange(grobs = lst, layout_matrix = rbind(c(1, 2), c(3, NA)), top = '')
# gene sets
p1 = lapply(ensembl_features1_plot, function(goi) {
  if(goi %in% rownames(seurat1)){
    fp1 = FeaturePlot(object = seurat1, features = goi, coord.fixed = TRUE, label=F, pt.size = 0.25) + ggtitle(rownames(ensembl_features1)[ensembl_features1$ensembl_features1 == goi]) + theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.text=element_text(size=8))
    return(fp1 + FontSize(x.title = 10, y.title = 10))
  }
})
grid.arrange(grobs = p1, layout_matrix = rbind(c(1, 2, 3), c(4,5,6), c(7,8,9), c(10, 11, 12)), top = '')
VlnPlot(seurat1, features = ensembl_features3_plot, cols = ccAFv2_colors[sub1], ncol = 2, group.by = 'ccAFv2')
dev.off()

# Save as h5ad or rds file
#seurat1[["RNA"]] = as(object = seurat1[["RNA"]], Class = "Assay")
#SaveH5Seurat(seurat1, file.path(data_dir, paste0(tag, '.h5Seurat')), overwrite = TRUE)
#Convert(file.path(data_dir, paste0(tag, '.h5Seurat')), dest = "h5ad", overwrite = TRUE)
saveRDS(seurat1, file.path(data_dir, paste0(tag, '.rds')))


#---------------------------------------------------
# Remove outlier clusters
#--------------------------------------------------

# Remove cluster 11 and redo downstream analysis
seurat_subset = subset(seurat1, subset = seurat_clusters != 11)
seurat_subset = SCTransform(seurat_subset, return.only.var.genes = FALSE, verbose = FALSE)
seurat_subset = PredictCellCycle(seurat_subset, do_sctransform=FALSE) # 854 ccAFv2 mgenes
seurat_subset = RunPCA(seurat_subset)
#pdf(file.path(resdir2, paste0(tag, '_subset_elbow_plot.pdf')))
#ElbowPlot(seurat_subset)
#dev.off()
seurat_subset = FindNeighbors(seurat_subset, dims = 1:10)
seurat_subset = FindClusters(seurat_subset, resolution = 0.6)
seurat_subset = RunUMAP(seurat_subset, dims = 1:10)

# Plot
d1 = DimPlot(seurat_subset, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
d2 = DimPlot(seurat_subset, reduction = "umap", group.by = "Phase", cols = ccSeurat_colors[sub2])
d3 = DimPlot(seurat_subset, reduction = "umap", group.by = "ccAFv2", cols = ccAFv2_colors[sub1])
lst = list(d1, d2, d3)
pdf(file.path(resdir2, paste0(tag, '_umap_without_outlier_cluster_102124.pdf')), height = 8, width = 8)
grid.arrange(grobs = lst, layout_matrix = rbind(c(1, 2), c(3, NA)), top = '')
# gene sets
p1 = lapply(ensembl_features1_plot, function(goi) {
  if(goi %in% rownames(seurat_subset)){
    fp1 = FeaturePlot(object = seurat_subset, features = goi, coord.fixed = TRUE, label=F, pt.size = 0.25) + ggtitle(rownames(ensembl_features1)[ensembl_features1$ensembl_features1 == goi]) + theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.text=element_text(size=8))
    return(fp1 + FontSize(x.title = 10, y.title = 10))
  }
})
grid.arrange(grobs = p1, layout_matrix = rbind(c(1, 2, 3), c(4,5,6), c(7,8,9), c(10, 11, 12)), top = '')
VlnPlot(seurat_subset, features = ensembl_features3_plot, cols = ccAFv2_colors[sub1], ncol = 2, group.by = 'ccAFv2')
dev.off()

# Functions
plotHeatmap = function(data, tag, group.by = 'seurat_clusters', save_dir = file.path(resdir2)){
  datas = data
  savedir = save_dir
  meta = group.by
  df = table(datas@meta.data[[meta]], datas@meta.data[['ccAFv2']])
  enrich = data.frame(matrix(NA, nrow = nrow(df), ncol = ncol(df)))
  rownames(enrich) = rownames(df)
  colnames(enrich) = colnames(df)
  N = sum(df)
  for (class1 in colnames(df)){
    k = sum(df[,class1])
    for (clust1 in rownames(df)){
      m = sum(df[clust1,])
      q = df[clust1, class1]
      n = N-m
      enrich[clust1,class1] = phyper(q, m, n, k, lower.tail = F)
    }
  }
  pm = -log10(as.matrix(enrich))
  pm[pm>20] = 20
  pdf(file.path(savedir, paste0(tag, '_ccAFv2_', meta, '_heatmap.pdf')))
  print(pheatmap(pm, cluster_cols = F, cluster_rows = F, colorRampPalette(c("white", "red"))(100), display_numbers = round(-log10(as.matrix(enrich)),2)))
  dev.off()
}
plotHeatmap(seurat_subset, tag = tag)

d1 = DimPlot(seurat_subset, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + theme(legend.position="none")
d2 = DimPlot(seurat_subset, reduction = "umap", group.by = "Phase", cols = ccSeurat_colors[sub2]) + theme(legend.position="none")
d3 = DimPlot(seurat_subset, reduction = "umap", group.by = "ccAFv2", cols = ccAFv2_colors[sub1]) + theme(legend.position="none")
lst = list(d1, d2, d3)
pdf(file.path(resdir2, paste0(tag, '_umaps.pdf')), height = 4, width = 12)
grid.arrange(grobs = lst, layout_matrix = rbind(c(1, 2, 3)), top = '')
dev.off()

pdf(file.path(resdir2, paste0(tag, '_ccAFv2_likelihood_feature_umaps_0.9.pdf')), height = 4, width = 10)
FeaturePlot(seurat_subset, c('qG0', 'S'), min.cutoff = 0.9)
dev.off()


# Save as h5ad or rds file
#seurat_subset[["RNA"]] = as(object = seurat_subset[["RNA"]], Class = "Assay")
#SaveH5Seurat(seurat_subset, file.path(data_dir, paste0(tag, '_subset.h5Seurat')), overwrite = TRUE)
#Convert(file.path(data_dir, paste0(tag, '_subset.h5Seurat')), dest = "h5ad", overwrite = TRUE)
#saveRDS(seurat_subset, file.path(data_dir, paste0(tag, '_subset.rds')))
# save out after add S phase clusters to metadata


data1 = readRDS(file.path(data_dir, paste0(tag, '_subset.rds')))
# Find marker genes
cluster_markers = FindAllMarkers(data1, only.pos = TRUE, logfc.threshold = 0.25)
cluster_markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
cluster_markers_genes = mapIds(org.Hs.eg.db, keys = cluster_markers$gene, keytype = 'ENSEMBL', column='SYMBOL', multiVals='first')
cluster_markers$gene = cluster_markers_genes
write.csv(cluster_markers, file.path(resdir2, paste0(tag, '_subset_SCTransform_markers_together.csv')))
top10_symbol = mapIds(org.Hs.eg.db, keys = top10$gene, keytype = 'ENSEMBL', column='SYMBOL', multiVals='first')


#---------------------------------------------------
# Discover DPs
#--------------------------------------------------

# Subset S phase cells for DP analyses
new_fold = 'subset_s'
dir.create(file.path(resdir2, new_fold), showWarnings = FALSE)
savedir = file.path(resdir2, new_fold)

# Subset to S phase cells
s_subset = subset(seurat_subset, subset = ccAFv2 == 'S') # 259 cells
DefaultAssay(s_subset) = 'RNA'
# Redo downstream analysis
s_subset = SCTransform(s_subset, return.only.var.genes = FALSE, verbose = FALSE)
s_subset = RunPCA(s_subset, features = VariableFeatures(s_subset))
#pdf(file.path(savedir, paste0(tag, '_s_elbow_plot_new.pdf')))
#ElbowPlot(s_subset)
#dev.off()
s_subset = FindNeighbors(s_subset, dims = 1:10)
s_subset = FindClusters(s_subset, resolution = 0.3)
s_subset = RunUMAP(s_subset, dims = 1:10)

# Test plot umap
d1 = DimPlot(s_subset, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
d2 = DimPlot(s_subset, reduction = "umap", group.by = "Phase", cols = ccSeurat_colors[sub2])
d3 = DimPlot(s_subset, reduction = "umap", group.by = "ccAFv2", cols = ccAFv2_colors[sub1])
lst = list(d1, d2, d3)
pdf(file.path(savedir, paste0(tag, '_umap_s.pdf')), width = 8, height = 8)
grid.arrange(grobs = lst, layout_matrix = rbind(c(1, 2), c(3, NA)), top = '')
dev.off()

# Save as an h5ad or rds file
#s_subset[["RNA"]] = as(object = s_subset[["RNA"]], Class = "Assay")
#SaveH5Seurat(s_subset, file.path(savedir, paste0(tag, '_s_subset_new.h5Seurat')), overwrite = TRUE)
#Convert(file.path(savedir, paste0(tag, '_s_subset_new.h5Seurat')), dest = "h5ad", overwrite = TRUE)
saveRDS(s_subset, file.path(savedir, paste0(tag, '_s.rds')))

# Find marker genes
cluster_markers = FindAllMarkers(s_subset, only.pos = TRUE, logfc.threshold = 0.25)
cluster_markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
cluster_markers_genes = mapIds(org.Hs.eg.db, keys = cluster_markers$gene, keytype = 'ENSEMBL', column='SYMBOL', multiVals='first')
cluster_markers$gene = cluster_markers_genes
write.csv(cluster_markers, file.path(savedir, paste0(tag, '_SCTransform_markers_together_s.csv')))
top10_symbol = mapIds(org.Hs.eg.db, keys = top10$gene, keytype = 'ENSEMBL', column='SYMBOL', multiVals='first')

# Subset to significant marker genes avg_log2FC > 0.5 and p_val_adj <= 0.05
sig_cluster_markers = cluster_markers %>% dplyr::filter(avg_log2FC>0.5) %>% dplyr::filter(p_val_adj <= 0.05)
write.csv(sig_cluster_markers, file.path(savedir, paste0(tag, '_signficant_SCTransform_markers_together_s.csv')))

# Plot
d1 = DimPlot(s_subset, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
d2 = DimPlot(s_subset, reduction = "umap", group.by = "Phase", cols = ccSeurat_colors[sub2])
d3 = DimPlot(s_subset, reduction = "umap", group.by = "ccAFv2", cols = ccAFv2_colors[sub1])
lst = list(d1, d2, d3)
pdf(file.path(savedir, paste0(tag,'_umap_s.pdf')), height = 10, width = 8)
grid.arrange(grobs = lst, layout_matrix = rbind(c(1, 2), c(3, NA)), top = '')
RidgePlot(s_subset, features = ensembl_features3_plot, ncol=2)
DoHeatmap(object = s_subset, features = names(top10_symbol), size = 4) + scale_y_discrete(labels = top10_symbol)
dev.off()

d5 = DimPlot(s_subset, reduction = "umap", group.by = "seurat_clusters") + theme(legend.position="none")
pdf(file.path(savedir, paste0(tag,'_umap_s_clusters.pdf')), height = 5, width = 5)
d5
dev.off()


# Colorize all cells by new S phase clusters
seurat_subset@meta.data$overlay_S = NA
seurat_subset@meta.data[rownames(s_subset@meta.data),'overlay_S'] = s_subset@meta.data$SCT_snn_res.0.3

# Prepare for plotting
# Order ccAFv2 calls
sub1 = ccAFv2_order %in% factor(seurat_subset$ccAFv2)
seurat_subset$ccAFv2 = factor(seurat_subset$ccAFv2, levels = ccAFv2_order[sub1])

# Order ccSeurat calls
sub2 = ccSeurat_order %in% factor(seurat_subset$Phase)
seurat_subset$Phase = factor(seurat_subset$Phase, levels = ccSeurat_order[sub2])

d1 = DimPlot(seurat_subset, reduction = "umap", group.by = "seurat_clusters")
d2 = DimPlot(seurat_subset, reduction = "umap", group.by = "Phase", cols = ccSeurat_colors[sub2])
d3 = DimPlot(seurat_subset, reduction = "umap", group.by = "ccAFv2", cols = ccAFv2_colors[sub1])
d4 = DimPlot(seurat_subset, reduction = "umap", group.by = "overlay_S")
pdf(file.path(resdir2, paste0(tag,'_umap_s_overlay.pdf')), height = 8, width = 8)
lst = list(d1, d2, d3, d4)
grid.arrange(grobs = lst, layout_matrix = rbind(c(1, 2), c(3, 4)), top = '')
dev.off()

# Save out all data with S phase metadata
saveRDS(seurat_subset, file.path(data_dir, paste0(tag, '_subset.rds')))
