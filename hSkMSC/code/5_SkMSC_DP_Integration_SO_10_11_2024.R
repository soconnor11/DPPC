# Run docker
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
library(presto)
library(keras)

# Set working directory
setwd("/files")

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

tag = 'HSkMSC'
data_dir = 'Qiu/SkMSC/usftp21.novogene.com/01.RawData/HSkMSC_cr/seurat_objects'
resdir2 = 'Qiu/SkMSC/usftp21.novogene.com/01.RawData/HSkMSC_cr/analysis_output'

# Read in data
seurat_subset = readRDS(file.path(data_dir, paste0(tag, '_with_dp_predictions.rds')))

### Integrate DPs
seurat2 = seurat_subset
seurat2[['RNA']] = as(object=seurat2[['RNA']], Class='Assay5')
seurat2[['RNA']] = split(seurat2[["RNA"]], f = seurat2$DP_Predictions)
seurat2 = SCTransform(seurat2, verbose = FALSE)
#seurat2 = PredictCellCycle(seurat2)
seurat2 = RunPCA(seurat2, verbose=FALSE)
seurat2 = FindNeighbors(seurat2, dims = 1:20, reduction = "pca")
seurat2 = FindClusters(seurat2, resolution = 0.8, cluster.name = "unintegrated_clusters")
seurat2 = RunUMAP(seurat2, dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated")
# visualize by batch and cell type annotation

d1 = DimPlot(seurat2, reduction = 'umap.unintegrated', label=F, group.by = 'unintegrated_clusters') + ggtitle('seurat_clusters')
d2 = DimPlot(seurat2, reduction = 'umap.unintegrated', label=F, group.by = 'Phase', cols = ccSeurat_colors[sub2]) + ggtitle('Phase')
d3 = DimPlot(seurat2, reduction = 'umap.unintegrated', label=F, group.by = 'ccAFv2', cols = ccAFv2_colors[sub1]) + ggtitle('ccAFv2')
d4 = DimPlot(seurat2, reduction = 'umap.unintegrated', label=F, group.by = 'DP_Predictions')

lst = list(d1, d2, d3, d4)
pdf(file.path(resdir2, paste0(tag, '_umap_dp_unintegrated_new.pdf')), height = 5, width = 20)
grid.arrange(grobs = lst, layout_matrix = rbind(c(1, 2, 3, 4)), top = '')
dev.off()

seurat_unintegrated = seurat2

# CCA
seurat2 = IntegrateLayers(object = seurat_unintegrated, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", normalization.method = "SCT", verbose = FALSE)
seurat2 = FindNeighbors(seurat2, reduction = "integrated.cca", dims = 1:20)
seurat2 = FindClusters(seurat2, resolution = 0.8, cluster.name = "cca_clusters")
seurat2 = RunUMAP(seurat2, reduction = "integrated.cca", dims = 1:20, reduction.name = "umap.cca")
# Plot
d1 = DimPlot(seurat2, reduction = 'umap.cca', label=F, group.by = 'cca_clusters') + ggtitle('integrated_clusters')
d2 = DimPlot(seurat2, reduction = 'umap.cca', label=F, group.by = 'Phase', cols = ccSeurat_colors[sub2]) + ggtitle('Phase')
d3 = DimPlot(seurat2, reduction = 'umap.cca', label=F, group.by = 'ccAFv2', cols = ccAFv2_colors[sub1]) + ggtitle('ccAFv2')
d4 = DimPlot(seurat2, reduction = 'umap.cca', label=F, group.by = 'DP_Predictions')
lst = list(d1, d2, d3, d4)
pdf(file.path(resdir2, paste0(tag, '_umap_dp_integrated_cca_new.pdf')), height = 5, width = 20)
grid.arrange(grobs = lst, layout_matrix = rbind(c(1, 2, 3, 4)), top = '')
dev.off()













# Old code


# Load in DPPCs
#dppc_predictions = read.csv('old_home/home/soconnor/ccNN/ccAFv3/results/HSkMSC/HSkMSC_DPPC_classifier_predictions_051714_0.5.csv', row.names ='X')
#dppc_predictions = read.csv('old_home/home/soconnor/ccNN/ccAFv3/results/HSkMSC_DPPC_classifier_predictions_093024_0.5.csv', row.names ='X')
dppc_predictions = read.csv('old_home/home/soconnor/ccNN/ccAFv3/results/HSkMSC_DPPC_classifier_predictions_093024_0.5_new.csv', row.names ='X')
colnames(dppc_predictions) = 'x'
#dppc_prob = read.csv('old_home/home/soconnor/ccNN/ccAFv3/results/SkMSC_DPPC_classifier_predictions_093024_0.5_probabilities.csv', row.names ='X')
dppc_prob = read.csv('old_home/home/soconnor/ccNN/ccAFv3/results/SkMSC_DPPC_classifier_predictions_093024_0.5_probabilities_new.csv', row.names ='X')

seurat_subset$dppc_predictions = dppc_predictions$x
seurat_subset$DPPC1 = dppc_prob$DPPC1
seurat_subset$DPPC2 = dppc_prob$DPPC2
#seurat_subset$DPPC3 = dppc_prob$DPPC3

# Colorize by DPPC
seurat2 = seurat_subset
seurat2 = RunPCA(seurat2, verbose=FALSE)
seurat2 = FindNeighbors(seurat2, dims = 1:12, verbose=FALSE)
seurat2 = FindClusters(seurat2, verbose=FALSE)
seurat2 = RunUMAP(seurat2, dims=1:12, verbose=FALSE)

d1 = DimPlot(seurat2, reduction = 'umap', label=F, group.by = 'seurat_clusters') + ggtitle('seurat_clusters')
d3 = DimPlot(seurat2, reduction = 'umap', label=F, group.by = 'ccAFv2', cols = ccAFv2_colors[sub1]) + ggtitle('ccAFv2')
d4 = DimPlot(seurat2, reduction = 'umap', label=F, group.by = 'dppc_predictions')

lst = list(d1, d3, d4)
pdf(file.path(resdir2, 'dppc_predictions_093024_new.pdf'), width = 16, height = 4)
grid.arrange(grobs = lst, layout_matrix = rbind(c(1, 2, 3)), top = '')
dev.off()

#### Regress out DPPCs
seurat2 = seurat_subset
#seurat2 = SCTransform(seurat2, return.only.var.genes = FALSE, vars.to.regress = c('DPPC1', 'DPPC2', 'DPPC3'))
seurat2 = SCTransform(seurat2, return.only.var.genes = FALSE, vars.to.regress = c('DPPC1', 'DPPC2'))
seurat2 = PredictCellCycle(seurat2, do_sctransform=FALSE)
seurat2 = RunPCA(seurat2, verbose=FALSE)

pdf(file.path(resdir2, 'regressed_elbow_plot.pdf'))
ElbowPlot(seurat2)
dev.off()

seurat2 = FindNeighbors(seurat2, dims = 1:20, verbose=FALSE)
seurat2 = FindClusters(seurat2, verbose=FALSE, resolution = 0.6)
seurat2 = RunUMAP(seurat2, dims=1:20, verbose=FALSE)

d1 = DimPlot(seurat2, reduction = 'umap', label=F, group.by = 'seurat_clusters') + ggtitle('seurat_clusters')
d3 = DimPlot(seurat2, reduction = 'umap', label=F, group.by = 'ccAFv2', cols = ccAFv2_colors[sub1]) + ggtitle('ccAFv2')
d4 = DimPlot(seurat2, reduction = 'umap', label=F, group.by = 'dppc_predictions')
lst = list(d1, d3, d4)

pdf(file.path(resdir2, 'dppc_predictions_regressed_093024.pdf'), width = 10, height = 10)
grid.arrange(grobs = lst, layout_matrix = rbind(c(1, 2), c(3, NA)), top = '')
# gene sets
p1 = lapply(ensembl_features1_plot, function(goi) {
  if(goi %in% rownames(seurat2)){
    fp1 = FeaturePlot(object = seurat2, features = goi, coord.fixed = TRUE, label=F, pt.size = 0.25) + ggtitle(rownames(ensembl_features1)[ensembl_features1$ensembl_features1 == goi]) + theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.text=element_text(size=8))
    return(fp1 + FontSize(x.title = 10, y.title = 10))
  }
})
grid.arrange(grobs = p1, layout_matrix = rbind(c(1, 2, 3), c(4,5,6), c(7,8,9)), top = '')
p2 = lapply(myogenic_to_plot, function(goi) {
  if(goi %in% rownames(seurat2)){
    fp1 = FeaturePlot(object = seurat2, features = goi, coord.fixed = TRUE, label=F, pt.size = 0.25) + ggtitle(rownames(myogenic2)[myogenic2$myogenic2 == goi]) + theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.text=element_text(size=8))
    return(fp1 + FontSize(x.title = 10, y.title = 10))
  }
})
grid.arrange(grobs = p2, layout_matrix = rbind(c(1, 2, 3), c(4,5,6), c(NA,NA,NA)), top = '')
VlnPlot(seurat2, features = ensembl_features3_plot, cols = ccAFv2_colors[sub1], ncol = 2, group.by = 'ccAFv2')
dev.off()

plotHeatmap(seurat2, 'hSkMSC')

saveRDS(seurat2, file.path(data_dir, paste0(tag, '_dppc_regressed.rds')))

# Find marker genes
seurat2 = readRDS(file.path(data_dir, paste0(tag, '_dppc_regressed.rds')))
# Find marker genes
cluster_markers= FindAllMarkers(seurat2, only.pos = TRUE, logfc.threshold = 0.25)
cluster_markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
cluster_markers_genes = mapIds(org.Hs.eg.db, keys = cluster_markers$gene, keytype = 'ENSEMBL', column='SYMBOL', multiVals='first')
cluster_markers$gene = cluster_markers_genes
write.csv(cluster_markers, file.path(resdir2, paste0(tag, '_SCTransform_markers_together_dppc_regressed.csv')))
top10_symbol = mapIds(org.Hs.eg.db, keys = top10$gene, keytype = 'ENSEMBL', column='SYMBOL', multiVals='first')


### Integrate DPPCs

seurat2 = seurat_subset
seurat2[["RNA"]] <- split(seurat2[["RNA"]], f = seurat2$dppc_predictions)
seurat2 = SCTransform(seurat2, verbose = FALSE)
#seurat2 = PredictCellCycle(seurat2)
seurat2 = RunPCA(seurat2, verbose=FALSE)
seurat2 <- FindNeighbors(seurat2, dims = 1:20, reduction = "pca")
seurat2 <- FindClusters(seurat2, resolution = 0.8, cluster.name = "unintegrated_clusters")
seurat2 <- RunUMAP(seurat2, dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated")
# visualize by batch and cell type annotation

d1 = DimPlot(seurat2, reduction = 'umap.unintegrated', label=F, group.by = 'unintegrated_clusters') + ggtitle('seurat_clusters')
d3 = DimPlot(seurat2, reduction = 'umap.unintegrated', label=F, group.by = 'ccAFv2', cols = ccAFv2_colors[sub1]) + ggtitle('ccAFv2')
d4 = DimPlot(seurat2, reduction = 'umap.unintegrated', label=F, group.by = 'dppc_predictions')

lst = list(d1, d3, d4)
pdf(file.path(resdir2, 'umap_unintegrated.pdf'), width = 8, height = 8)
grid.arrange(grobs = lst, layout_matrix = rbind(c(1, 2), c(3, NA)), top = '')
dev.off()

seurat_unintegrated = seurat2

# CCA
seurat2 <- IntegrateLayers(object = seurat_unintegrated, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", normalization.method = "SCT", verbose = FALSE)
seurat2 <- FindNeighbors(seurat2, reduction = "integrated.cca", dims = 1:20)
seurat2 <- FindClusters(seurat2, resolution = 0.8, cluster.name = "cca_clusters")
seurat2 <- RunUMAP(seurat2, reduction = "integrated.cca", dims = 1:20, reduction.name = "umap.cca")
# Plot
d1 = DimPlot(seurat2, reduction = 'umap.cca', label=F, group.by = 'unintegrated_clusters') + ggtitle('unintegrated_clusters')
d2 = DimPlot(seurat2, reduction = 'umap.cca', label=F, group.by = 'cca_clusters') + ggtitle('integrated_clusters')
d3 = DimPlot(seurat2, reduction = 'umap.cca', label=F, group.by = 'ccAFv2', cols = ccAFv2_colors[sub1]) + ggtitle('ccAFv2')
d4 = DimPlot(seurat2, reduction = 'umap.cca', label=F, group.by = 'dppc_predictions')
lst = list(d1, d2, d3, d4)
pdf(file.path(resdir2, 'umap_integrated_cca.pdf'), width = 8, height = 8)
grid.arrange(grobs = lst, layout_matrix = rbind(c(1, 2), c(3, 4)), top = '')
dev.off()

# RPCA
seurat2 <- IntegrateLayers(object = seurat_unintegrated, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca", normalization.method = "SCT", verbose = FALSE)
seurat2 <- FindNeighbors(seurat2, reduction = "integrated.rpca", dims = 1:20)
seurat2 <- FindClusters(seurat2, resolution = 0.8, cluster.name = "rpca_clusters")
seurat2 <- RunUMAP(seurat2, reduction = "integrated.rpca", dims = 1:20, reduction.name = "umap.rpca")
# Plot
d1 = DimPlot(seurat2, reduction = 'umap.rpca', label=F, group.by = 'unintegrated_clusters') + ggtitle('unintegrated_seurat_clusters')
d2 = DimPlot(seurat2, reduction = 'umap.rpca', label=F, group.by = 'rpca_clusters') + ggtitle('integrated_clusters')
d3 = DimPlot(seurat2, reduction = 'umap.rpca', label=F, group.by = 'ccAFv2', cols = ccAFv2_colors[sub1]) + ggtitle('ccAFv2')
d4 = DimPlot(seurat2, reduction = 'umap.rpca', label=F, group.by = 'dppc_predictions')
lst = list(d1, d2, d3, d4)
pdf(file.path(resdir2, 'umap_integrated_rpca.pdf'), width = 8, height = 8)
grid.arrange(grobs = lst, layout_matrix = rbind(c(1, 2), c(3, 4)), top = '')
dev.off()


# Harmony
seurat2 <- IntegrateLayers(
  object = seurat2, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
seurat2 <- FindNeighbors(seurat2, reduction = "harmony", dims = 1:20)
seurat2 <- FindClusters(seurat2, resolution = 0.8, cluster.name = "harmony_clusters")
seurat2 <- RunUMAP(seurat2, reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony")
# Plot
d1 = DimPlot(seurat2, reduction = 'umap.harmony', label=F, group.by = 'unintegrated_clusters') + ggtitle('unintegrated_clusters')
d2 = DimPlot(seurat2, reduction = 'umap.harmony', label=F, group.by = 'harmony_clusters') + ggtitle('integrated_clusters')
d3 = DimPlot(seurat2, reduction = 'umap.harmony', label=F, group.by = 'ccAFv2', cols = ccAFv2_colors[sub1]) + ggtitle('ccAFv2')
d4 = DimPlot(seurat2, reduction = 'umap.harmony', label=F, group.by = 'dppc_predictions')
lst = list(d1, d2, d3, d4)
pdf(file.path(resdir2, 'umap_integrated_harmony.pdf'), width = 8, height = 8)
grid.arrange(grobs = lst, layout_matrix = rbind(c(1, 2), c(3, 4)), top = '')
dev.off()


# MNN
devtools::install_github("satijalab/seurat-wrappers", ref = "seurat5", force = TRUE, upgrade = "never")
SeuratWrappers::scVIIntegration()
SeuratWrappers::FastMNNIntegration()

seurat2 <- IntegrateLayers(
  object = seurat_unintegrated, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)




#----- Recluster based off DPPCs -----#

# Set up
tag = 'HSkMSC'
data_dir = 'Qiu/SkMSC/usftp21.novogene.com/01.RawData/HSkMSC_cr/seurat_objects'
resdir2 = 'Qiu/SkMSC/usftp21.novogene.com/01.RawData/HSkMSC_cr/analysis_output'

# Load data
seurat1 = readRDS(file.path(data_dir, 'HSkMSC_cr_normalized_ensembl.rds'))

# Apply ccAFv2
seurat1 = PredictCellCycle(seurat1) # 854 ccAFv2 mgenes

# Order ccAFv2 calls
sub1 = ccAFv2_order %in% factor(seurat1$ccAFv2)
seurat1$ccAFv2 <- factor(seurat1$ccAFv2, levels = ccAFv2_order[sub1])

phase_calls = read.csv(file.path(resdir2, paste0('HSkMSC_cr_ccSeurat_calls.csv')), row.names = 'X')
seurat1$Phase = phase_calls$x

sub2 = ccSeurat_order %in% factor(seurat1$Phase)
seurat1$Phase <- factor(seurat1$Phase, levels = ccSeurat_order[sub2])

dppc_predictions = read.csv('old_home/home/soconnor/ccNN/ccAFv3/results/HSkMSC/HSkMSC_DPPC_classifier_predictions_051714_0.5.csv', row.names ='X')
colnames(dppc_predictions) = 'x'

seurat1$dppc_predictions = dppc_predictions$x

# Add all module scores
seurat1 = AddModuleScore(seurat1, features = list(osteogenic_to_plot), name = 'osteogenic')
seurat1 = AddModuleScore(seurat1, features = list(chondrogenic_to_plot), name = 'chronogenic')
seurat1 = AddModuleScore(seurat1, features = list(adipogenic_to_plot), name = 'adipogenic')
seurat1 = AddModuleScore(seurat1, features = list(myogenic_to_plot), name = 'myogenic')
seurat1 = AddModuleScore(seurat1, features = list(neurogenic_to_plot), name = 'neurogenic')
seurat1 = AddModuleScore(seurat1, features = list(mgenes_neuralg0), name = 'neuralg0')
seurat1 = AddModuleScore(seurat1, features = list(mgenes_g1), name = 'g1')
seurat1 = AddModuleScore(seurat1, features = list(mgenes_mg1), name = 'mg1')
seurat1 = AddModuleScore(seurat1, features = list(mes_progen_ensembl), name = 'mes progenitor')



#seurat2 = seurat1
#seurat2 = RunPCA(seurat2, dims = 1:30, verbose=FALSE)
#seurat2 = FindNeighbors(seurat2, dims = 1:30, verbose=FALSE)
#seurat2 = FindClusters(seurat2, verbose=FALSE)
#seurat2 = RunUMAP(seurat2, dims=1:30, verbose=FALSE)


# Split into DPPCS
dppc1_subset = subset(seurat1, subset = dppc_predictions == 'DPPC1')
dppc2_subset = subset(seurat1, subset = dppc_predictions == 'DPPC2')

# Put them together in a list
datas = list()
datas[['dppc1']] = dppc1_subset
datas[['dppc2']] = dppc2_subset

# Recluster DPPCs
#for (data1 in c('dppc1', 'dppc2')){
for (data1 in c('dppc2')){
  # Subset to specific dppc
  cat('\nReclustering', data1,'\n')
  dppc_subset = datas[[data1]]
  dir.create(file.path(resdir2, data1), showWarnings = F)
  savedir2 = file.path(resdir2, data1)
  #dir.create(file.path(savedir2, 'regress_counts'), showWarnings = F)
  #savedir2 = file.path(savedir2, 'regress_counts')
  # Recluster
  cat('\nPreprocessing...\n')
  DefaultAssay(dppc_subset) = 'RNA'
  dppc_subset = SCTransform(dppc_subset, return.only.var.genes = FALSE, verbose = FALSE)
  #dppc_subset = SCTransform(dppc_subset, return.only.var.genes = FALSE, verbose = FALSE, vars.to.regress = c('nCount_RNA'))
  dppc_subset = RunPCA(dppc_subset, features = VariableFeatures(dppc_subset), verbose = FALSE)
  dppc_subset = FindNeighbors(dppc_subset, dims = 1:10, verbose = FALSE)
  cat('\nTesting different resolutions for FindClusters function\n')
  for (res1 in c(0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)){
  #for (res1 in c(1.1)){
    cat('\n', as.character(res1), '\n')
    dppc_subset = FindClusters(dppc_subset, resolution = res1)
    dppc_subset = RunUMAP(dppc_subset, dims = 1:10)
    #saveRDS(dppc_subset, file.path(savedir2, paste0(tag, '_', data1, '_res_', as.character(res1), '.rds')))
    #write.csv(dppc_subset@reductions$umap@cell.embeddings, file.path(savedir2, paste0(tag, "_", as.character(res1), "_UMAP_Coordinates_for_scvelo.csv")))
    # Save out data for scvelo
    #data_loom <- as.loom(dppc_subset, file.path(savedir2, paste0(tag, '_', as.character(res1), '_data_052124.loom')), verbose = FALSE, overwrite=TRUE)
    #data_loom
    #data_loom$close_all()
    # Plot heatmap of hypergeometrics (cells in a ccAFv2 state cells vs. cells in a cluster)
    plotHeatmap(dppc_subset, tag = res1, savedir2) # saved out as a pdf
    cluster_markers= FindAllMarkers(dppc_subset, only.pos = TRUE, logfc.threshold = 0.25)
    # Plot heatmap of hypergeometrics (ccAFv2 state marker genes vs. marker genes from a cluster)
    #plotHeatmap_mgenes(dppc_subset, tag = res1, markers = cluster_markers, log2fc_cutoff = 0.8, save_dir = savedir2) # saved out as a pdf
    cluster_markers %>%
        group_by(cluster) %>%
        dplyr::filter(avg_log2FC > 0.5) %>%
        slice_head(n = 10) %>%
        ungroup() -> top10
    cluster_markers_genes = mapIds(org.Hs.eg.db, keys = cluster_markers$gene, keytype = 'ENSEMBL', column='SYMBOL', multiVals='first')
    cluster_markers$gene = cluster_markers_genes
    write.csv(cluster_markers, file.path(savedir2, paste0(tag, '_', data1, '_res_', as.character(res1), '_SCTransform_markers_together.csv')))
    top10_symbol = mapIds(org.Hs.eg.db, keys = top10$gene, keytype = 'ENSEMBL', column='SYMBOL', multiVals='first')
    d1 = DimPlot(dppc_subset, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
    d2 = DimPlot(dppc_subset, reduction = "umap", group.by = "Phase", cols = ccSeurat_colors[sub2])
    d3 = DimPlot(dppc_subset, reduction = "umap", group.by = "ccAFv2", cols = ccAFv2_colors[sub1])
    d4 = DimPlot(dppc_subset, reduction = "umap", group.by = "dppc_predictions")
    #d5 = DimPlot(dppc_subset, reduction = "umap", group.by = "bm_predictions")
    lst = list(d1, d2, d3, d4)
    # Plot
    cat('\nPlotting\n')
    pdf(file.path(savedir2, paste0(tag,'_', data1, '_res_', as.character(res1), '_test.pdf')), height = 10, width = 10)
    grid.arrange(grobs = lst, layout_matrix = rbind(c(1, 2), c(3, 4)), top = '')
    print(VlnPlot(dppc_subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 2))
    print(FeaturePlot(dppc_subset, features = c('osteogenic1', 'chronogenic1', 'adipogenic1', 'myogenic1', 'neurogenic1', 'mes progenitor1')))
    print(VlnPlot(object = dppc_subset, features = c('osteogenic1', 'chronogenic1', 'adipogenic1', 'myogenic1', 'neurogenic1', 'mes progenitor1')))
    print(FeaturePlot(dppc_subset, features = c('neuralg01', 'g11', 'mg11', sk_quiescence_to_plot), reduction = 'umap', ncol=2))
    print(VlnPlot(dppc_subset, features = c('neuralg01', 'g11', 'mg11', sk_quiescence_to_plot), ncol=2))
    #--- ccAFv2 vs. cluster ids stacked barplot ---#
    cf <- table(dppc_subset$ccAFv2, dppc_subset$seurat_clusters)
    totals <- colSums(cf)
    cnewdf <- rbind(cf, totals)
    cf_1 = matrix(ncol=length(unique(dppc_subset$seurat_clusters)), nrow=length(unique(dppc_subset$ccAFv2)))
    for(i in c(1:length(unique(dppc_subset$seurat_clusters)))){
      for(n in c(1:length(unique(dppc_subset$ccAFv2)))) {
        cf_1[n,i] = cnewdf[n,i]/cnewdf[length(unique(dppc_subset$ccAFv2))+1, i]
      }
    }
    colnames(cf_1) = colnames(cf)
    rownames(cf_1) = rownames(cf)
    sub4 = rownames(data.frame(ccAFv2_colors)) %in% rownames(cf_1)
    par(mar = c(8, 8, 8, 8) + 2.0)
    barplot(cf_1, xlab = '', ylab = 'Cell Percentage', las=2, legend.text = rownames(cf_1),  col = ccAFv2_colors[sub1], args.legend=list(x=ncol(cf_1) + 4.5, y=max(colSums(cf_1)), bty = 'n'))
    print(RidgePlot(dppc_subset, features = ensembl_features3_plot, ncol=2))
    print(DoHeatmap(object = dppc_subset, features = names(top10_symbol), size = 4) + scale_y_discrete(labels = top10_symbol))
    # ccAF Heatmap
    #print(DoHeatmap(object = dppc_subset, features=goi_lst_ensembl2) + scale_y_discrete(labels = names(goi_lst_ensembl2)))
    dev.off()
  }
}
