
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
# Read in data
tag = 'HSkMSC'
data_dir = 'Qiu/SkMSC/usftp21.novogene.com/01.RawData/HSkMSC_cr/seurat_objects'
resdir2 = 'Qiu/SkMSC/usftp21.novogene.com/01.RawData/HSkMSC_cr/analysis_output/subset_s'
savedir = 'Qiu/SkMSC/usftp21.novogene.com/01.RawData/HSkMSC_cr/analysis_output'


# Load data
seurat_subset = readRDS(file.path(data_dir, paste0(tag, '_subset.rds')))
s_subset = readRDS(file.path(resdir2, paste0(tag, '_s.rds')))
sig_cluster_markers = read.csv(file.path(resdir2, paste0(tag, '_signficant_SCTransform_markers_together_s.csv')), row.names = 'X')

# List of marker genes to train on
sig_cluster_marker_genes = intersect(unique(rownames(sig_cluster_markers)),rownames(s_subset))

## 1. Set up data
# TODO:  Generalize the cluster information for y_train
x_train = scale(t(as.matrix(s_subset@assays$SCT$scale.data)[sig_cluster_marker_genes,]))
DPs = s_subset@meta.data$SCT_snn_res.0.3
numDPs = length(unique(DPs))
y_train = to_categorical(DPs, numDPs)


## 2. Bulid classifier based on the S phase cells
model = keras_model_sequential()
model %>%
    layer_dense(units=length(sig_cluster_marker_genes), input_shape=length(sig_cluster_marker_genes), activation='relu') %>% # Input layer
    layer_dense(units=200, activation='relu') %>% # First hidden layer
    layer_dropout(0.5) %>% # First dropout regularization layer
    layer_dense(units=100, activation='relu') %>% # Second hidden layer
    layer_dropout(0.5) %>% # Second dropout regularization layer
    layer_dense(numDPs, activation='softmax')

model %>% compile(
  optimizer = 'sgd',
  loss = 'categorical_crossentropy',
  metrics = list('accuracy')
)

history = model %>% fit(
  x_train, y_train,
  epochs = 10, batch_size = 5,
  validation_split = 0.2
)

pdf(file.path(resdir2, paste0(tag,'_keras_training.pdf')), height = 10, width = 8)
plot(history)
dev.off()


## 3. Apply to all the other cells in the dataset
x_all = scale(t(as.matrix(seurat_subset@assays$SCT$scale.data)[sig_cluster_marker_genes,]))

# Get likelihoods
likelihoods = model %>% predict(x_all)
for(i in 1:ncol(likelihoods)) {
    seurat_subset@meta.data[[paste0('DP',i)]] = likelihoods[,i]
}

# Get classes
seurat_subset@meta.data[['DP_Predictions']] = paste0('DP',1+as.numeric(likelihoods %>% k_argmax()))


# Plot
sub1 = ccAFv2_order %in% factor(seurat_subset$ccAFv2)
sub2 = ccSeurat_order %in% factor(seurat_subset$Phase)
d1 = DimPlot(seurat_subset, reduction = "umap", group.by = "orig.ident", cols = "#D3D3D3") + theme(legend.position="none")
#d2 = DimPlot(seurat_subset, reduction = "umap", group.by = "seurat_clusters")
#d3 = DimPlot(seurat_subset, reduction = "umap", group.by = "Phase", cols = ccSeurat_colors[sub2])
d4 = DimPlot(seurat_subset, reduction = "umap", group.by = "ccAFv2", cols = ccAFv2_colors[sub1]) + theme(legend.position="none")
d5 = DimPlot(seurat_subset, reduction = "umap", group.by = "overlay_S") + theme(legend.position="none")
d6 = DimPlot(seurat_subset, reduction = "umap", group.by = "DP_Predictions") + theme(legend.position="none")
lst = list(d1, d4, d5, d6)
pdf(file.path(savedir, paste0(tag,'_umap_colorized_all.pdf')), height = 5, width = 20)
grid.arrange(grobs = lst, layout_matrix = rbind(c(1, 2, 3, 4)), top = '')
dev.off()

# Save out RDS object
saveRDS(seurat_subset, file.path(data_dir, paste0(tag, '_with_dp_predictions.rds')))
