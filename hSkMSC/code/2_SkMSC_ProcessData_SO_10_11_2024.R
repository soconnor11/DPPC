# docker run -it -v '/media:/files' cplaisier/ccafv2_seurat5_presto

#--------------------------------
# Set up section / load packages
#--------------------------------

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
library(scales)
library(pheatmap)
#devtools::install_github('immunogenomics/presto')
library(presto)

# Set working directory
setwd("files/")

mgenes = read.csv(system.file('extdata', 'ccAFv2_genes.csv', package='ccAFv2'), header=TRUE, row.names=1)[,paste0('human_ensembl')]

# Some features to investigate
features1 = c('S100B', 'SOX2', 'SOX4', 'MKI67', 'PAX7','PAX3', 'MYOD1', 'ACT1', 'CD82', 'EMD')
# convert to ensembl IDs
ensembl_features1 = mapIds(org.Hs.eg.db, keys = features1, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
ensembl_features1 = na.omit(data.frame(ensembl_features1))
ensembl_features1_plot = ensembl_features1$ensembl_features1

features2 <- c("MBP", "PLP1", "ETNPPL", "CD14","CX3CR1","PTPRC", "RBFOX3")
# convert to ensembl IDs
ensembl_features2 = mapIds(org.Hs.eg.db, keys = features2, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
ensembl_features2 = na.omit(data.frame(ensembl_features2))
ensembl_features2_plot = ensembl_features2$ensembl_features2

features3 <- c("CCND1", "CCNE2", "CCNA2", "CCNB1", "CDK1", "CDK2")
# convert to ensembl IDs
ensembl_features3 = mapIds(org.Hs.eg.db, keys = features3, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
ensembl_features3 = na.omit(data.frame(ensembl_features3))
ensembl_features3_plot = ensembl_features3$ensembl_features3

# Plotting order & colors
ccSeurat_order = c("G1", "S", "G2M")
ccSeurat_colors = c("G1" = "#f37f73", "S" = "#8571b2", "G2M" = "#3db270")
ccAFv2_order = c('qG0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1', 'Unknown')
ccAFv2_colors = c("qG0" = "#d9a428", "G1" = "#f37f73", "Late G1" = "#1fb1a9",  "S" = "#8571b2", "S/G2" = "#db7092", "G2/M" = "#3db270" ,"M/Early G1" = "#6d90ca",  "Unknown" = "#D3D3D3")

osteogenic = c('RUNX2','SP7', 'ALPL', 'BGLAP', 'COL1A1', 'OGN', 'ALPI', 'ASPN')
osteogenic2 = mapIds(org.Hs.eg.db, keys = osteogenic, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
osteogenic2 = na.omit(data.frame(osteogenic2))
osteogenic_to_plot = osteogenic2$osteogenic2

chondrogenic = c('SOX5', 'SOX6', 'COL10A1', 'COL11A1', 'ACAN', 'COMP', 'SOX9', 'EPYC', 'HAPLN1', 'COL2A1')
chondrogenic2 = mapIds(org.Hs.eg.db, keys = chondrogenic, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
chondrogenic2 = na.omit(data.frame(chondrogenic2))
chondrogenic_to_plot = chondrogenic2$chondrogenic2

adipogenic = c('PPARG', 'CIDEC', 'ADIPOQ', 'CD36', 'FABP4', 'AOC3', 'PGC1A')
adipogenic2 = mapIds(org.Hs.eg.db, keys = adipogenic, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
adipogenic2 = na.omit(data.frame(adipogenic2))
adipogenic_to_plot = adipogenic2$adipogenic2

myogenic = c('MYF4', 'MYOD1', 'MYF5', 'ACTA2', 'TNNT3', 'MYLPF', 'MCAM1', 'SOX10', 'MYOG')
myogenic2 = mapIds(org.Hs.eg.db, keys = myogenic, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
myogenic2 = na.omit(data.frame(myogenic2))
myogenic_to_plot = myogenic2$myogenic2

sk_quiescence = c('PAX7')
sk_quiescence2 = mapIds(org.Hs.eg.db, keys = sk_quiescence, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
sk_quiescence2 = na.omit(data.frame(sk_quiescence2))
sk_quiescence_to_plot = sk_quiescence2$sk_quiescence2

neurogenic = c('NES', 'GFAP', 'PAX6', 'SOX2', 'NCAM1', 'VIM', 'NEUROD1')
neurogenic2 = mapIds(org.Hs.eg.db, keys = neurogenic, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
neurogenic2 = na.omit(data.frame(neurogenic2))
neurogenic_to_plot = neurogenic2$neurogenic2

paracrine_signals = c('TGFBA', 'VEGFA', 'HGF', 'IGF1', 'CXCL12', 'IL6', 'FGF2', 'PDGFA')
paracrine2 = mapIds(org.Hs.eg.db, keys = paracrine_signals, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
paracrine2 = na.omit(data.frame(paracrine2))
paracrine_to_plot = paracrine2$paracrine2

quiescence_markers = c('NT5E', 'ENG', 'NGFR', 'HIC1', 'CYP1B1', 'MMP13')
quiescence2 = mapIds(org.Hs.eg.db, keys = quiescence_markers, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
quiescence2 = na.omit(data.frame(quiescence2))
quiescence_to_plot = quiescence2$quiescence2

stemness_markers = c('SOX4', 'GAS1', 'DPP4')
stemness2 = mapIds(org.Hs.eg.db, keys = stemness_markers, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
stemness2 = na.omit(data.frame(stemness2))
stemness_to_plot = stemness2$stemness2

early_bmp = c('PRDM1', 'HOXD1', 'HIVEP2', 'CEBPB', 'NR1D2', 'IRX3', 'KLF9')
early_bmp2 = mapIds(org.Hs.eg.db, keys = early_bmp, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
early_bmp2 = na.omit(data.frame(early_bmp2))
early_bmp_to_plot = early_bmp2$early_bmp2

late_bmp = c('LEF1', 'IRF2', 'SOX9', 'TBL1XR1', 'FOXC1', 'FOXP1', 'FOXN3', 'NR3C1')
late_bmp2 = mapIds(org.Hs.eg.db, keys = late_bmp, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
late_bmp2 = na.omit(data.frame(late_bmp2))
late_bmp_to_plot = late_bmp2$late_bmp2

# Woods et al., 2021
mes_progen = c('CD200', 'OSTN', 'POSTN', 'IBSP', 'ADIPQ', 'CLEC2D', 'SERPINE2', 'MMP13', 'TIMP1', 'RSPO2', 'DKK3', 'CXCL12', 'KITLG')
mes_progen_ensembl = mapIds(org.Hs.eg.db, keys = mes_progen, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
mes_progen_ensembl = mes_progen_ensembl[!is.na(mes_progen_ensembl)]

mgenes = read.csv(system.file('extdata', 'ccAFv2_genes.csv', package='ccAFv2'), header=TRUE)
mgenes_neuralg0 = mgenes[mgenes$Neural.G0 == 1,]$human_ensembl
mgenes_g1 = mgenes[mgenes$G1 == 1,]$human_ensembl
mgenes_mg1 = mgenes[mgenes$M.Early.G1 == 1,]$human_ensembl

# Functions
plotHeatmap = function(data, tag, save_dir = file.path(resdir2)){
  datas = data
  savedir = save_dir
  df = table(datas$seurat_clusters, datas$ccAFv2)
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
  pdf(file.path(savedir, paste0(tag, '_ccAFv2_seurat_clusters_heatmap.pdf')))
  print(pheatmap(pm, cluster_cols = F, cluster_rows = F, colorRampPalette(c("white", "red"))(100), display_numbers = round(-log10(as.matrix(enrich)),2)))
  dev.off()
}


# Process data individually
scProcessData = function(res_dir, tag, cutoff = 0.5, assay = 'SCT', do_sctransform = TRUE, species= 'human', resolution = 0.8, save_dir = 'analysis_output', obj_dir = 'seurat_objects', clusters_to_remove = F, symbol = F, norm_regress = F){
  cat('\n',tag,'\n')
  # Set up folders
  resdir1 = file.path(res_dir, tag)
  resdir2 = file.path(resdir1, save_dir)
  resdir3 = file.path(resdir1, obj_dir)
  #---------------------------
  # Load filtered / normalized data
  #---------------------------
  gene_id = 'ensembl'
  seurat2 = readRDS(file.path(resdir3, paste0(tag, '_normalized_', paste0(gene_id),'.rds')))
  cat('\n', dim(seurat2), '\n')
  # Load ccSeurat calls
  ccseurat_calls = read.csv(file.path(resdir2, paste0(tag, '_ccSeurat_calls.csv')), row.names = 'X')
  seurat2 <- AddMetaData(seurat2, ccseurat_calls, col.name = 'Phase')
  # Order ccSeurat calls
  sub1 = ccSeurat_order %in% factor(seurat2$Phase)
  seurat2$Phase <- factor(seurat2$Phase, levels = ccSeurat_order[sub1])
  #---------------------------
  # Classify with ccAFv2
  #---------------------------
  seurat2 = PredictCellCycle(seurat2, cutoff = cutoff, assay = assay, do_sctransform = do_sctransform, species = species, gene_id = gene_id)
  # Order ccAFv2 calls
  sub2 = ccAFv2_order %in% factor(seurat2$ccAFv2)
  seurat2$ccAFv2 <- factor(seurat2$ccAFv2, levels = ccAFv2_order[sub2])
  tmp = data.frame(table(seurat2$ccAFv2))
  rownames(tmp) = tmp$Var1
  write.csv(seurat2$ccAFv2, file.path(resdir2, paste0(tag, '_ccAFv2_calls.csv')))
  write.csv(data.frame((tmp['Freq']/dim(seurat2)[2])*100), file.path(resdir2, paste0(tag, '_ccAFv2_call_frequency.csv')))
  #---------------------------
  # Normalize
  #---------------------------
  #seurat2 = SCTransform(seurat2, verbose = FALSE) # started with normalized object; do not need to redo
  # Add marker gene counts for each cell
  #seurat_subset = seurat2[mgenes]
  tmp = GetAssayData(seurat2, slot= 'counts')
  tmp2 = tmp[rownames(tmp) %in% mgenes,]
  non_zero_mgenes = colSums(tmp2 > 0)
  cat('Cells that have non-zero ccAFv2 genes: ', length(non_zero_mgenes), '\n')
  # Add as meta data column
  seurat2 = AddMetaData(seurat2, non_zero_mgenes, col.name = 'ccAFv2_mgene_counts')
  seurat2 = RunPCA(seurat2, verbose=FALSE)
  seurat2 = FindNeighbors(seurat2, dims = 1:10, verbose=FALSE)
  seurat2 = FindClusters(seurat2, verbose=FALSE, resolution = resolution)
  seurat2 = RunUMAP(seurat2, dims=1:10, verbose=FALSE)
  #Idents(seurat2) = seurat2$ccAFv2
  # Find cluster marker genes for each ccAFv2 class
  cluster_markers = FindAllMarkers(seurat2, only.pos = TRUE)
  cluster_markers %>%
      group_by(cluster) %>%
      dplyr::filter(avg_log2FC > 0.5) %>%
      slice_head(n = 10) %>%
      ungroup() -> top10
  cluster_markers_genes = mapIds(org.Hs.eg.db, keys = cluster_markers$gene, keytype = 'ENSEMBL', column='SYMBOL', multiVals='first')
  cluster_markers$gene = cluster_markers_genes
  write.csv(cluster_markers, file.path(resdir2, paste0(tag,'_scTransform_Markers_together.csv')))
  top10_symbol = mapIds(org.Hs.eg.db, keys = top10$gene, keytype = 'ENSEMBL', column='SYMBOL', multiVals='first')
  #------------------------------------------------------
  # Plotting
  #---------------------------------------------------
  cat('Plotting UMAPs and ccAFv2 barplot \n')
  d1 = DimPlot(seurat2, reduction = 'umap', label=F, group.by = 'seurat_clusters') + ggtitle('seurat_clusters')
  d2 = DimPlot(seurat2, reduction = 'umap', label=F, group.by = 'Phase', cols = ccSeurat_colors)  + ggtitle('Phase')
  d3 = DimPlot(seurat2, reduction = 'umap', label=F, group.by = 'ccAFv2', cols = ccAFv2_colors[sub2]) + ggtitle('ccAFv2')
  pdf(file.path(resdir2, paste0(tag, '.pdf')), width = 10, height = 8)
  lst = list(d1, d2, d3)
  grid.arrange(grobs = lst, layout_matrix = rbind(c(1, 2), c(3, NA)), top = '')
  # gene sets
  p1 = lapply(ensembl_features1_plot, function(goi) {
    if(goi %in% rownames(seurat2)){
      fp1 = FeaturePlot(object = seurat2, features = goi, coord.fixed = TRUE, label=F, pt.size = 0.25) + ggtitle(rownames(ensembl_features1)[ensembl_features1$ensembl_features1 == goi]) + theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.text=element_text(size=8))
      return(fp1 + FontSize(x.title = 10, y.title = 10))
    }
  })
  grid.arrange(grobs = p1, layout_matrix = rbind(c(1, 2, 3, 4), c(5,6,7,8), c(9,10,11, 12)), top = '')
  #------- ccAFv2 vs. seurat clusters - cell percentages -------#
  par(mar = c(8, 8, 8, 8) + 2.0)
  barplot(table(seurat2$ccAFv2, seurat2$seurat_clusters), beside = FALSE, col = ccAFv2_colors[sub2], xlab = 'clusters ID', ylab = 'Cell count', legend.text = rownames(table(seurat2$ccAFv2, seurat2$seurat_clusters)), args.legend=list(title='ccAFv2 classification'))
  #--- ccAFv2 vs. cluster ids stacked barplot ---#
  cf <- table(seurat2$ccAFv2, seurat2$seurat_clusters)
  totals <- colSums(cf)
  data.frame(totals)
  cnewdf <- rbind(cf, totals)
  cf_1 = matrix(ncol=length(unique(seurat2$seurat_clusters)), nrow=length(unique(seurat2$ccAFv2)))
  for(i in c(1:length(unique(seurat2$seurat_clusters)))){
    for(n in c(1:length(unique(seurat2$ccAFv2)))) {
      cf_1[n,i] = cnewdf[n,i]/cnewdf[length(unique(seurat2$ccAFv2))+1, i]
    }
  }
  colnames(cf_1) = colnames(cf)
  rownames(cf_1) = rownames(cf)
  sub4 = rownames(data.frame(ccAFv2_colors)) %in% rownames(cf_1)
  par(mar = c(8, 8, 8, 8) + 2.0)
  barplot(cf_1, xlab = '', ylab = 'Cell Percentage', las=2, legend.text = rownames(cf_1),  col = ccAFv2_colors[sub4], args.legend=list(x=ncol(cf_1) + 4.5, y=max(colSums(cf_1)), bty = 'n'))
  #--- ccAFv2 vs. ccSeurat stacked barplot ---#
  cf <- table(seurat2$ccAFv2, seurat2$Phase)
  totals <- colSums(cf)
  data.frame(totals)
  cnewdf <- rbind(cf, totals)
  cf_1 = matrix(ncol=length(unique(seurat2$Phase)), nrow=length(unique(seurat2$ccAFv2)))
  for(i in c(1:length(unique(seurat2$Phase)))){
    for(n in c(1:length(unique(seurat2$ccAFv2)))) {
      cf_1[n,i] = cnewdf[n,i]/cnewdf[length(unique(seurat2$ccAFv2))+1, i]
    }
  }
  colnames(cf_1) = colnames(cf)
  rownames(cf_1) = rownames(cf)
  sub4 = rownames(data.frame(ccAFv2_colors)) %in% rownames(cf_1)
  par(mar = c(8, 8, 8, 8) + 2.0)
  barplot(cf_1, xlab = '', ylab = 'Cell Percentage', las=2, legend.text = rownames(cf_1),  col = ccAFv2_colors[sub4], args.legend=list(x=ncol(cf_1) + 1.5, y=max(colSums(cf_1)), bty = 'n'))
  #--- ccAFv2 and number of marker genes boxplot ---#
  v1 = VlnPlot(seurat2, features = 'ccAFv2_mgene_counts', group.by = 'ccAFv2', cols = ccAFv2_colors[sub4]) + theme(legend.position = 'none') + xlab('ccAFv2') + ylab('ccAFv2 marker gene counts')
  v2 = VlnPlot(seurat2, features = 'nCount_RNA', group.by = 'ccAFv2', cols = ccAFv2_colors[sub4]) + theme(legend.position = 'none') + xlab('ccAFv2') + ylab('nCount_RNA')
  v3 = VlnPlot(seurat2, features = 'nFeature_RNA', group.by = 'ccAFv2', cols = ccAFv2_colors[sub4]) + theme(legend.position = 'none') + xlab('ccAFv2') + ylab('nFeature_RNA')
  lst2 = list(v2, v3, v1)
  grid.arrange(grobs = lst2, layout_matrix = rbind(c(1, 2), c(3, NA)), top = '')
  print(DoHeatmap(object = seurat2, features = names(top10_symbol), group.colors = ccAFv2_colors[sub4], size = 4) + scale_y_discrete(labels = top10_symbol))
  print(RidgePlot(seurat2, features = ensembl_features3_plot, ncol=2, cols = ccAFv2_colors[sub4]))
  dev.off()
  # Change factored metadata to characters
  #seurat2$ccAFv2 = as.character(seurat2$ccAFv2)
  #seurat2$Phase = as.character(seurat2$Phase)
  #cat('saving processed data as loom and rds...\n')
  #data_loom_2 <- as.loom(seurat2, file.path(resdir3, paste0(tag, '_processed.loom')), verbose = FALSE, overwrite = TRUE)
  #data_loom_2$close_all()
  #saveRDS(seurat2, file.path(resdir3, paste0(tag, '_processed.rds')))
  return(seurat2)
}

#---------------------------------------------------
# Downstream analysis & plotting (all together)
#--------------------------------------------------

scProcessData(res_dir = 'Qiu/SkMSC/usftp21.novogene.com/01.RawData', tag = 'HSkMSC_cr', resolution = 0.8)
