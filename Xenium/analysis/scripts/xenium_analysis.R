library(Seurat)
library(sctransform)
# library(future)
# plan("multisession", workers = 4)
# plan('sequential')
# library(ggplot2)
library(tidyverse)
library(scuttle)
library(scater)
library(ggpubr)
library(ggthemes)
library(viridis)
library(Polychrome)
library(MetBrewer)
library(ggdark)
library(glmGamPoi)
library(scales)
library(SpatialExperiment)
library(SpaceTrooper)
library(Banksy)
library(SeuratWrappers)
library(gridExtra)
library(pals)

options(future.seed = TRUE)
set.seed(281330800)

# 1. ** SET WORKING DIRECTORY WHERE PROJECT OUPUTS WILL SAVE TO ** ----
# results_folder = 'z:/aguilada/xenium/results'
# figs_folder = 'z:/aguilada/xenium/results'
# data_folder = 'z:/aguilada/xenium/data'

results_folder = "c:/Users/david/Documents/Research/PhD/Spatial_Omics/Xenium/results/"
figs_folder = "c:/Users/david/Documents/Research/PhD/Spatial_Omics/Xenium/figures/"
data_folder = "c:/Users/david/Documents/Research/PhD/Spatial_Omics/Xenium/data/"
annotations_folder = "c:/Users/david/Documents/Research/PhD/Spatial_Omics/Xenium/annotations/"

dirs = c(results_folder,figs_folder,data_folder , annotations_folder)

for (dir in dirs) {
  if (!dir.exists(dir)) { dir.create(dir,
                                     recursive=TRUE) }
}

# ** SET PATH TO FOLDER CONTAINING XENIUM DATA ** ----
xenium_folder = "d:/aguilada/xenium/10190-JC/20240322__211909__10190-JC/20240322__211909__10190-JC/output-XETG00077__0014270__10190-JC-4_ROI_B__20240322__211926/"

# general files (some are supplemental files)
settings_path = paste0(xenium_folder, 'experiment.xenium')
he_img_path = "d:/aguilada/xenium/10190-JC/Post-Xenium_H&E/0014270_24-401_20x.ome.tif"
# if_img_path = paste0(xenium_folder, 'Xenium_FFPE_Human_Breast_Cancer_Rep1_if_image.tif')
# panel_meta_path = paste0(xenium_folder, 'Xenium_FFPE_Human_Breast_Cancer_Rep1_panel.tsv') # (optional)

# files (SUBCELLULAR): (tutorial focuses on working with these files)
cell_bound_path = paste0(xenium_folder, 'cell_boundaries.csv.gz')
nuc_bound_path = paste0(xenium_folder, 'nucleus_boundaries.csv.gz')
tx_path = paste0(xenium_folder, 'transcripts.csv.gz')
feat_meta_path = paste0(xenium_folder, 'cell_feature_matrix/features.tsv.gz') # (also used in aggregate)

# files (AGGREGATE):
expr_mat_path = paste0(xenium_folder, 'cell_feature_matrix')
cell_meta_path = paste0(xenium_folder, 'cells.csv.gz') # contains spatlocs

# load data ----
# xenium.obj <- LoadXenium(xenium_folder,
#                          fov = "fov",
#                          flip.xy=TRUE,
#                          mols.qv.threshold = 20)# Remove transcript molecules with a QV less than this threshold. QV >= 20 is the standard
# saveRDS(xenium.obj, paste0(data_folder,'xenium.rds'))
# xenium.obj <- LoadXenium(xenium_folder,
#                          fov = "fov",
#                          segmentations = 'nucleus',
#                          flip.xy=TRUE,
#                          mols.qv.threshold = 20) # Remove transcript molecules with a QV less than this threshold. QV >= 20 is the standard
# saveRDS(xenium.obj, paste0(data_folder,'xenium_nucleus_segmented.rds'))
# xenium.obj <- LoadXenium(xenium_folder,
#                          fov = "fov",
#                          segmentations = 'cell',
#                          flip.xy=TRUE,
#                          mols.qv.threshold = 20) # R
# saveRDS(xenium.obj, paste0(data_folder,'xenium_cell_expansion_segmented.rds'))

xenium.obj = readRDS(paste0(data_folder, 'xenium_nucleus_segmented.rds'))

# Add metrics to metadata ----
## add cells ids and coordinates to metadata ----
xenium.obj@meta.data$cells = colnames(xenium.obj)
xenium.obj = AddMetaData(xenium.obj, 
                         metadata=xenium.obj@images[["fov"]]@boundaries[["centroids"]]@coords)

summary(xenium.obj@meta.data$x)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 377    2976    4067    3897    4879    6160


summary(xenium.obj@meta.data$y)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 156.7  2288.3  3244.4  3429.0  4780.2  6428.4 

## calculate log10 z-scores and complexity ----
calc_metrics = function(x) {
  
  log10genes = log10(x$nFeature_Xenium)
  log10counts = log10(x$nCount_Xenium)
  
  # mean and standard deviation of log10 ngenes and ncount
  mean_ngene = mean(log10genes)
  mean_ncount = mean(log10counts)
  std_ngene = sd(log10genes)
  std_ncount = sd(log10counts)
  
  # z-scores of log10 ngenes and ncounts
  x$log10_genes_zscored = (log10genes - mean_ngene) / std_ngene 
  x$log10_counts_zscored = (log10counts - mean_ncount) / std_ncount
  
  # calculate log10 genes per log10 UMI for each cell 
  # acts as a novelty score to measure complexity
  x$complexity = log10genes/log10counts
  
  return(x)
}
# add metrics to metadata
xenium.obj = calc_metrics(xenium.obj)


# proposed QC cutoffs ----
# gene_cutoff_low = 150
# gene_cutoff_high = 5000
# umi_cutoff_low = 150
# umi_cutoff_high = 20000 
complexity_cutoff = 0.8 # log10genes / log10UMIs
min_cells_per_gene = 3 # filtered genes must be expressed in at least this # of cells. Keep low, at 3, for rare cell types


# qc plots ----
VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), 
        ncol = 2, pt.size = 0, assay = 'Xenium', group.by = 'orig.ident')
ggsave(paste0(figs_folder,'Vln_nCounts_and_nFeatures.png'), dpi=320, width=20, height=20)

png(paste0(figs_folder,'histogram_complexity.png'), res=320, units='in', width=10, height=10)
hist(xenium.obj$complexity)
dev.off()

summary(xenium.obj$nFeature_Xenium)
summary(xenium.obj$nCount_Xenium)
summary(xenium.obj$complexity)

# spatial plots for datasets acquired through
plotCoords(xenium.obj, point_size=0, y_reverse = FALSE) + ggtitle('Xenium')


# filter data ----
## cell level filtering  ----
# xenium.obj = readRDS(paste0(data_folder, 'xenium.rds'))
ncells = length(xenium.obj$orig.ident) # for tracking # cells before filtering

xenium.obj = subset(xenium.obj, subset=(nCount_Xenium > 0) &
                      (nFeature_Xenium > 0))

### number of cell before and after cell-level filtering ----
ncells_post_cell_filter = length(xenium.obj$orig.ident)
cat(paste0(ncells,' cells before cell-level filtering\n',
           ncells_post_cell_filter,' cells left after cell-level filtering\n',
           ncells - ncells_post_cell_filter,' cells removed after cell-level filtering'))

# 37086 cells before cell-level filtering
# 37084 cells left after cell-level filtering
# 2 cells removed after cell-level filtering

## gene level filtering ----
counts = LayerData(xenium.obj, assay='Xenium', layer='counts')

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero = counts > 0

# Sums all TRUE values and returns TRUE if more than 20 TRUE values per gene
keep_genes = Matrix::rowSums(nonzero) >= min_cells_per_gene # gene expressed in at least 3 cells

# Only keeping those genes expressed in more than 20 cells
filtered_counts = counts[keep_genes, ]

# Reassign to filtered Seurat object and filtered metadata
xenium.obj@assays[["Xenium"]]@layers[["counts"]] = filtered_counts

## number of genes before and after filtering
ngenes = sum(Matrix::rowSums(nonzero) > 0)
ngenes_post_filtering = sum(Matrix::rowSums(nonzero) >min_cells_per_gene)
cat(paste0(ngenes,' genes before gene-level filtering.\n',
           ngenes_post_filtering,' genes left after gene-level filtering.\n',
           ngenes - ngenes_post_filtering,' genes removed after cell-level filtering.'))
# 380 genes before gene-level filtering.
# 380 genes left after gene-level filtering.
# 0 genes removed after cell-level filtering.

# normalize data ----
options(future.globals.maxSize = 2 * 1024^3)  # 2 GB
xenium.obj <- SCTransform(xenium.obj, assay = "Xenium",)
# saveRDS(xenium.obj, paste0(data_folder, 'xenium_obj_sct.rds'))
# xenium.obj = readRDS(paste0(data_folder, 'xenium_obj_sct.rds'))

xenium.obj <- NormalizeData(xenium.obj, assay = "Xenium",
                            normalization.method='LogNormalize')

# clustering ----
## clustering with seurat vignette ---
xenium.obj <- RunPCA(xenium.obj, npcs = 50, 
                     features = rownames(xenium.obj),
                     assay='SCT',
                     layer='data')


ElbowPlot(object = xenium.obj, reduction='pca',
          ndims = 50)+
  ggtitle('Elbow plot of PCs')+
  theme_classic2()+
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))
ggsave(paste0(figs_folder,'elbow_plot_seurat.png'), dpi=320, width=20, height= 10)

### Determine percent of variation associated with each PC
pct = xenium.obj[["pca"]]@stdev / sum(xenium.obj[["pca"]]@stdev) * 100
pct
#
# # Calculate cumulative percents for each PC
cumu = cumsum(pct)
cumu
#
# # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 = which(cumu >= 90 & pct <= 5)[1]
co1 # 43
#
# # Determine the difference between variation of PC and subsequent PC.
# # This is the last point where change of % of variation is more than 0.1%. Afterwards, % of variation change is < 0.1%
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2 # 18
#
# # Minimum of the two calculations
pcs = min(co1, co2)
pcs 
# [1] 19


resolutions = c(0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4)

xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:50) # use all pcs since used SCTransform
xenium.obj <- FindClusters(xenium.obj, resolution = resolutions)
xenium.obj <- RunUMAP(xenium.obj, 
                      reduction = 'pca', 
                      dims = 1:50)
# xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:50)
# xenium.obj <- FindClusters(xenium.obj, resolution = resolutions)

# SCT_snn_res.0.3 clusters: 11
# SCT_snn_res.0.4 clusters: 12
# SCT_snn_res.0.6 clusters: 18
# SCT_snn_res.0.7 clusters: 20
# SCT_snn_res.0.8 clusters: 22
# SCT_snn_res.0.9 clusters: 22
# SCT_snn_res.1 clusters: 22
# SCT_snn_res.1.2 clusters: 26
# SCT_snn_res.1.4 clusters: 29


### visualize clustering results with clustree ----
library(clustree)

clustree(xenium.obj, prefix='SCT_snn_res.')
ggsave(paste0(figs_folder,'clustree.png'), dpi=320, width=20, height= 10)

clustree(seurat_umap_sct, prefix = 'SCT_snn_res.',
         node_colour = '', node_colour_aggr = 'mean')

### UMAPs ----


Idents(xenium.obj) = 'SCT_snn_res.0.4'
# Idents(xenium.obj) = 'SCT_snn_res.0.8'
# Idents(xenium.obj) = 'SCT_snn_res.1.4'
# Idents(xenium.obj) = 'SCT_snn_res.1.2'
DimPlot(xenium.obj,
        reduction = "umap",
        label=T,
        label.size=6)

DefaultAssay(xenium.obj) = 'Xenium'
# xenium.obj = NormalizeData(xenium.obj)

xenium.obj = NormalizeData(xenium.obj, 
                           assay='Xenium',
                           normalization.method = 'RC',
                           scale.factor=100)

# DefaultAssay(xenium.obj) = 'SCT'
FeaturePlot(xenium.obj, features = c('EPCAM'),
            reduction='umap',
            min.cutoff='q50',
            max.cutoff='q90',
            cols=viridis(256)) &
  DarkTheme()

DefaultAssay(xenium.obj) = 'Xenium'
DotPlot(xenium.obj, features = c('PECAM1'),
        assay='Xenium')






### spatial plots ----

ImageDimPlot(xenium.obj, cols = "polychrome", size = 1.5, 
             group.by="SCT_snn_res.0.8", dark.background = FALSE)

ImageDimPlot(xenium.obj, cols = "polychrome", size = 1.5, 
             group.by="SCT_snn_res.1.4", dark.background = TRUE)

FeatureScatter(xenium.obj, 'x', 'y', pt.size = 0.75, cols= hue_pal()(22),
               group.by='SCT_snn_res.1,4')+
  guides(color = guide_legend(override.aes = list(size = 3)))

cropped.coords <- Crop(xenium.obj[["fov"]], x = c(3500, 5250), y = c(2500, 3750), coords = "plot")
xenium.obj[["zoom"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentations"
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = FALSE, molecules = c('EPCAM','ADIPOQ','PECAM1','CD8A'), nmols = 10000)

DefaultAssay(xenium.obj) = 'SCT'
ImageFeaturePlot(xenium.obj, features = c('EPCAM','ADIPOQ','PECAM1','CD8A'), 
                 max.cutoff = rep('q50',4), 
                 size = 0.75, fov = 'zoom',
                 cols = c('white','red'))





SpatialDimPlot(xenium.obj, group.by = "SCT_snn_res.0.8", 
               label = T, 
               repel = T, 
               label.size = 4)


# DefaultAssay(xenium.obj) = 'SCT'

FeatureScatter(xenium.obj, 'x', 'y', pt.size = 0.75, cols= hue_pal()(22),
               group.by='SCT_snn_res.0.8')+
  guides(color = guide_legend(override.aes = list(size = 3)))

FeatureScatter(xenium.obj, 'x', 'y', pt.size = 0.75, cols= hue_pal()(22),
               group.by='SCT_snn_res.0.8')+
  guides(color = guide_legend(override.aes = list(size = 3)))


## function to plot initial umaps ----
plot_umap = function (seurat_umap, resolution, reduction='umap', directory) {
  
  DefaultAssay(seurat_umap) = 'SCT'
  
  # Assign umap resolution
  Idents(seurat_umap) = resolution
  
  # create directory if it doesn't exist
  if (!dir.exists(paste0(directory,resolution))) { dir.create(paste0(directory,resolution),
                                                                          recursive=TRUE) }
  
  # Plot the UMAP
  DimPlot(seurat_umap,
          reduction = reduction,
          label=TRUE,
          repel=TRUE,
          label.size=6)
  ggsave(paste0(directory,resolution,'/UMAP.png'), dpi=320, width=10, height=10)
}

## UMAP clustering QC ----
metrics =  c("nCount_Xenium", "nFeature_Xenium", "complexity") # metrics to overlay on umap

### create function to plot all QC plots and umap overlays ----
umap_clustering_qc = function(seurat_umap, resolution, reduction='umap', directory) {
  
  DefaultAssay(seurat_umap) = 'SCT'
  
  # create directory if it doesn't exist
  if (!dir.exists(paste0(directory,resolution))) { dir.create(paste0(directory,resolution),
                                                              recursive=TRUE) }
  
  # set umap resolution
  Idents(seurat_umap) = resolution
  seurat_umap@meta.data$resolution = Idents(seurat_umap)
  
  # Metrics to plot on UMAP
  FeaturePlot(seurat_umap, 
              reduction = reduction, 
              features = metrics,
              pt.size = 0.4, 
              order = TRUE,
              min.cutoff = 'q10',
              label = TRUE)
  ggsave(paste0(directory,resolution,'/counts_features_metrics.png'), dpi=320, width=20, height=20)
  
  # Boxplot of nUMIs per cluster
  ggplot(seurat_umap@meta.data) +
    geom_boxplot(aes(x=resolution, y=nCount_Xenium, fill=resolution))+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    scale_y_log10(n.breaks=10, limits=c(100,NA),labels = comma_format())+
    annotation_logticks(sides = "l")+
    NoLegend()+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5))
  ggsave(paste0(directory,resolution,'/nUMIs_per_cluster.png'), dpi=320, width=20, height= 10)
  
  # Boxplot of nGenes per cluster. 
  ggplot(seurat_umap@meta.data) +
    geom_boxplot(aes(x=resolution, y=nFeature_Xenium, fill=resolution))+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    scale_y_log10(n.breaks=10, limits=c(100,NA),labels = comma_format())+
    annotation_logticks(sides = "l")+
    NoLegend()+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5))
  ggsave(paste0(directory,resolution,'/nGenes_per_cluster.png'), dpi=320, width=20, height= 10)
  
}
DefaultAssay(xenium.obj) = 'SCT'
plot_umap(xenium.obj,resolution ='SCT_snn_res.0.8', directory=figs_folder)
plot_umap(xenium.obj,resolution ='SCT_snn_res.1.4', directory=figs_folder)

umap_clustering_qc(xenium.obj,resolution ='SCT_snn_res.0.8', directory=figs_folder)
umap_clustering_qc(xenium.obj,resolution ='SCT_snn_res.1.4', directory=figs_folder)



### create csv for xenium explorer

df = data.frame(
  cell_id = names(xenium.obj$SCT_snn_res.1.4),
  group = xenium.obj$SCT_snn_res.1.4
)

write.csv(df, paste0(annotations_folder, 'custom_cell_groups_SCT_snn_res.1.4.csv'))


## clustering with BANKSY using SCTransform ----

# Kelly palette for visualization
mypal <- kelly()[-1]

DefaultAssay(xenium.obj) = 'Xenium'
xenium.obj = NormalizeData(xenium.obj,
                           assay = 'Xenium',
                           normalization.method = 'RC',
                           scale.factor = 100)
xenium.obj = FindVariableFeatures(xenium.obj, nfeatures = 50)
# Identify the 10 most highly variable genes
top20 <- head(VariableFeatures(xenium.obj), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(xenium.obj)

plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot2


# xenium.obj <- RunBanksy(xenium.obj, lambda = 0.2, verbose=TRUE,
#                         assay = 'Xenium', slot = 'data', features = 'variable',
# k_geom = 15)
xenium.obj <- RunBanksy(xenium.obj, lambda = 0.2, verbose=TRUE,
                        assay = 'SCT', slot = 'data', features = 'variable',
                        k_geom = 15)

xenium.obj <- RunPCA(xenium.obj, 
                     assay='BANKSY',
                     npcs = 50,
                     reduction.name='pca.banksy',
                     features = rownames(xenium.obj))


ElbowPlot(object = xenium.obj, reduction='pca.banksy',
          ndims = 50)+
  ggtitle('Elbow plot of PCs')+
  theme_classic2()+
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))
ggsave(paste0(figs_folder,'elbow_plot_banksy.png'), dpi=320, width=20, height= 10)

### Determine percent of variation associated with each PC
pct = xenium.obj[["pca.banksy"]]@stdev / sum(xenium.obj[["pca.banksy"]]@stdev) * 100
pct
#
# # Calculate cumulative percents for each PC
cumu = cumsum(pct)
cumu
#
# # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 = which(cumu >= 90 & pct <= 5)[1]
co1 # 44
#
# # Determine the difference between variation of PC and subsequent PC.
# # This is the last point where change of % of variation is more than 0.1%. Afterwards, % of variation change is < 0.1%
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2 # 11
#
# # Minimum of the two calculations
pcs = min(co1, co2)
pcs 
# [1] 11


xenium.obj <- FindNeighbors(xenium.obj, 
                            reduction = "pca.banksy", 
                            dims = 1:50)

xenium.obj <- FindClusters(xenium.obj, 
                           cluster.name = 'banksy_cluster',
                           # resolution = 0.5,
                           resolution = 2.0,
                           graph.name='BANKSY_snn')

xenium.obj <- RunUMAP(xenium.obj, 
                      reduction = 'pca.banksy', 
                      reduction.name='umap.banksy',
                      dims = 1:50,
                      )
# SCT_snn_res.0.3 clusters: 14
# SCT_snn_res.0.4 clusters: 17
# SCT_snn_res.0.5 clusters: 
# SCT_snn_res.0.6 clusters: 17    4 singletons identified. 13 final clusters.
# SCT_snn_res.0.7 clusters: 19
# SCT_snn_res.0.8 clusters: 21    4 singletons identified. 17 final clusters.
# SCT_snn_res.0.9 clusters: 21
# SCT_snn_res.1 clusters: 23      4 singletons identified. 19 final clusters.
# SCT_snn_res.1.2 clusters: 26    4 singletons identified. 22 final clusters.
# SCT_snn_res.1.4 clusters: 28    4 singletons identified. 24 final clusters.

# SCT_snn_res.1.8 clusters: 28    4 singletons identified. 26 final clusters.
# SCT_snn_res.2.0 clusters: 33    4 singletons identified. 29 final clusters.

saveRDS(xenium.obj, paste0(data_folder,'xenium_sct_clustered_seurat_banksy.rds'))


## UMAPs Banksy ----



## viz spatial ----
grid.arrange(
  DimPlot(xenium.obj, group.by = 'banksy_cluster', label = TRUE, 
          reduction='umap.banksy', label.size = 6, shuffle=TRUE),
  SpatialDimPlot(xenium.obj, stroke = NA, label = TRUE, label.size = 3, 
                 repel = TRUE, alpha = 0.5, pt.size.factor = 2),
  ncol = 2
)

DefaultAssay(xenium.obj) = 'Xenium'
# DefaultAssay(xenium.obj) = 'SCT'
FeaturePlot(xenium.obj, features = c('CD8A'),
            reduction='umap.banksy',
            min.cutoff='q50',
            max.cutoff='q90',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(xenium.obj, features = c('PECAM1'),
        assay='Xenium')

Idents(xenium.obj) = 'banksy_cluster'
# Idents(xenium.obj) = 'SCT_snn_res.1.4'
# Idents(xenium.obj) = 'SCT_snn_res.1.2'
DimPlot(xenium.obj,
        reduction = "umap",
        group.by = 'banksy_cluster',
        label=T,
        label.size=6)

plot_umap(xenium.obj,resolution ='banksy_cluster', 
          reduction='umap.banksy',directory=figs_folder)
umap_clustering_qc(xenium.obj,resolution ='banksy_cluster', 
                   reduction='umap.banksy',directory=figs_folder)


### banksy spatial plots ----
ImageDimPlot(xenium.obj, cols = "polychrome", size = 1.5,
             group.by = 'banksy_cluster',dark.background = FALSE)

ImageDimPlot(xenium.obj, cols = "polychrome", size = 1.5,
             group.by = 'banksy_cluster',dark.background = TRUE)

cropped.coords <- Crop(xenium.obj[["fov"]], x = c(4250, 4500), y = c(3000, 3200), coords = "plot")
xenium.obj[["zoom"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentations"
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", 
             border.size = 0.1, cols = "polychrome",coord.fixed = FALSE, 
             molecules = c('EPCAM','ADIPOQ','PECAM1','CD8A'), nmols = 10000,
             group.by = 'banksy_cluster', dark.background=TRUE)

DefaultAssay(xenium.obj) = 'SCT'
ImageFeaturePlot(xenium.obj, features = c('EPCAM','ADIPOQ','PECAM1','CD8A'), 
                 max.cutoff = rep('q50',4), 
                 size = 0.75, fov = 'zoom',
                 cols = c('white','red'))


SpatialDimPlot(xenium.obj, group.by = "SCT_snn_res.0.8", 
               label = T, 
               repel = T, 
               label.size = 4)



## clustering with Banksy using NormalizeData RC method ----
xenium.obj = NormalizeData(xenium.obj,
                           assay = 'Xenium',
                           normalization.method = 'RC',
                           scale.factor = 100)
xenium.obj = FindVariableFeatures(xenium.obj, nfeatures = 50)
# Identify the 10 most highly variable genes
top20 <- head(VariableFeatures(xenium.obj), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(xenium.obj)

plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot2


xenium.obj <- RunBanksy(xenium.obj, lambda = 0.2, verbose=TRUE,
                        assay = 'Xenium', slot = 'data', features = 'variable',
                        k_geom = 15)
# xenium.obj <- RunBanksy(xenium.obj, lambda = 0.2, verbose=TRUE, 
#                         assay = 'SCT', slot = 'data', features = 'variable',
#                         k_geom = 15)

xenium.obj <- RunPCA(xenium.obj, 
                     assay='BANKSY',
                     npcs = 50,
                     reduction.name='pca.banksy',
                     features = rownames(xenium.obj))


ElbowPlot(object = xenium.obj, reduction='pca.banksy',
          ndims = 50)+
  ggtitle('Elbow plot of PCs')+
  theme_classic2()+
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))
ggsave(paste0(figs_folder,'elbow_plot_banksy.png'), dpi=320, width=20, height= 10)

### Determine percent of variation associated with each PC
pct = xenium.obj[["pca.banksy"]]@stdev / sum(xenium.obj[["pca.banksy"]]@stdev) * 100
pct
#
# # Calculate cumulative percents for each PC
cumu = cumsum(pct)
cumu
#
# # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 = which(cumu >= 90 & pct <= 5)[1]
co1 # 44
#
# # Determine the difference between variation of PC and subsequent PC.
# # This is the last point where change of % of variation is more than 0.1%. Afterwards, % of variation change is < 0.1%
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2 # 11
#
# # Minimum of the two calculations
pcs = min(co1, co2)
pcs 
# [1] 11


xenium.obj <- FindNeighbors(xenium.obj, 
                            reduction = "pca.banksy", 
                            dims = 1:50)

xenium.obj <- FindClusters(xenium.obj, 
                           cluster.name = 'banksy_cluster',
                           # resolution = 0.5,
                           resolution = 2.0,
                           graph.name='BANKSY_snn')

xenium.obj <- RunUMAP(xenium.obj, 
                      reduction = 'pca.banksy', 
                      reduction.name='umap.banksy',
                      dims = 1:50,
)
# SCT_snn_res.0.3 clusters: 14
# SCT_snn_res.0.4 clusters: 17
# SCT_snn_res.0.5 clusters: 
# SCT_snn_res.0.6 clusters: 17    4 singletons identified. 13 final clusters.
# SCT_snn_res.0.7 clusters: 19
# SCT_snn_res.0.8 clusters: 21    4 singletons identified. 17 final clusters.
# SCT_snn_res.0.9 clusters: 21
# SCT_snn_res.1 clusters: 23      4 singletons identified. 19 final clusters.
# SCT_snn_res.1.2 clusters: 26    4 singletons identified. 22 final clusters.
# SCT_snn_res.1.4 clusters: 28    4 singletons identified. 24 final clusters.

# SCT_snn_res.1.8 clusters: 28    4 singletons identified. 26 final clusters.
# SCT_snn_res.2.0 clusters: 33    4 singletons identified. 29 final clusters.

saveRDS(xenium.obj, paste0(data_folder,'xenium_sct_clustered_seurat_banksy.rds'))


## UMAPs Banksy ----



## viz spatial ----
grid.arrange(
  DimPlot(xenium.obj, group.by = 'banksy_cluster', label = TRUE, 
          reduction='umap.banksy', label.size = 6, shuffle=TRUE),
  SpatialDimPlot(xenium.obj, stroke = NA, label = TRUE, label.size = 3, 
                 repel = TRUE, alpha = 0.5, pt.size.factor = 2),
  ncol = 2
)

DefaultAssay(xenium.obj) = 'Xenium'
# DefaultAssay(xenium.obj) = 'SCT'
FeaturePlot(xenium.obj, features = c('CD8A'),
            reduction='umap.banksy',
            min.cutoff='q50',
            max.cutoff='q90',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(xenium.obj, features = c('PECAM1'),
        assay='Xenium')

Idents(xenium.obj) = 'banksy_cluster'
# Idents(xenium.obj) = 'SCT_snn_res.1.4'
# Idents(xenium.obj) = 'SCT_snn_res.1.2'
DimPlot(xenium.obj,
        reduction = "umap",
        group.by = 'banksy_cluster',
        label=T,
        label.size=6)

plot_umap(xenium.obj,resolution ='banksy_cluster', 
          reduction='umap.banksy',directory=figs_folder)
umap_clustering_qc(xenium.obj,resolution ='banksy_cluster', 
                   reduction='umap.banksy',directory=figs_folder)


### banksy spatial plots ----
ImageDimPlot(xenium.obj, cols = "polychrome", size = 1.5,
             group.by = 'banksy_cluster',dark.background = FALSE)

ImageDimPlot(xenium.obj, cols = "polychrome", size = 1.5,
             group.by = 'banksy_cluster',dark.background = TRUE)

cropped.coords <- Crop(xenium.obj[["fov"]], x = c(4250, 4500), y = c(3000, 3200), coords = "plot")
xenium.obj[["zoom"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentations"
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", 
             border.size = 0.1, cols = "polychrome",coord.fixed = FALSE, 
             molecules = c('EPCAM','ADIPOQ','PECAM1','CD8A'), nmols = 10000,
             group.by = 'banksy_cluster', dark.background=TRUE)

DefaultAssay(xenium.obj) = 'SCT'
ImageFeaturePlot(xenium.obj, features = c('EPCAM','ADIPOQ','PECAM1','CD8A'), 
                 max.cutoff = rep('q50',4), 
                 size = 0.75, fov = 'zoom',
                 cols = c('white','red'))


SpatialDimPlot(xenium.obj, group.by = "SCT_snn_res.0.8", 
               label = T, 
               repel = T, 
               label.size = 4)


































## clustering with banksy and locations ----
xenium.obj = NormalizeData(xenium.obj, scale.factor=100,
                           normalization.method='RC',assay='Xenium')

# xenium.obj <- RunBanksy(xenium.obj, lambda = 0.2, verbose=TRUE, dimx='x',dimy='y',
#                         assay = 'Xenium', slot = 'data', features = 'all',
#                         k_geom = 50, assay_name = 'BANKSY_loc')

xenium.obj <- RunBanksy(xenium.obj, lambda = 0.8, verbose=TRUE, dimx='x',dimy='y',
                        assay = 'Xenium', slot = 'data', features = 'all',
                        k_geom = 30, assay_name = 'BANKSY_loc')

xenium.obj <- RunPCA(xenium.obj, 
                     assay='BANKSY_loc',
                     npcs = 50,
                     reduction.name='pca.banksy_loc',
                     features = rownames(xenium.obj))


ElbowPlot(object = xenium.obj, reduction='pca.banksy_loc', 
          ndims = 50)+
  ggtitle('Elbow plot of PCs')+
  theme_classic2()+
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))
ggsave(paste0(figs_folder,'elbow_plot_banksy_loc.png'), dpi=320, width=20, height= 10)

### Determine percent of variation associated with each PC
pct = xenium.obj[["pca.banksy_loc"]]@stdev / sum(xenium.obj[["pca.banksy_loc"]]@stdev) * 100
pct
#
# # Calculate cumulative percents for each PC
cumu = cumsum(pct)
cumu
#
# # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 = which(cumu >= 90 & pct <= 5)[1]
co1 # 43
#
# # Determine the difference between variation of PC and subsequent PC.
# # This is the last point where change of % of variation is more than 0.1%. Afterwards, % of variation change is < 0.1%
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2 # 13
#
# # Minimum of the two calculations
pcs = min(co1, co2)
pcs 
# [1] 13


xenium.obj <- FindNeighbors(xenium.obj, 
                            reduction = "pca.banksy_loc", 
                            dims = 1:20)

xenium.obj <- FindClusters(xenium.obj, 
                           cluster.name = 'banksy_loc_cluster',
                           resolution = 0.6,
                           graph.name='BANKSY_loc_snn')


xenium.obj <- RunUMAP(xenium.obj, 
                      reduction = 'pca.banksy_loc', 
                      reduction.name='umap.banksy_loc',
                      dims = 1:20,
)

saveRDS(xenium.obj, paste0(data_folder,'xenium_sct_clustered_seurat_banksy_loc.rds'))


# SCT_snn_res.0.5 clusters: 13


### UMAPs banksy with locations ----
# Idents(xenium.obj) = 'SCT_snn_res.1.4'
# Idents(xenium.obj) = 'SCT_snn_res.1.2'
  DimPlot(xenium.obj, cols= hue_pal()(14), order=13:0,
        reduction = "umap.banksy_loc",
        group.by = 'banksy_loc_cluster',
        label=T,
        label.size=6)

DefaultAssay(xenium.obj) = 'Xenium'
# DefaultAssay(xenium.obj) = 'SCT'
FeaturePlot(xenium.obj, features = c('CD8A'),
            reduction='umap.banksy_loc',
            min.cutoff='q50',
            max.cutoff='q90',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(xenium.obj, features = c('PECAM1'),
        assay='Xenium')




plot_umap(xenium.obj,resolution ='banksy_loc_cluster', 
          reduction='umap.banksy_loc',directory=figs_folder)
umap_clustering_qc(xenium.obj,resolution ='banksy_loc_cluster', 
                   reduction='umap.banksy_loc',directory=figs_folder)





## Viz spatial domain niche clustering ----
FeatureScatter(xenium.obj, 'x', 'y', pt.size = 0.75, cols= hue_pal()(14))+
  guides(color = guide_legend(override.aes = list(size = 3)))

FeatureScatter(xenium.obj, 'x', 'y', pt.size = 0.75, cols= hue_pal()(14))+
  guides(color = guide_legend(override.aes = list(size = 3)))


FeatureScatter(xenium.obj, 'x', 'y', cols= hue_pal()(14), pt.size = 0.1) + facet_wrap(~ colors)+
  guides(color = guide_legend(override.aes = list(size = 3)))



# annotate clusters ----
# will annotate at a resolution of 0.8 and subcluster as needed
# Idents(xenium.obj) = 'SCT_snn_res.0.8' # set cluster #s as ident
# 
Idents(xenium.obj) = 'SCT_snn_res.1.4' # set cluster #s as ident

# Idents(xenium.obj) = 'banksy_cluster' # set cluster #s as ident

DefaultAssay(xenium.obj) = "Xenium"

## Find marker genes ----

## Find all positive markers ----
# join layers to be able to run FindAllMarkers()
# combined_seurat = JoinLayers(xenium.obj)

all.pos.markers = FindAllMarkers(xenium.obj,
                                 only.pos=TRUE,
                                 logfc.threshold = 0.25,
                                 min.pct=0.1,
                                 min.diff.pct = -Inf) 

# Combine markers with gene descriptions 
annotations = read.csv('C:/Users/david/Documents/Research/PhD/Spatial_Omics/Xenium/annotations/annotationhub_human_genes.csv')
ann.pos.markers = left_join(x = all.pos.markers, 
                            y = annotations[, c("gene_name", "description")],
                            by = c("gene" = "gene_name")) %>%
  unique()

# add pct.diff col
ann.pos.markers$pct.diff = ann.pos.markers$pct.1 - ann.pos.markers$pct.2

# Rearrange the columns to be more intuitive
ann.pos.markers = ann.pos.markers[ , c(6, 7, 2:4, 9, 1, 5,8)]

# Order the rows by p-adjusted values
ann.pos.markers = ann.pos.markers %>%
  dplyr::arrange(cluster, p_val_adj)

# save markers
# write.csv(ann.pos.markers,paste0(output_annot_dir,'annotated_pos_markers_res.0.8_df.csv'),row.names=FALSE)
write.csv(ann.pos.markers,paste0(annotations_folder,'annotated_pos_markers_res.1.4.csv'),row.names=FALSE)
# write.csv(ann.pos.markers,paste0(annotations_folder,'annotated_pos_markers_res.banksy_cluster.csv'),row.names=FALSE)


# Extract top 20 markers per cluster
top = ann.pos.markers %>% 
  group_by(cluster) %>%
  slice_max(order_by=tibble(avg_log2FC,pct.diff,p_val_adj),n=20)
write.csv(ann.pos.markers,paste0(annotations_folder,'top20_annotated_pos_markers_res.1.4.csv'),row.names=FALSE)

rm(combined_seurat) # free up RAM
gc()



# celltype assignments ----
DefaultAssay(xenium.obj) = 'Xenium'
Idents(xenium.obj) = xenium.obj$SCT_snn_res.1.4

xenium.obj = NormalizeData(xenium.obj, assay='Xenium', normalization.method='LogNormalize')
xenium.obj = NormalizeData(xenium.obj, assay='Xenium', normalization.method='RC', scale.factor = 1e6)


xenium.obj = ScaleData(xenium.obj)
# DoHeatmap(xenium.obj,
#           assay='SCT',
#           # features = VariableFeatures(xenium.obj,assay='SCT')[50:100], 
#           size = 2.5,
#           angle = 90, 
#           group.by='SCT_snn_res.1.4')
# ggsave(paste0(annotations_folder,'heheatmap.png'), dpi=320, width=30, height= 10)


features = c('TGFB1','POSTN','ALDH1A1', 
             'MMP12','MMP1','MMP2',
             'PDGFRA','CXCL12','ACTA2') # fibroblasts
# features = c('CCN1','MAL2','CDH1','ADLH1A3')
features=c('CD24',)

VlnPlot(xenium.obj, features = features, ncol = 3,slot = 'data',pt.size = 0)
ggsave(paste0(annotations_folder,'vlnplot.png'), dpi=320, width=20, height= 20)

FeaturePlot(xenium.obj, features = features,
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()
ggsave(paste0(annotations_folder,'featureplot.png'), dpi=320, width=20, height= 20)

DotPlot(xenium.obj, features = features,assay='Xenium') + RotatedAxis()
ggsave(paste0(annotations_folder,'dotplot.png'), dpi=320, width=20, height= 10)

DoHeatmap(xenium.obj, features=features,size=3,assay='SCT')
ggsave(paste0(annotations_folder,'heatmap.png'), dpi=320, width=20, height= 20)



## Assign major cell types ----
xenium.obj = readRDS(paste0(data_folder,'annotated_xenium_sctype.rds'))

Idents(xenium.obj) = xenium.obj$SCT_snn_res.1.4
xenium.obj = RenameIdents(xenium.obj, 
                          "0" = "Cancer-associated fibroblasts", # Markers such as LUM (lumican), DCN (decorin), POSTN (periostin), FN1 (fibronectin), ACTA2 (smooth muscle actin), PDGFRB, and MYH11 indicate a stromal, contractile, ECM‐producing fibroblast phenotype.
                          "1" = "Tumor epithelial", # The presence of keratins (KRT23, KRT14, KRT7, KRT5) and basal markers (TAGLN, CSRP1) along with TGFB2 suggest basal-like epithelial tumor cells.
                          "2" = "Tumor associated macrophages", # Expression of myeloid markers such as LYZ, CD68, CD14, FCER1G, ITGAX, and complement components (C3, C1QC) points to macrophage/dendritic cell lineage.
                          "3" = "Tumor epithelial", # Epithelial markers (EPCAM, various keratins, DSC2, CD24) combined with proliferation indicators (CENPF, SMS) suggest proliferative basal-like tumor cells.
                          "4" = "Tumor epithelial", # Markers such as MAL2, CDH1, KRT8, EPCAM, ALDH1A3, and luminal-associated genes (TACSTD2, SCD) indicate a luminal epithelial phenotype.
                          "5" = "Cancer-associated fibroblasts", # The expression of extracellular matrix genes (DCN, LUM, POSTN, FN1), PDGF receptors (PDGFRA, PDGFRB), and remodeling inhibitors (TIMP3) is typical of activated fibroblasts in the tumor stroma.
                          "6" = "Tumor epithelial", # High levels of EPCAM, CD24, MAL2, KRT7, and luminal regulators (GATA3, ERBB2) point to a luminal/secretory epithelial cell type.
                          "7" = "Tumor epithelial", # In addition to epithelial markers (EPCAM, CD24, CDH1), the expression of KIT along with secretory markers (LTF, SCGB2A1) defines a KIT-positive epithelial subset.
                          "8" = "Tumor associated macrophages", # Co-expression of myeloid markers (CD68, CD14, ITGAX, ITGAM, LYZ, CD163) and immune regulators (HAVCR2, FCER1G) indicates a macrophage population.
                          "9" = "Endothelial", # Classic endothelial markers such as PECAM1, VWF, AQP1, KDR, ANGPT2, and CLEC14A clearly mark vascular endothelial cells.
                          "10" = "Tumor epithelial", # The combination of multiple keratins (KRT5, KRT7, KRT23, KRT14), epithelial adhesion molecules (EPCAM, DSC2), and basal markers (CSRP1, TAGLN) supports a basal-like epithelial phenotype.
                          "11" = "Tumor epithelial", # Epithelial markers (EPCAM, CD24, DSC2, KRT7) with luminal transcription factors (GATA3, TFAP2A) and proliferation marker (CCND1) indicate a luminal-like epithelial subtype.
                          "12" = "Cancer-associated fibroblasts", # High expression of ECM genes (DCN, LUM, FBLN1), growth factor receptors (PDGFRA), and matrix remodeling enzymes (MMP2) suggests matrix-producing fibroblasts in the stroma.
                          "13" = "Myoepithelial", # Expression of smooth muscle/myoepithelial markers (MYH11, MYLK, ACTA2, TAGLN, DSP) alongside basal keratins (KRT14, KRT5, KRT7) is indicative of myoepithelial cells.
                          "14" = "Cancer-associated fibroblasts", # The co‐expression of contractile proteins (ACTA2, TAGLN), ECM components (FN1, POSTN), and epithelial markers (KRT14, KRT5) points to a myoepithelial phenotype.
                          "15" = "Cytoxic T cells", # T cell receptor components (CD3E, TRAC, CD3G, CD247) together with cytotoxic molecules (GZMA, PRF1, NKG7) and co-stimulatory markers (TIGIT, SLAMF1) denote cytotoxic T cells.
                          "16" = "Adipocytes", # Markers such as ADIPOQ, LPL, PPARG, SCD, and metabolic enzymes are characteristic of adipocyte lineage and fat metabolism.
                          "17" = "Tumor associated macrophages", # High expression of M2 macrophage markers (MRC1, CD163, CD14, CD68, AIF1, ITGAM, LYZ) along with IL2RA and PDK4 suggests an immunosuppressive, M2-like macrophage phenotype.
                          "18" = "Tumor epithelial", # Expression of luminal markers including FOXA1, XBP1, CD24, KRT7, ERBB2, and GATA3 is typical of luminal epithelial cells in breast tissue.
                          "19" = "Tumor epithelial", # The presence of SPP1, TACSTD2, TGFB2, TPD52, MAL2, CD24, SFRP1, CCN2, SMS, and basal keratins (KRT14, KRT7) point to a basal-like epithelial phenotype.
                          "20" = "Cancer-associated fibroblasts", # ECM remodeling enzymes (MMP1, MMP14, MMP2), matrix components (LUM, DCN, POSTN), and PDGF receptors (C1R, PDGFRB) are indicative of activated fibroblasts in the tumor stroma.
                          "21" = "Tumor epithelial", # Epithelial markers (CD24, DSC2, SYNE2, KRT7) with luminal features (ESR1, TFAP2A, KLF5) and KIT expression suggest a luminal epithelial cell population.
                          "22" = "Smooth muscle", # The co-expression of smooth muscle genes (MYH11, ACTA2, MYLK, ACTG2, TAGLN) along with pericyte markers (PDGFRB, ANGPT2, CXCL12) defines smooth muscle cells or pericytes.
                          "23" = "Endothelial", # Markers such as AQP1, VWF, MMRN2, PECAM1, and CLEC14A are definitive for endothelial cells lining blood vessels.
                          "24" = "Mast cells", # The presence of CPA3, TPSAB1, and CTSG are classic mast cell markers; KIT is also highly expressed in mast cells.
                          "25" = "Endothelial", # Endothelial markers like RAMP2, PECAM1, and CAVIN2, along with FOXC2 and EGFL7, define a specialized vascular endothelial subset.
                          "26" = "Cancer-associated fibroblasts", # Expression of ECM-related genes (SFRP4, DPT, FBLN1, CCDC80, LUM) and growth factors (IGF1, PDGFRA) supports a matrix-producing fibroblast identity.
                          "27" = "Tumor associated macrophages", # The high expression of MMP12, ITGAX, CD68, CD163, ITGAM, and LYZ, along with other myeloid markers, identifies a macrophage population within the tumor.
                          "28" = "Plasma cells" # B/plasma cell markers such as MZB1, CD79A, TENT5C, SLAMF7, TNFRSF17, PRDM1, and XBP1 are hallmarks of plasma cells in the tumor microenvironment.
                          )

xenium.obj$cell_type = Idents(xenium.obj)

library(presto)

clusters = c(0,5,12,14,20,26) # CAFs
clusters = c(1,3,4,6,7,10,11,18,19,21) # tumor epithelial
clusters = c(2,8,17,27) # macrophages
clusters = c(9,23,25) # endothelial
clusters = c()

for (clust in clusters) {
  cat('\nTesting cluster',clust, 'against clusters:',clusters[clusters != clust],'\n\n')
  test = FindMarkers(xenium.obj,ident.1 = clust, ident.2 = clusters[clusters != clust],
            group.by = 'SCT_snn_res.1.4',
            assay = 'Xenium',
            slot='counts',
            logfc.threshold = 1,
            only.pos=TRUE) %>%
  filter(p_val_adj < 0.05)
  test$pct.diff = test$pct.1 - test$pct.2
  test = test[ , c(2,6,3:4,6,1,5)]
  print(test)
}
test2 = FindMarkers(xenium.obj,ident.1 = 19, ident.2 = 10,
                   group.by = 'SCT_snn_res.1.4',
                   assay = 'Xenium',
                   slot='counts',
                   logfc.threshold = 0.25,
                   only.pos=F) %>%
  filter(p_val_adj < 0.05)
test2$pct.diff = test2$pct.1 - test2$pct.2
test2 = test2[ , c(2,6,3:4,6,1,5)]
test2


## Assign cell subtypes, long version ----
Idents(xenium.obj) = xenium.obj$SCT_snn_res.1.4
xenium.obj = RenameIdents(xenium.obj, 
                          "0" = "CAFs - matrix-producing", # Markers such as LUM (lumican), DCN (decorin), POSTN (periostin), FN1 (fibronectin), ACTA2 (smooth muscle actin), PDGFRB, and MYH11 indicate a stromal, contractile, ECM‐producing fibroblast phenotype.
                          "1" = "Tumor epithelial - luminal-basal-like angiogenic", # The presence of keratins (KRT23, KRT14, KRT7, KRT5) and basal markers (TAGLN, CSRP1) along with TGFB2 suggest basal-like epithelial tumor cells.
                          "2" = "TAMs - mixed M1 M2 tumor infiltrating", # Expression of myeloid markers such as LYZ, CD68, CD14, FCER1G, ITGAX, and complement components (C3, C1QC) points to macrophage/dendritic cell lineage.
                          "3" = "Tumor epithelilal - basal-luminal-like proliferating", # Epithelial markers (EPCAM, various keratins, DSC2, CD24) combined with proliferation indicators (CENPF, SMS) suggest proliferative basal-like tumor cells.
                          "4" = "Tumor epithelial - luminal-basal-like EMT", # Markers such as MAL2, CDH1, KRT8, EPCAM, ALDH1A3, and luminal-associated genes (TACSTD2, SCD) indicate a luminal epithelial phenotype. MAL2 mediates T cell evasion
                          "5" = "CAFs - matrix-remodeling.invasive", # The expression of extracellular matrix genes (DCN, LUM, POSTN, FN1), PDGF receptors (PDGFRA, PDGFRB), and remodeling inhibitors (TIMP3) is typical of activated fibroblasts in the tumor stroma.
                          "6" = "Tumor epithelial - luminal-secretory-like", # High levels of EPCAM, CD24, MAL2, KRT7, and luminal regulators (GATA3, ERBB2) point to a luminal/secretory epithelial cell type. Luminal-A like?
                          "7" = "Tumor epithelial - luminal-progenitor-like HER2+", # In addition to epithelial markers (EPCAM, CD24, CDH1), the expression of KIT along with secretory markers (LTF, SCGB2A1) defines a KIT-positive epithelial subset.
                          "8" = "TAMs - fibrosis-associated", # Co-expression of myeloid markers (CD68, CD14, ITGAX, ITGAM, LYZ, CD163) and immune regulators (HAVCR2, FCER1G) indicates a macrophage population.
                          "9" = "Endothelial - periepithelial", # Classic endothelial markers such as PECAM1, VWF, AQP1, KDR, ANGPT2, and CLEC14A clearly mark vascular endothelial cells.
                          "10" = "Tumor epithelial - basal-like SERPINA3.PIGR+", # The combination of multiple keratins (KRT5, KRT7, KRT23, KRT14), epithelial adhesion molecules (EPCAM, DSC2), and basal markers (CSRP1, TAGLN) supports a basal-like epithelial phenotype.
                          "11" = "Tumor epithelial - luminal-like HER2.GATA3+)", # Epithelial markers (EPCAM, CD24, DSC2, KRT7) with luminal transcription factors (GATA3, TFAP2A) and proliferation marker (CCND1) indicate a luminal-like epithelial subtype.
                          "12" = "CAFs - stromal - angiogenic", # High expression of ECM genes (DCN, LUM, FBLN1), growth factor receptors (PDGFRA), and matrix remodeling enzymes (MMP2) suggests matrix-producing fibroblasts in the stroma.
                          "13" = "Myoepithelial", # Expression of smooth muscle/myoepithelial markers (MYH11, MYLK, ACTA2, TAGLN, DSP) alongside basal keratins (KRT14, KRT5, KRT7) is indicative of myoepithelial cells.
                          "14" = "CAFs - Matrix stiffening - TGFB signaling", # The co‐expression of contractile proteins (ACTA2, TAGLN), ECM components (FN1, POSTN), and epithelial markers (KRT14, KRT5) points to a myoepithelial phenotype.
                          "15" = "Cytotoxic T cells", # T cell receptor components (CD3E, TRAC, CD3G, CD247) together with cytotoxic molecules (GZMA, PRF1, NKG7) and co-stimulatory markers (TIGIT, SLAMF1) denote cytotoxic T cells.
                          "16" = "Adipocytes", # Markers such as ADIPOQ, LPL, PPARG, SCD, and metabolic enzymes are characteristic of adipocyte lineage and fat metabolism.
                          "17" = "TAMs - adipose tissue resident - M2-like", # High expression of M2 macrophage markers (MRC1, CD163, CD14, CD68, AIF1, ITGAM, LYZ) along with IL2RA and PDK4 suggests an immunosuppressive, M2-like macrophage phenotype.
                          "18" = "Tumor epithelial cells - luminal-like FOXA1+)", # Expression of luminal markers including FOXA1, XBP1, CD24, KRT7, ERBB2, and GATA3 is typical of luminal epithelial cells in breast tissue.
                          "19" = "Tumor epithelial cells - basal-like fibrogenic", # The presence of SPP1, TACSTD2, TGFB2, TPD52, MAL2, CD24, SFRP1, CCN2, SMS, and basal keratins (KRT14, KRT7) point to a basal-like epithelial phenotype.
                          "20" = "CAFs - matrix remodeling", # ECM remodeling enzymes (MMP1, MMP14, MMP2), matrix components (LUM, DCN, POSTN), and PDGF receptors (C1R, PDGFRB) are indicative of activated fibroblasts in the tumor stroma.
                          "21" = "Tumor epithelial cells - luminal-secratory-like SERPINA3.PIGR+)", # Epithelial markers (CD24, DSC2, SYNE2, KRT7) with luminal features (ESR1, TFAP2A, KLF5) and KIT expression suggest a luminal epithelial cell population.
                          "22" = "Smooth muscle", # The co-expression of smooth muscle genes (MYH11, ACTA2, MYLK, ACTG2, TAGLN) along with pericyte markers (PDGFRB, ANGPT2, CXCL12) defines smooth muscle cells or pericytes.
                          "23" = "Endothelial - adipocyte-associated", # metabolic. Markers such as AQP1, VWF, MMRN2, PECAM1, and CLEC14A are definitive for endothelial cells lining blood vessels.
                          "24" = "Mast cells", # The presence of CPA3, TPSAB1, and CTSG are classic mast cell markers; KIT is also highly expressed in mast cells.
                          "25" = "Endothelial cells - vascular", # Endothelial markers like RAMP2, PECAM1, and CAVIN2, along with FOXC2 and EGFL7, define a specialized vascular endothelial subset.
                          "26" = "CAFs - stromal - metabolic-immune", # Expression of ECM-related genes (SFRP4, DPT, FBLN1, CCDC80, LUM) and growth factors (IGF1, PDGFRA) supports a matrix-producing fibroblast identity.
                          "27" = "TAMs - periluminal ", # The high expression of MMP12, ITGAX, CD68, CD163, ITGAM, and LYZ, along with other myeloid markers, identifies a macrophage population within the tumor.
                          "28" = "Plasma cells" # B/plasma cell markers such as MZB1, CD79A, TENT5C, SLAMF7, TNFRSF17, PRDM1, and XBP1 are hallmarks of plasma cells in the tumor microenvironment.
)

xenium.obj$cell_subtype_long = Idents(xenium.obj)

xenium.obj$cell_subtype = Idents(xenium.obj)


## Assign cell subtypes, short version ----
Idents(xenium.obj) = xenium.obj$SCT_snn_res.1.4
xenium.obj = RenameIdents(xenium.obj,
                          "0" = "mCAFs", # Markers such as LUM (lumican), DCN (decorin), POSTN (periostin), FN1 (fibronectin), ACTA2 (smooth muscle actin), PDGFRB, and MYH11 indicate a stromal, contractile, ECM‐producing fibroblast phenotype.
                          "1" = "Tumor epithelial - luminal-basal-like angiogenic", # The presence of keratins (KRT23, KRT14, KRT7, KRT5) and basal markers (TAGLN, CSRP1) along with TGFB2 suggest basal-like epithelial tumor cells.
                          "2" = "TAMs - mixed M1 M2 tumor infiltrating", # Expression of myeloid markers such as LYZ, CD68, CD14, FCER1G, ITGAX, and complement components (C3, C1QC) points to macrophage/dendritic cell lineage.
                          "3" = "Tumor epithelilal - basal-luminal-like proliferating", # Epithelial markers (EPCAM, various keratins, DSC2, CD24) combined with proliferation indicators (CENPF, SMS) suggest proliferative basal-like tumor cells.
                          "4" = "Tumor epithelial - luminal-basal-like EMT", # Markers such as MAL2, CDH1, KRT8, EPCAM, ALDH1A3, and luminal-associated genes (TACSTD2, SCD) indicate a luminal epithelial phenotype. MAL2 mediates T cell evasion
                          "5" = "CAFs - matrix-remodeling.invasive", # The expression of extracellular matrix genes (DCN, LUM, POSTN, FN1), PDGF receptors (PDGFRA, PDGFRB), and remodeling inhibitors (TIMP3) is typical of activated fibroblasts in the tumor stroma.
                          "6" = "Tumor epithelial - luminal-secretory-like", # High levels of EPCAM, CD24, MAL2, KRT7, and luminal regulators (GATA3, ERBB2) point to a luminal/secretory epithelial cell type. Luminal-A like?
                          "7" = "Tumor epithelial - luminal-progenitor-like HER2+", # In addition to epithelial markers (EPCAM, CD24, CDH1), the expression of KIT along with secretory markers (LTF, SCGB2A1) defines a KIT-positive epithelial subset.
                          "8" = "TAMs - fibrosis-associated", # Co-expression of myeloid markers (CD68, CD14, ITGAX, ITGAM, LYZ, CD163) and immune regulators (HAVCR2, FCER1G) indicates a macrophage population.
                          "9" = "Endothelial - periepithelial", # Classic endothelial markers such as PECAM1, VWF, AQP1, KDR, ANGPT2, and CLEC14A clearly mark vascular endothelial cells.
                          "10" = "Tumor epithelial - basal-like metaplastic link", # The combination of multiple keratins (KRT5, KRT7, KRT23, KRT14), epithelial adhesion molecules (EPCAM, DSC2), and basal markers (CSRP1, TAGLN) supports a basal-like epithelial phenotype.
                          "11" = "Tumor epithelial - luminal-like HER2.GATA3+)", # Epithelial markers (EPCAM, CD24, DSC2, KRT7) with luminal transcription factors (GATA3, TFAP2A) and proliferation marker (CCND1) indicate a luminal-like epithelial subtype.
                          "12" = "CAFs - stromal - angiogenic", # High expression of ECM genes (DCN, LUM, FBLN1), growth factor receptors (PDGFRA), and matrix remodeling enzymes (MMP2) suggests matrix-producing fibroblasts in the stroma.
                          "13" = "Myoepithelial cells", # Expression of smooth muscle/myoepithelial markers (MYH11, MYLK, ACTA2, TAGLN, DSP) alongside basal keratins (KRT14, KRT5, KRT7) is indicative of myoepithelial cells.
                          "14" = "CAFs - Matrix stiffening - TGFB signaling", # The co‐expression of contractile proteins (ACTA2, TAGLN), ECM components (FN1, POSTN), and epithelial markers (KRT14, KRT5) points to a myoepithelial phenotype.
                          "15" = "Cytotoxic T cells", # T cell receptor components (CD3E, TRAC, CD3G, CD247) together with cytotoxic molecules (GZMA, PRF1, NKG7) and co-stimulatory markers (TIGIT, SLAMF1) denote cytotoxic T cells.
                          "16" = "Adipocytes", # Markers such as ADIPOQ, LPL, PPARG, SCD, and metabolic enzymes are characteristic of adipocyte lineage and fat metabolism.
                          "17" = "TAMs - adipose tissue resident - M2-like", # High expression of M2 macrophage markers (MRC1, CD163, CD14, CD68, AIF1, ITGAM, LYZ) along with IL2RA and PDK4 suggests an immunosuppressive, M2-like macrophage phenotype.
                          "18" = "Tumor epithelial cells - luminal-hormone-sensing-like FOXA1+)", # Expression of luminal markers including FOXA1, XBP1, CD24, KRT7, ERBB2, and GATA3 is typical of luminal epithelial cells in breast tissue.
                          "19" = "Tumor epithelial cells - basal-like", # The presence of SPP1, TACSTD2, TGFB2, TPD52, MAL2, CD24, SFRP1, CCN2, SMS, and basal keratins (KRT14, KRT7) point to a basal-like epithelial phenotype.
                          "20" = "CAFs - Matrix remodeling", # ECM remodeling enzymes (MMP1, MMP14, MMP2), matrix components (LUM, DCN, POSTN), and PDGF receptors (C1R, PDGFRB) are indicative of activated fibroblasts in the tumor stroma.
                          "21" = "Tumor epithelial cells - luminal-secratory like SERPINA3.PIGR+)", # Epithelial markers (CD24, DSC2, SYNE2, KRT7) with luminal features (ESR1, TFAP2A, KLF5) and KIT expression suggest a luminal epithelial cell population.
                          "22" = "Smooth muscle", # The co-expression of smooth muscle genes (MYH11, ACTA2, MYLK, ACTG2, TAGLN) along with pericyte markers (PDGFRB, ANGPT2, CXCL12) defines smooth muscle cells or pericytes.
                          "23" = "Endothelial - adipocyte associated - metabolic", # Markers such as AQP1, VWF, MMRN2, PECAM1, and CLEC14A are definitive for endothelial cells lining blood vessels.
                          "24" = "Mast cells", # The presence of CPA3, TPSAB1, and CTSG are classic mast cell markers; KIT is also highly expressed in mast cells.
                          "25" = "Endothelial cells - vascular", # Endothelial markers like RAMP2, PECAM1, and CAVIN2, along with FOXC2 and EGFL7, define a specialized vascular endothelial subset.
                          "26" = "iCAFs", # Expression of ECM-related genes (SFRP4, DPT, FBLN1, CCDC80, LUM) and growth factors (IGF1, PDGFRA) supports a matrix-producing fibroblast identity.
                          "27" = "TAMs - periluminal - TMEM", # The high expression of MMP12, ITGAX, CD68, CD163, ITGAM, and LYZ, along with other myeloid markers, identifies a macrophage population within the tumor.
                          "28" = "Plasma cells", # B/plasma cell markers such as MZB1, CD79A, TENT5C, SLAMF7, TNFRSF17, PRDM1, and XBP1 are hallmarks of plasma cells in the tumor microenvironment.
)

# xenium.obj$cell_subtype = Idents(xenium.obj)


# paste(shQuote(colnames(xenium.obj@meta.data),type='sh'), collapse=", ")

# export_metadata = xenium.obj@meta.data[,c('orig.ident', 'nCount_Xenium', 'nFeature_Xenium', 
#                        'nCount_BlankCodeword', 'nFeature_BlankCodeword', 
#                        'nCount_ControlCodeword', 'nFeature_ControlCodeword', 
#                        'nCount_ControlProbe', 'nFeature_ControlProbe', 
#                        'cells', 'x', 'y', 'complexity',  'SCT_snn_res.1.4',
#                        'seurat_clusters', 'cell_type', 'cell_subtype')] # need to add long subtype later if used
# 
# export = DietSeurat(xenium.obj, 
#                     layers=c('counts'),
#                     assays=c('Xenium',"BlankCodeword",'ControlCodeword','ControlProbe'),
#                     dimreducs=c('pca','umap'),
#                     graphs=c('SCT_nn','SCT_snn')
# )
# export@meta.data = export_metadata
# 
# saveRDS(export, paste0(data_folder,'xenium_annot.rds'))


table(xenium.obj$cell_type)
table(xenium.obj$cell_subtype)


saveRDS(xenium.obj, paste0(data_folder,'xenium_annotated.rds'))


### create csv for xenium explorer
df_clusters = data.frame(
  cell_id = names(xenium.obj$SCT_snn_res.1.4),
  group = xenium.obj$SCT_snn_res.1.4
)

df_celltype = data.frame(
  cell_id = names(xenium.obj$cell_type),
  group = xenium.obj$cell_type
)

df_cellsubtype = data.frame(
  cell_id = names(xenium.obj$cell_subtype),
  group = xenium.obj$cell_subtype
)



write.csv(df_clusters, paste0(annotations_folder, 'custom_cell_groups_SCT_snn_res.1.4.csv'))
write.csv(df_celltype, paste0(annotations_folder, 'custom_cell_groups_cell_type.csv'))
write.csv(df_cellsubtype, paste0(annotations_folder, 'custom_cell_groups_cell_subtype.csv'))


