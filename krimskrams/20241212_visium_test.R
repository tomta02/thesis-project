#============================================================================#
# 2024-12-12
# Visium trial analysis script 
# (familiarizing with visium data, differential gene expression)

#============================================================================#
# loading necessary libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# defining Seurat object
testdata = Load10X_Spatial(data.dir = '/path/to/visium/outs',
                           filename = "filtered_feature_bc_matrix.h5",
                           assay = "Spatial",
                           slice = "pat1_r3")

testdata # An object of class Seurat

#============================================================================#
# writing counts to .tsv                                                      
#============================================================================#

# before proceeding with the analysis, we want to write the gene counts 
# to a .tsv file for running tool "spatialinferCNV" (siCNV) later.
count_data <- (as.matrix(GetAssayData(object = testdata,
                                  layer = "counts")))

# check
count_data[1:5, 1:5]

write.table(count_data, "pat1_R3_counts.tsv",
            sep = "\t",
            quote = FALSE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# analysis continue

# raw data plots
# making violin and spatial feature plot as shown in seurat vignette
plot1 <- VlnPlot(testdata, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(testdata, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)


#============================================================================#
# data processing                                                           
#============================================================================#
# 1.) Normalization: using log normalization vs. scTransform? 
# read scTransform paper
testdata_norm_sct <- SCTransform(testdata, assay = "Spatial", verbose = TRUE)
testdata_lognorm <- NormalizeData(testdata, assay = "Spatial", verbose = TRUE)

# making violin plot 
# log-norm vs. scTransform-normalized data
plot3 <- VlnPlot(testdata_lognorm, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot4 <- VlnPlot(testdata_norm_sct, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
wrap_plots(plot3, plot4)

# making spatialfeature plot 
# log-norm vs. scTransform normalized data
plot5 <- SpatialFeaturePlot(testdata_lognorm, features = "nCount_Spatial") + theme(legend.position = "right")
plot6 <- SpatialFeaturePlot(testdata_norm_sct, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot5, plot6)

# proceed with scTransform-normalized data

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 2.) remove spots with count level/spot < 500

# extract normalized data from the SCT assay
normalized_data <- GetAssayData(testdata_norm_sct, assay = "SCT", layer = "data")

# calculate total RNA counts per spot
rnacounts_per_spot <- colSums(normalized_data)
head(rnacounts_per_spot)

# extract only spots wheree RNA count > 500
keepspots <- names(rnacounts_per_spot[rnacounts_per_spot >= 500])
# subset seurat object to only keep these spots
testdata_norm_subs <- subset(testdata_norm_sct, cells = keepspots)

# turns out there are no spots with RNA counts < 500 :)

#============================================================================#
# PCA                                                     
#============================================================================#

# 3.) run PCA
testdata_norm_sct <- RunPCA(testdata_norm_sct, assay = "SCT", verbose = FALSE)

# look at elbow plot to determine optimal nr of principal components to use
# to do this, calculate the variance explained by each of the PCs
ElbowPlot(testdata_norm_sct, ndims = 30, reduction = "pca")

# calculate thestandard deviations of each PC
pca_stdev <- Stdev(testdata_norm_sct, reduction = "pca")

# calculate the variance explained by each PC (in %! )
var_explained <- (pca_stdev^2) / sum(pca_stdev^2) * 100  

# calculate cumulative variance explained (by each additional PC)
cumulative_variance <- cumsum(var_explained)

# find the number of PCs that explain at least 90% of the variance in data
selected_pcs <- which(cumulative_variance >= 90)[1]  
print(selected_pcs) # 29
# at this point, just use all 30 PCs that were calculated in the first place ...

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 4.) calculate nearest neighbors, clusters; run UMAP
testdata_norm_sct <- FindNeighbors(testdata_norm_sct, reduction = "pca", dims = 1:30)

# use Lovain alg for clustering, need extra python package for Leiden?
# ask Kristy how to deal with this
# try out different resolutions, settled for resolution 0.3 to reduce the nr 
# of clusters we get. Since on each spot, we are already averaging across 
# ca. 20 - 30 cells, does not make sense to look for fine grained clusters? think about

testdata_norm_sct <- FindClusters(testdata_norm_sct, resolution = 0.3,
                              verbose = FALSE)
testdata_norm_sct <- RunUMAP(testdata_norm_sct, reduction = "pca", dims = 1:30)

# plotting
plot7 <- DimPlot(testdata_norm, reduction = "umap", label = TRUE)
plot8 <- SpatialDimPlot(testdata_norm, label = TRUE, label.size = 3)
plot7 + plot8

## do differential gene expression analysis for the 6 different clusters 
## plus look at spatially differentially expressed genes -> maybe can use
## them already for pathway analysis? Look at pathway analysis tools by Saez lab


#============================================================================#
# spatially differentially expressed genes                                                        
#============================================================================#

# two different options: 
# 1.) look for SDE genes based on previously annotated regions,
# e.g. clusters that were calculated
# 2.) look for SDE genes without previous knowledge
# try both here 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Option 1: look for SDE genes based on pre-annotated anatomical regions 
# within the tissue (= e.g. calculated clusters)

# as trial: look at differentially expr. genes in cluster 0 vs. cluster 2
# correspond to follicular (cl. 2) vs diffuse (cl. 0)
de_markers <- FindMarkers(testdata_norm_sct, ident.1 = 0, ident.2 = 2)
de_markers$gene <- rownames(de_markers)
SpatialFeaturePlot(object = testdata_norm_sct, 
                  features = rownames(de_markers)[1:6], 
                  alpha = c(0.1, 1), ncol = 3)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Option 2: look for spatially different gene expression patterns in
# absence of prior annotation. Based on 
# (https://www.nature.com/articles/nmeth.4634)

# look for differentially expressed genes
testdata_norm_sct <- FindSpatiallyVariableFeatures(testdata_norm_sct, 
                                               assay = "SCT", 
                                               features = VariableFeatures(testdata_norm)[1:1000],
                                               selection.method = "moransi")

# Note: the Seurat function "SpatiallyVariableFeatures" leads to an 
# error message here that is known when using it with Morans I alg
# followed workaround as proposed by github issue seurat #7422

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Workaround function

SpatiallyVariableFeatures_workaround <- function(object, assay="SCT", 
                                                 selection.method = "moransi") {
  # This is work around function to replace SeuratObject::SpatiallyVariableFeatures function.
  # return ranked list of Spatially Variable Features
  
  # Check if object is a Seurat object
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }
  
  # Check if assay is a valid assay
  if (!assay %in% names(object@assays)) {
    stop("assay must be a valid assay")
  }
  
  # Extract meta.features from the specified object and assay
  data <- object@assays[[assay]]@meta.features
  
  # Select columns starting with the provided col_prefix
  moransi_cols <- grep(paste0("^", selection.method), colnames(data), value = TRUE)
  
  # Filter rows where "moransi.spatially.variable" is TRUE
  filtered_data <- data[data[[paste0(selection.method, ".spatially.variable")]], moransi_cols]
  
  # Sort filtered data by "moransi.spatially.variable.rank" column in ascending order
  sorted_data <- filtered_data[order(filtered_data[[paste0(selection.method, ".spatially.variable.rank")]]), ]
  
  # Return row names of the sorted data frame
  rownames(sorted_data)
}

# calculating spatially differentially expressed genes
top_features <- head(SpatiallyVariableFeatures_workaround(testdata_norm_sct, selection.method = "moransi"), 6)
SpatialFeaturePlot(testdata_norm_sct, features = top_features, ncol = 3, alpha = c(0.1, 1))


#============================================================================#
# spatially differentially expressed genes - plots                                                       
#============================================================================#

# go ahead with SDE genes calculated using "option 1" - calc. genes based
# on predetermined clusters.

# cluster 0 -> DLBCL
# cluster 2 -> FL

# look at differentially expressed genes between cluster 0 and 2
# convert p-vals to neg. log10 p-values
de_markers$neglog10_p_val_adj <- -log10(de_markers$p_val_adj)

# determine threshold for significance 
thresh_log2FC <- 1 # =expression must be at least doubled or halved to be sig
                   # maybe too conservative, adjust?
thresh_p_val <- 0.05

# determine for each gene whether it got upregulated or downregulated 
# (if sigificant) or no significant change (not sig)
de_markers$significance <- with(de_markers, ifelse(
  (avg_log2FC > thresh_log2FC & p_val_adj < thresh_p_val), "upreg",
  ifelse(avg_log2FC < -thresh_log2FC & p_val_adj < thresh_p_val, "downreg", "not sig")
))

# make volcano plot of spatially differentially expressed genes 
# between cluster 0 and cluster 2
volc_plot <- ggplot(de_markers, aes(x = avg_log2FC, y = neglog10_p_val_adj)) +
  # color code points by sig level!
  geom_point(aes(color = significance), size = 2, alpha = 0.8) +  
  scale_color_manual(values = c("upreg" = "red", "downreg" = "blue", "not sig" = "gray")) +
  theme_minimal() +  
  labs(title = "Volcano Plot of SDE genes - cluster 0 vs. cluster 2",
       x = "log2 fold change",
       y = "-log10 adjusted p-value",
       color = "gene regulation") +
  theme(text = element_text(size = 12),
        legend.position = "top")

print(volc_plot)

# highliting most significant SDE genes
highlight_genes <- subset(de_markers, p_val_adj < 1e-5 & abs(avg_log2FC) > 2)

library(ggrepel)
volc_plot <- volc_plot + geom_text_repel(data = highlight_genes,
                  aes(label = gene),
                  size = 3,
                  max.overlaps = 10)
print(volc_plot)



#============================================================================#
# Unrelated: use defined clusters to create annotation file for running siCNV                                                     
#============================================================================#
# for siCNV, we need two annotation files containing the spots we want to
# compare w.r.t. copy number variations. Annotation file is simply .cnv 
# file containing spatial barcode + cluster they belong to.

# let's continue comparing spots belonging to cluster 0 and spots of cl. 2 
# cluster info stored in metadata
head(testdata_norm_sct@meta.data)

# make dataframe containing spatial id + cluster it belongs to
annot_df <- data.frame(
  Barcode = rownames(testdata_norm_sct@meta.data),
  cluster = testdata_norm_sct@meta.data$seurat_clusters
)

# filter annot_df for spots belonging to cluster 0 & cluster 2
cluster_0 = annot_df[annot_df$cluster == "0", ]
cluster_2 = annot_df[annot_df$cluster == "2", ]

# save the files
write.table(cluster_0, file = "20241212_pat1_r3_annot_file_cluster0.tsv",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
write.table(cluster_2, file = "20241212_pat1_r3_annot_file_cluster2.tsv",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
