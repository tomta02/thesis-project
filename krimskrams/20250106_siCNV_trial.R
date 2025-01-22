# 2025-01-06
# trying out the tool "spatial infer CNV" (siCNV):
# (https://github.com/aerickso/SpatialInferCNV)
# on one of the visium samples 
#============================================================================#

# importing necessary packagesm
library(infercnv)
library(tidyverse)
library(Seurat)
library(phylogram)
library(ape)
library(hdf5r)
library(SpatialInferCNV)

options(scipen = 100)

#============================================================================#
# Importing data for running siCNV                                                  
#============================================================================#

# loading count data from Visium sample, which we had saved as tsv
counts <- read.delim("pat1_R3_counts.tsv", sep = "\t")

# check
counts[1:5, 1:5] 


# prepend patient nr to rownames
rownames(counts) <- paste0("patient1", rownames(counts))

# the "ImportHistologicalAnnotations" function from the siCNV package which we use
# later will import barcode names such that "-" gets replaced by ".". 
# Therefore adjust the barcode names in the count data, so that they match later
rownames(counts) <- gsub("-", ".", rownames(counts))
counts$Barcode <- rownames(counts)

# importing annotation file for the two clusters we want to compare 
# file needed: csv file with 2 rows, first row: Visium spot ID, second row:
# histological annotation of spot (eg. "healthy", "cancer")
# in our case: hist annot. = cluster identity (either "0", "2")
cluster_0_annot <- ImportHistologicalAnnotations("patient1", 
                        "20241212_pat1_r3_annot_file_cluster0.tsv")
cluster_2_annot <- ImportHistologicalAnnotations("patient1", 
                        "20241212_pat1_r3_annot_file_cluster2.tsv")

# edit: inferCNV function gives some weird error because our spot annotations are 
# numbers (from clustering), not strings. Replace clust "0" -> "diffuse", "2" -> "follicular"
cluster_0_annot <- cluster_0_annot %>%
  mutate(Histology = ifelse(Histology == 0, "diffuse", Histology))

cluster_2_annot <- cluster_2_annot %>%
  mutate(Histology = ifelse(Histology == 2, "follicular", Histology))

# append cluster annotations together
cluster_annot <- rbind(cluster_0_annot, cluster_2_annot)



# Now we have count data, as well as the histological annotations that
# we want to compute CNVs for. Last thing we need is "GeneOrderFile"
# is essentially a .tsv file containing 4 cols: 
# gene name, chromosome #, start, stop (i.e. pos of gene on chromosom)
# we can import one from the Trinity CTAT project for hg38

#url <- "https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt"
#dest_file <- "hg38_gencode_v27.txt"  

# download file
#download.file(url, destfile = dest_file)

# loading the gene order file 
geneorderfile_txt <- read.table("hg38_gencode_v27.txt", 
                            sep = "\t", 
                            header = FALSE)



#============================================================================#
# preparations before running siCNV                                                
#============================================================================#
# Before we can run siCNV:
# 1.) use siCNVs MergingCountAndAnnotationData" func
#   -> simply subsets "counts" df (i.e. spot vs gene matrix of whole slide) 
#      by the barcodes that will go into the analysis
#      i.e. barcodes that we have in "cluster_annot"   
# 2.) create an "inferCNV object" to run siCNV on

# merge read counts with the spot barcode information:
joined_counts_annot <- MergingCountAndAnnotationData("patient1",
                                              cluster_annot, 
                                              counts)
# setting "Genes" column as the rownames, instead of separate colum
joined_counts_annot <- joined_counts_annot %>% column_to_rownames(var = "Genes")

# save counts to directory
write.table(joined_counts_annot, "./pat1_r3_joined_counts.tsv",
            sep = "\t")

# save the histological annotations to directory
write.table(cluster_annot, "./pat1_r3_histo_annot.tsv", 
            sep = "\t",
            quote = FALSE, 
            col.names = FALSE, 
            row.names = FALSE)


# create the inferCNV object from the three files
infCNV_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix = "pat1_r3_joined_counts.tsv",
                                                   gene_order_file = "hg38_gencode_v27.txt",
                                                   annotations_file = "pat1_r3_histo_annot.tsv",
                                                   delim = "\t",
                                                   ref_group_name = "follicular")


#============================================================================#
# run siCNV                                                
#============================================================================#

# siCNV not (yet) compatible with Seurat v5, therefore add this
options("Seurat.object.assay.version" = "v3")


infCNV_run <- infercnv::run(infCNV_obj,
                            cutoff = 0.1,
                            out_dir = "./20250109_pat1_siCNV_outputs",
                            #cluster_by_groups = TRUE,
                            HMM = TRUE,
                            denoise = TRUE)
                                            



