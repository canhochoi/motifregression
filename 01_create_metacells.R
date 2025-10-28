library(Seurat)
library(Signac)
library(SuperCell)

#for motif information
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)


load(file = '~/luan/Juno/TF_gene_expression/TF_gene_Rusland_multiome.RData')

#remove NA cells in SCT assay----
DefaultAssay(multiome_seurat) <- "SCT"
Idents(multiome_seurat) <- "projection"

multiome_seurat_noNA <- subset(multiome_seurat, idents = c(1:10))

#integrate two models----
#follow Ruslan and this link
#https://satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette
#https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis

#check that the wnn.umap is from the previous one
all(multiome_seurat_noNA@reductions$wnn.umap@cell.embeddings[colnames(multiome_seurat_noNA), ] == multiome_seurat@reductions$wnn.umap@cell.embeddings[colnames(multiome_seurat_noNA), ])

#need to do pca and umap and stuffs
DefaultAssay(multiome_seurat_noNA) <- "RNA"

multiome_seurat_noNA <- SCTransform(multiome_seurat_noNA, verbose = FALSE, vst.flavor = "v2") %>%
  RunPCA() %>% RunUMAP(dims = 1:50,
                       reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

### ATAC pre-processing
DefaultAssay(multiome_seurat_noNA) <- "ATAC"
multiome_seurat_noNA <- RunTFIDF(multiome_seurat_noNA)
multiome_seurat_noNA <- FindTopFeatures(multiome_seurat_noNA, min.cutoff = 'q0')
multiome_seurat_noNA <- RunSVD(multiome_seurat_noNA)
multiome_seurat_noNA <- RunUMAP(multiome_seurat_noNA, reduction = 'lsi', dims = 2:50,
                        reduction.name = "umap.atac", reduction.key = "atacUMAP_")
### Estimate gene activities using accessibility of promoters
fragment_dir <- "/home/lel2/luan/Juno/multiome_Rusland/multiome/data/raw"
fragment_list <- c("atac_fragments_C1.tsv.gz", "atac_fragments_C2.tsv.gz", "atac_fragments_D.tsv.gz", "atac_fragments_T.tsv.gz")

#remember to generate fragment index file
#https://stackoverflow.com/questions/72524686/fragment-file-is-not-indexed-in-signac-r
#https://github.com/stuart-lab/signac/discussions/1504
#tabix -p bed fragments.tsv.gz

for (i in seq_along(fragment_list)){
  multiome_seurat_noNA@assays[["ATAC"]]@fragments[[i]]@path <- paste0(fragment_dir, "/", fragment_list[i])

}
#a customed function without index file
#https://www.biostars.org/p/9525109/

gene.activities <- GeneActivity(multiome_seurat_noNA)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
multiome_seurat_noNA[['Promoter']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(multiome_seurat_noNA) <- 'Promoter'
multiome_seurat_noNA <- SCTransform(multiome_seurat_noNA,
                            assay = 'Promoter', new.assay.name = 'SCT_prom',
                            verbose = FALSE, vst.flavor = "v2")
multiome_seurat_noNA <- RunPCA(multiome_seurat_noNA ,reduction.name = 'pca.prom', assay = "SCT_prom")
multiome_seurat_noNA <- RunUMAP(multiome_seurat_noNA, dims = 1:50, reduction = 'pca.prom',
                        reduction.name = 'umap.promoter', reduction.key = 'promUMAP_')


### calculate a joint ATAC-RNA WNN graph
## since we have matched cells --> use FindMultiModalNeighbors
#https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis

multiome_seurat_noNA <- FindMultiModalNeighbors(multiome_seurat_noNA, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
multiome_seurat_noNA <- RunUMAP(multiome_seurat_noNA, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
multiome_seurat_noNA <- FindClusters(multiome_seurat_noNA, graph.name = "wsnn",
                             algorithm = 3, resolution = 1, verbose = FALSE)

#get the SCTModel similar to multiome_seurat
multiome_seurat_noNA@assays[["SCT"]]@SCTModel.list <- multiome_seurat@assays[["SCT"]]@SCTModel.list
#remove NA cells

for (i in as.character(c(1:4))){
  cell_id <- rownames(multiome_seurat_noNA@assays[["SCT"]]@SCTModel.list[[i]]@cell.attributes[, ])

  cell_id <- cell_id[is.na(multiome_seurat_noNA$projection[cell_id]) == FALSE]

  multiome_seurat_noNA@assays[["SCT"]]@SCTModel.list[[i]]@cell.attributes <- multiome_seurat_noNA@assays[["SCT"]]@SCTModel.list[[i]]@cell.attributes[cell_id, ]

}

#get celltype colors
multiome_seurat_noNA@misc[["celltype_colors"]] <- old.scrna$clust.final.cols


DefaultAssay(multiome_seurat_noNA) <- "SCT"
DimPlot(multiome_seurat_noNA,
        reduction = "wnn.umap",
        group.by = "projection",
        label = TRUE,
        cols = multiome_seurat_noNA@misc[["celltype_colors"]])

DimPlot(multiome_seurat_noNA,
        reduction = "umap.rna",
        group.by = "projection",
        label = TRUE,
        cols = multiome_seurat_noNA@misc[["celltype_colors"]])

DimPlot(multiome_seurat_noNA,
        reduction = "umap.atac",
        group.by = "projection",
        label = TRUE,
        cols = multiome_seurat_noNA@misc[["celltype_colors"]])


#need to correct for multiple SCT models
#especially for find markers
#https://satijalab.org/seurat/reference/prepsctfindmarkers#examples

multiome_seurat_noNA <- PrepSCTFindMarkers(object = multiome_seurat_noNA)


#obtain metacell for SCT ----
##follow https://gfellerlab.github.io/MetacellAnalysisTutorial/

gamma = 50 # the requested graining level.
k_knn = 30 # the number of neighbors considered to build the knn network.
#nb_var_genes = 2000 # number of the top variable genes to use for dimensionality reduction
nb_pc = 50 # the number of principal components to use.

#to find membership: which cells belong to which metacells
DefaultAssay(multiome_seurat_noNA) <- "SCT"
MC <- SuperCell::SCimplify(Seurat::GetAssayData(multiome_seurat_noNA, assay = "SCT", slot = "data"),  # single-cell log-normalized gene expression data
                           k.knn = k_knn,
                           gamma = gamma,
                           # n.var.genes = nb_var_genes,
                           n.pc = nb_pc,
                           genes.use = Seurat::VariableFeatures(multiome_seurat_noNA, nfeatures = 3000)
)

MC$annotation <- supercell_assign(clusters = multiome_seurat_noNA$projection, # single-cell annotation
                                  supercell_membership = MC$membership, # single-cell assignment to metacells
                                  method = "jaccard"
)


#construct metacells
MC.GE <- supercell_GE(Seurat::GetAssayData(multiome_seurat_noNA, assay = "SCT", slot = "counts"),
                      MC$membership,
                      mode =  "sum"
)
dim(MC.GE)
#do this so that cell names in MC.GE is the same as cell names in MC$annotation
colnames(MC.GE) <- paste0(rep(1:dim(MC.GE)[2]))

#get ATAC counts
MC_ATAC.GE <- supercell_GE(Seurat::GetAssayData(multiome_seurat_noNA, assay = "ATAC", slot = "counts"),
                           MC$membership,
                           mode =  "sum"
)
colnames(MC_ATAC.GE) <- colnames(MC.GE)
dim(MC_ATAC.GE)

#https://stuartlab.org/signac/articles/data_structures
chromatinassay <- CreateChromatinAssay(counts = MC_ATAC.GE, genome = "hg38")
df_size <- data.frame(size = as.vector(table(MC$membership)))
rownames(df_size) <- colnames(MC.GE)

MC_seurat <- CreateSeuratObject(counts = chromatinassay, assay = "ATAC", meta.data = df_size)
MC_seurat[["ATAC"]]@annotation <- multiome_seurat@assays[["ATAC"]]@annotation

#get motif information
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
# add motif information
# for similar commands, see 3_check_marker_genes_peaks.R
MC_seurat <- AddMotifs(MC_seurat, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm_set)


MC_seurat[["RNA"]] <- CreateAssayObject(counts = MC.GE, assay = "RNA")

MC_seurat@meta.data[["cluster_cell"]] <- MC$annotation
#get color for each metacells
MC_seurat@misc[["celltype_colors"]] <- old.scrna$clust.final.cols
#membership of each cell to metacells
# membership_df <- data.frame(row.names = names(MC$membership), membership = MC$membership, celltype = multiome_seurat_noNA@meta.data[["projection"]])
membership_df <- data.frame(row.names = names(MC$membership), membership = MC$membership)
membership_df$cluster_cell <- MC$annotation[membership_df$membership]
head(membership_df)
membership_df <- membership_df %>% tibble::rownames_to_column("cell")
rownames(membership_df) <- membership_df$cell
MC_seurat@misc[["membership_df"]] <- membership_df

#metacell normalization----

### RNA pre-processing

DefaultAssay(MC_seurat) <- "RNA"
library(dplyr)

MC_seurat <- SCTransform(MC_seurat, verbose = FALSE, vst.flavor = "v2") %>%
  RunPCA() %>% RunUMAP(dims = 1:50,
                       reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

#
# MC_seurat <- NormalizeData(MC_seurat)
# MC_seurat <- FindVariableFeatures(MC_seurat, nfeatures = 2000)
# MC_seurat <- ScaleData(MC_seurat)
# #> Centering and scaling data matrix
# MC_seurat <- RunPCA(MC_seurat, npcs = 50)
# MC_seurat <- RunUMAP(MC_seurat, reduction = "pca", dims = 1:50, verbose = F)
# MC_seurat <- FindNeighbors(MC_seurat, reduction = "pca", dims = 1:50)
# MC_seurat <- FindClusters(MC_seurat, resolution = 0.5)

### ATAC pre-processing
DefaultAssay(MC_seurat) <- "ATAC"
MC_seurat <- RunTFIDF(MC_seurat)
MC_seurat <- FindTopFeatures(MC_seurat, min.cutoff = 'q0')
MC_seurat <- RunSVD(MC_seurat)
MC_seurat <- RunUMAP(MC_seurat, reduction = 'lsi', dims = 2:50,
                     reduction.name = "umap.atac", reduction.key = "atacUMAP_")

MC_seurat@active.ident <- as.factor(MC_seurat$cluster_cell)

### calculate a joint ATAC-RNA WNN graph
## since we have matched cells --> use FindMultiModalNeighbors
#https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis
DefaultAssay(MC_seurat) <- "SCT"
MC_seurat <- FindMultiModalNeighbors(MC_seurat, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
MC_seurat <- RunUMAP(MC_seurat, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
MC_seurat <- FindClusters(MC_seurat, graph.name = "wsnn",
                                     algorithm = 3, resolution = 1, verbose = FALSE)

DimPlot(MC_seurat,
        reduction = "wnn.umap",
        label = TRUE,
        group.by = "cluster_cell")

DimPlot(MC_seurat,
        reduction = "umap.rna",
        label = TRUE,
        group.by = "cluster_cell")

DimPlot(MC_seurat,
        reduction = "umap.atac",
        label = TRUE,
        group.by = "cluster_cell")

#if needed to save memory
file_name <- ls()
file_name_remove <- setdiff(file_name, c("MC_seurat", "multiome_seurat_noNA", "MC", "old.scrna"))
rm(list = file_name_remove)
gc()

save.image(file='/data1/soldatr/luan/projects/motifregression/result/CH/multiome_metacell.RData')
