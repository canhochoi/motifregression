library(Seurat)
library(Signac)
library(SuperCell)
library(scattermore)
library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library(parallel)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

library(monocle3)
library(SeuratWrappers)


# has multiome_seurat_noNA
load(file='/data1/soldatr/luan/projects/motifregression/result/CH/multiome_metacell.RData')

multiome_seurat <- readRDS("~/luan/Juno/multiome_Rusland/multiome/data/processed/multiome_seurat.rds")


# try slingshot
# https://nbisweden.github.io/workshop-archive/workshop-scRNAseq/2020-01-27/labs/compiled/slingshot/slingshot.html
# https://nbisweden.github.io/workshop-scRNAseq/labs/seurat/seurat_07_trajectory.html#preparing-data

library(slingshot)

# Save the objects as separate matrices for input in slingshot
DefaultAssay(multiome_seurat) <- "SCT"
dimred <- multiome_seurat@reductions$wnn.umap@cell.embeddings
#select only ct 1, 6, 7, 8
dimred <- dimred[Idents(multiome_seurat) %in% c(1, 6, 7, 8), ]

clustering <- multiome_seurat$projection

clustering <- clustering[Idents(multiome_seurat) %in% c(1, 6, 7, 8)]
#remove NAs cell
# dimred <- dimred[!is.na(clustering), ]
#remove NAs
# clustering <- clustering[!is.na(clustering)]

counts <- as.matrix(multiome_seurat[["SCT"]]$counts[VariableFeatures(multiome_seurat),])


set.seed(1)
lineages <- getLineages(data = dimred,
                        clusterLabels = clustering,
                        end.clus = "8",
                        start.clus = "1")

lineages <- as.SlingshotDataSet(lineages)

lineages

plot(dimred, col = pal[clustering],  pch = 16)
lines(lineages, lwd = 3, col = 'black')


#define principal curves

curves <- getCurves(lineages, approx_points = NULL, thresh = 0.001, stretch = 0, allow.breaks = FALSE, shrink = 0.99)
curves <- as.SlingshotDataSet(curves)

pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))

plot(dimred, col = pal[clustering], asp = 1, pch = 16)
lines(curves, lwd = 3, col = "black")

pseudotime <- slingPseudotime(curves, na = FALSE)
cellWeights <- slingCurveWeights(curves)

plot(pseudotime,
     multiome_seurat[["SCT"]]$data["GATA2", rownames(pseudotime)],
     col = pal[clustering])

#find a way to smooth

#try single cell data
Idents(multiome_seurat_noNA) <- multiome_seurat_noNA$projection
multiome_seurat_noNA_subset <- subset(multiome_seurat_noNA, idents = c(1,6,7,8))
DimPlot(multiome_seurat_noNA_subset, reduction = "wnn.umap")
Reductions(multiome_seurat_noNA_subset)

cds <- as.cell_data_set(multiome_seurat_noNA_subset)
reducedDims(cds)
#change WNN.UMAP to UMAP
names(cds@int_colData@listData[["reducedDims"]])[5] <- "UMAP"

clusters <- multiome_seurat_noNA$projection[colnames(cds)]

partitions <- factor(rep(1, length(clusters)))
names(partitions) <- colnames(cds)

cds@clusters$UMAP <- list(
  clusters = factor(clusters),
  partitions = partitions,
  louvain_res = NA,
  k = NA,
  cluster_info = data.frame(
    cluster = levels(factor(clusters)),
    size = as.numeric(table(clusters))
  ),
  full_observation_matrix = NULL
)

#plotting
plot_cells(cds, show_trajectory_graph = FALSE,
           color_cells_by = "cluster",
           group_cells_by = "cluster",
           cell_size = 1)

plot_cells(cds, show_trajectory_graph = FALSE,
           color_cells_by = "partition",
           cell_size = 1)

#build principle graph
#https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/

cds <- learn_graph(cds, use_partition = TRUE)

plot_cells(cds, label_branch_points = TRUE, label_leaves = TRUE)

#https://stuartlab.org/signac/articles/monocle
#https://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html

# plot the pseudotime graph:
p1 <- plot_cells(
  cds = cds,
  color_cells_by = "cluster",
  group_cells_by = "cluster",
  show_trajectory_graph = TRUE,
  label_principal_points = TRUE,
  cell_size = 1
)

# plot the UMAP partitions from the clustering algorithm
p2 <- plot_cells(
  cds = cds,
  color_cells_by = "partition",
  show_trajectory_graph = FALSE
)

p1 + p2

#pick origin

cell_ids <- which(colData(cds)[, "projection"] == 1)

closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names                                                            (which.max(table(closest_vertex[cell_ids,]))))]



cds <- order_cells(cds, root_pr_nodes= root_pr_nodes)
# cds <- order_cells(cds, root_cells = rownames(colData(cds))[cell_ids])

p1 <- plot_cells(
  cds = cds,
  color_cells_by = "ident",
  group_cells_by = "ident",
  show_trajectory_graph = FALSE,
  label_principal_points = FALSE,
  label_branch_points = FALSE,
  label_leaves = FALSE,
  label_cell_groups = FALSE,
  cell_size = 1
)

p2 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_cell_groups=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 graph_label_size=1.5,
                 cell_size = 1)
p1 + p2


#use metacell ----
#https://stuartlab.org/signac/articles/monocle
Idents(MC_seurat) <- MC_seurat$cluster_cell
MC_seurat_subset <- subset(MC_seurat, idents = c(1,6,7,8))
Reductions(MC_seurat_subset)
cds <- as.cell_data_set(MC_seurat_subset)
reducedDims(cds)
#change WNN.UMAP to UMAP
names(cds@int_colData@listData[["reducedDims"]])[4] <- "UMAP"

clusters <- MC_seurat_subset$cluster_cell[colnames(cds)]

partitions <- factor(rep(1, length(clusters)))
names(partitions) <- colnames(cds)

cds@clusters$UMAP <- list(
  clusters = factor(clusters),
  partitions = partitions,
  louvain_res = NA,
  k = NA,
  cluster_info = data.frame(
    cluster = levels(factor(clusters)),
    size = as.numeric(table(clusters))
  ),
  full_observation_matrix = NULL
)

#plotting
plot_cells(cds, show_trajectory_graph = FALSE,
           color_cells_by = "cluster_cell",
           group_cells_by = "cluster_cell",
           cell_size = 1)

plot_cells(cds, show_trajectory_graph = FALSE,
           color_cells_by = "partition",
           cell_size = 1)

#build principle graph
#https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/

cds <- learn_graph(cds, use_partition = TRUE)

#https://stuartlab.org/signac/articles/monocle
#https://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html

# plot the pseudotime graph:
p1 <- plot_cells(
  cds = cds,
  color_cells_by = "cluster_cell",
  group_cells_by = "cluster_cell",
  show_trajectory_graph = TRUE,
  label_principal_points = TRUE,
  cell_size = 1
)

# plot the UMAP partitions from the clustering algorithm
p2 <- plot_cells(
  cds = cds,
  color_cells_by = "partition",
  show_trajectory_graph = FALSE
)

p1 + p2

#pick origin

cell_ids <- which(colData(cds)[, "cluster_cell"] == 1)

closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                            (which.max(table(closest_vertex[cell_ids,]))))]
cds <- order_cells(cds, root_pr_nodes= root_pr_nodes)

p1 <- plot_cells(
  cds = cds,
  color_cells_by = "cluster_cell",
  group_cells_by = "cluster_cell",
  show_trajectory_graph = TRUE,
  label_principal_points = FALSE,
  cell_size = 1
)

p2 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_cell_groups=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 graph_label_size=1.5,
                 cell_size = 1)
p1 + p2


#option 2: try to do PCA again ----
#https://cole-trapnell-lab.github.io/monocle3/docs/clustering/
cds <- preprocess_cds(cds, num_dim = 50)

plot_pc_variance_explained(cds)

#use UMAP by default
cds <- reduce_dimension(cds)

cds <- cluster_cells(cds = cds, reduction_method = "UMAP")


