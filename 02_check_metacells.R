
#remotes::install_github("GfellerLab/MetacellAnalysisToolkit",upgrade = "never")

library(MetacellAnalysisToolkit)
library(dplyr)

# has multiome_seurat_noNA
load(file='/data1/soldatr/luan/projects/motifregression/result/CH/multiome_metacell.RData')

multiome_seurat <- readRDS("~/luan/Juno/multiome_Rusland/multiome/data/processed/multiome_seurat.rds")


#projection onto single cell

mc_p <- mc_projection(
  sc.obj = multiome_seurat_noNA,
  mc.obj = MC_seurat,
  cell.membership = MC_seurat@misc$membership_df,
  sc.reduction = "wnn.umap",
  sc.label = "projection",
  metacell.label = "cluster_cell"
)


# Extract UMAP coordinates from your Seurat object
umap_df <- Embeddings(multiome_seurat_noNA, "wnn.umap") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  left_join(MC_seurat@misc$membership_df, by = c("cell" = "cell"))

head(umap_df)

# Compute average UMAP positions per cluster_cell
label_positions <- umap_df %>%
  group_by(cluster_cell) %>%
  summarize(UMAP_1 = mean(wnnUMAP_1), UMAP_2 = mean(wnnUMAP_2))
# find dominant cell type for each metacell
# summarise(UMAP_1 = mean(wnnUMAP_1), UMAP_2 = mean(wnnUMAP_2), celltype = names(sort(table(celltype), decreasing = TRUE))[1])

head(label_positions)

#check with cluster_cell membership
# table(label_positions$celltype == MC_seurat$cluster_cell)

# label_positions$cluster_cell <- MC_seurat$cluster_cell

mc_p +
  geom_text(data = label_positions,
            aes(x = UMAP_1, y = UMAP_2, label = cluster_cell),
            size = 4, fontface = "bold", inherit.aes = FALSE)

##purity
MC_seurat$purity <- mc_purity(membership = MC_seurat@misc$membership_df$membership, annotation = multiome_seurat_noNA$projection)

qc_boxplot(mc.obj = MC_seurat, qc.metrics = "purity")
qc_boxplot(mc.obj = MC_seurat, qc.metrics = "purity", split.by = "cluster_cell")

#compactness
MC_seurat$compactness <- mc_compactness(cell.membership = MC_seurat@misc$membership_df,
                                        sc.obj = multiome_seurat_noNA,
                                        sc.reduction = "pca",
                                        dims = 1:30)
qc_boxplot(mc.obj = MC_seurat, qc.metrics = "compactness")
qc_boxplot(mc.obj = MC_seurat, qc.metrics = "compactness", split.by = "cluster_cell")

#separation
MC_seurat$separation <- mc_separation(cell.membership = MC_seurat@misc$membership_df,
                                      sc.obj = multiome_seurat_noNA,
                                      sc.reduction = "pca")
qc_boxplot(mc.obj = MC_seurat, qc.metrics = "separation")
qc_boxplot(mc.obj = MC_seurat, qc.metrics = "separation", split.by = "cluster_cell")

#INV

mc_INV_modified <- function (sc.obj, cell.membership, assay = "RNA", group.label = "membership")
{
  if (identical(multiome_seurat[["SCT"]]$counts@x, multiome_seurat[["SCT"]]$data@x)) {
    message("Counts and data slots are identical.")
    message("Normalizing data ...")
    sc.obj <- Seurat::NormalizeData(sc.obj, normalization.method = "LogNormalize")
  }
  message("Computing INV ...")
  memberships_without_outliers <- na.exclude(cell.membership)
  membership_vector <- memberships_without_outliers[, group.label]
  names(membership_vector) <- rownames(memberships_without_outliers)
  sc.obj <- sc.obj[, names(membership_vector)]
  INV_val <- stats::aggregate(Matrix::t(Seurat::GetAssayData(sc.obj,
                                                             slot = "data")), by = list(metacell = membership_vector),
                              FUN = function(x) (1/mean(x)) * var(x))
  rownames(INV_val) <- INV_val[, 1]
  INV_val <- INV_val[, -1]
  INV_val[is.na(INV_val)] <- 0
  INV_val_qt <- apply(INV_val, 1, function(x) quantile(x,
                                                       0.95, na.rm = TRUE))
  return(INV_val_qt)
}

MC_seurat$INV <- mc_INV_modified(cell.membership = MC_seurat@misc$membership_df, sc.obj = multiome_seurat_noNA, group.label = "membership")
qc_boxplot(mc.obj = MC_seurat, qc.metrics = "INV")

