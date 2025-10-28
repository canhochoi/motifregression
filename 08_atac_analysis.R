devtools::load_all("~/luan/projects/motifregression/")
devtools::document()

library(Seurat)
library(Signac)
library(SuperCell)

#for motif information
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(ggplot2)

#for annotating peaks
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

library(ComplexHeatmap)



load(file='/data1/soldatr/luan/projects/motifregression/result/CH/multiome_metacell.RData')

DefaultAssay(multiome_seurat) <- "ATAC"

# add motifs information ----
# get list of motif PWM matrix for human
pwm <- getMatrixSet(x = JASPAR2020,
                    opts = list(species = 9606, all_versions = FALSE)
                    )

#add motif object to assay
#https://stuartlab.org/signac/articles/motif_vignette
# DefaultAssay(multiome_seurat) <- "ATAC"
# multiome_seurat <- AddMotifs(object = multiome_seurat,
#                        genome = BSgenome.Hsapiens.UCSC.hg38,
#                        pfm = pwm)

# multiome_seurat <- RunChromVAR(multiome_seurat,
#                                genome = BSgenome.Hsapiens.UCSC.hg38)
#check CEBPA
# id_cell <- which(Idents(multiome_seurat) %in% c(1, 6, 7,8))

# mean(multiome_seurat[["chromvar"]]$data["MA0102.4", which(Idents(multiome_seurat) == 8)])

# max(multiome_seurat[["chromvar"]]$data["MA0102.4", which(Idents(multiome_seurat) == 8)])
# find differentially accessibile peaks for heatmap ----
# use 10 cell types
Idents(multiome_seurat) <- factor(multiome_seurat$projection, level = sort(as.numeric(levels(multiome_seurat))))
Idents(multiome_seurat)

#https://stuartlab.org/signac/articles/pbmc_vignette.html#:~:text=Find%20differentially%20accessible%20peaks%20between%20cell%20types&text=A%20simple%20approach%20is%20to,run%20on%20a%20Seurat%20object.

de_peaks_ct <- FindAllMarkers(multiome_seurat,
                              assay = "ATAC",
                              only.pos = FALSE,
                              min.pct = 0.1)

head(de_peaks_ct)

DA_threshold <- 1

de_peaks_sg <- de_peaks_ct %>%
                filter(p_val_adj <= 0.05) %>%
                mutate(regulation = case_when(
                  avg_log2FC > DA_threshold ~ "up",
                  avg_log2FC < -DA_threshold ~ "down",
                  TRUE ~ "neutral"
                ))

de_peaks_array <- unique(de_peaks_sg$gene)

celltype_array <- c(1, 6, 7, 8)

mean_de_peaks_ct <- lapply(celltype_array, function(ct){
  ct_id <- which(Idents(multiome_seurat) == ct)
  mean_peaks <- apply(multiome_seurat[["ATAC"]]$data[de_peaks_array, ct_id], 1, mean)
  return (mean_peaks)
  })

mean_de_peaks_ct <- do.call(cbind, mean_de_peaks_ct)
colnames(mean_de_peaks_ct) <- celltype_array


mean_de_peaks_ct_zscaled <- t(scale(t(mean_de_peaks_ct)))

all(rownames(mean_de_peaks_ct_zscaled) == de_peaks_array)

library(circlize)

library(RColorBrewer)

# Step 1: Assign colors
colors <- setNames(brewer.pal(n = length(celltype_array), name = "Set2"), celltype_array)


# Column annotation (just labels in this case)
col_annot <- HeatmapAnnotation(
  CellType = as.character(celltype_array),
  annotation_name_side = "left",
  col = list(CellType = colors)
)

peak_regulation <- de_peaks_sg[match(rownames(mean_de_peaks_ct_zscaled), de_peaks_sg$gene), "regulation"]

# Row annotation (up/down/neutral as colored bar)
row_annot <- rowAnnotation(
  Regulation = peak_regulation,
  col = list(Regulation = c("up" = "firebrick", "down" = "navy", "neutral" = "gray")),
  annotation_legend_param = list(title = "Peak Regulation")
)


Heatmap(mean_de_peaks_ct_zscaled,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        top_annotation = col_annot)

# find differentially accessible peaks ----
# use 10 cell types
Idents(multiome_seurat) <- factor(multiome_seurat$projection, level = sort(as.numeric(levels(multiome_seurat))))
Idents(multiome_seurat)

#https://stuartlab.org/signac/articles/pbmc_vignette.html#:~:text=Find%20differentially%20accessible%20peaks%20between%20cell%20types&text=A%20simple%20approach%20is%20to,run%20on%20a%20Seurat%20object.

de_peaks_ct <- FindAllMarkers(multiome_seurat,
                              assay = "ATAC",
                              only.pos = TRUE,
                              min.pct = 0.1)
de_peaks_ct <- de_peaks_ct %>%
                  filter(pct.1 > pct.2 & p_val_adj <= 0.05)

#plot distribution of log2FC across cell type
library(ggplot2)
ggplot(de_peaks_ct, aes(x = cluster, y = avg_log2FC, fill = cluster)) +
  geom_violin() +
  geom_jitter(position = position_jitter(seed = 1, width = 0.2), alpha = 0.1)

# summarize number of differentially accessible (DA) peaks per cell type
peaks_ct_summary <- de_peaks_ct %>% count(cluster)

peaks_ct_summary
# number of unique DA peaks

length(unique(de_peaks_ct$gene))


# change peaks name to GRange objects ----
de_peaks_gr <- de_peaks_ct %>%
  separate(gene, into = c("chr", "start", "end"), sep = "-", convert = TRUE) %>%
  drop_na(chr, start, end) %>%  # optional, clean malformed rows
  makeGRangesFromDataFrame(seqnames.field = "chr",
                           start.field = "start",
                           end.field = "end",
                           keep.extra.columns = TRUE)


# Annotate peaks
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peak_anno <- annotatePeak(de_peaks_gr, TxDb = txdb, annoDb = "org.Hs.eg.db")

# View the annotation categories
head(as.data.frame(peak_anno))$annotation

peak_anno_df <- as.data.frame(peak_anno)
peak_anno_df$peak_name <- paste0(peak_anno_df$seqnames, "-", peak_anno_df$start, "-", peak_anno_df$end)

#check difference
setdiff(de_peaks_ct$gene, peak_anno_df$peak_name)

head(peak_anno_df)

# check how many peaks shared between HSC and cell type 8 ----

# seems not to be useful
peaks_ct_1 <- de_peaks_ct %>%
  filter(cluster == 1) %>%
  pull("gene")

i <- 8
peaks_ct_others <- de_peaks_ct %>%
  filter(cluster == i) %>%
  pull("gene")

peaks_shared <- intersect(peaks_ct_1, peaks_ct_others)

de_peaks_bw_ct <- FindMarkers(multiome_seurat,
                              assay = "ATAC",
                              features = peaks_shared,
                              ident.1 = 1,
                              ident.2 = i,
                              min.pct = 0.1)

# try to compare directly via differential test ----
# cell type 8
i <- 8
de_peaks_bw_ct <- FindMarkers(multiome_seurat,
                              assay = "ATAC",
                              ident.1 = 1,
                              ident.2 = i,
                              min.pct = 0.1)
de_peaks_bw_ct[, c("pct.1", "pct.2")]


# get shared peaks with no significant change statistically
peak_shared <- de_peaks_bw_ct %>%
  filter(p_val_adj <= 0.05 & abs(avg_log2FC) <= 1) %>%
  rownames()


# calculate peak mean for each cell type
peaks_mean_ls <- lapply(c(1, 6, 7, 8), function(ct){
  cell_id <- which(Idents(multiome_seurat) == ct)
  #take mean of each peak across a cell type
  peaks_mean <- apply(multiome_seurat[["ATAC"]]$data[peaks_shared, cell_id], 1, mean)
  return (peaks_mean)
})

peaks_mean_ls <- do.call(cbind, peaks_mean_ls)
colnames(peaks_mean_ls) <- c(1, 6, 7, 8)

anno_col <- data.frame(celltype = c(1, 6, 7, 8))
anno_col$celltype <- as.factor(anno_col$celltype)
#important to get correct annotation
rownames(anno_col) <- colnames(peaks_mean_ls)


library(RColorBrewer)

# Step 1: Get unique groups
celltypes <- levels(anno_col$celltype)

# Step 2: Assign colors
colors <- setNames(brewer.pal(n = length(celltypes), name = "Set2"), celltypes)

# Step 3: Create annotation color list
anno_colors <- list(celltype = colors)


library(pheatmap)

pheatmap(peaks_mean_ls,
         scale = "none",
         annotation_col = anno_col,
         annotation_colors = anno_colors,
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_names_col = FALSE,
         annotation_names_row = FALSE)

library(ComplexHeatmap)
library(circlize)

# Create top annotation
col_anno <- HeatmapAnnotation(
  CellType = anno_col$celltype,
  col = list(CellType = colors),
  show_annotation_name = FALSE
)

# Plot heatmap with annotation
Heatmap(peaks_mean_ls,
        name = "Peak accessibility",  # This sets the legend title for the color scale
        top_annotation = col_anno,    # Add column annotation
        # col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        col = colorRampPalette(c("blue", "white", "red"))(100),
        show_row_names = FALSE,
        show_column_names = FALSE,
        cluster_columns = FALSE)

peaks_mean_ls_zscaled <- t(scale(t(peaks_mean_ls)))

Heatmap(peaks_mean_ls_zscaled,
        name = "z score of peak accessibility",  # This sets the legend title for the color scale
        top_annotation = col_anno,    # Add column annotation
        # col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        col = colorRampPalette(c("blue", "white", "red"))(100),
        show_row_names = FALSE,
        show_column_names = FALSE,
        cluster_columns = FALSE)

#find enriched motifs

enriched_motif_celltype <- function(multiome_seurat, celltype, peak_shared){
  multiome_seurat_subset <- subset(multiome_seurat, idents = celltype)

  open.peaks <- AccessiblePeaks(multiome_seurat_subset, idents = celltype)

  # match the overall GC content in the peak set
  # use all cells to search for background
  meta.feature <- GetAssayData(multiome_seurat, assay = "ATAC", layer = "meta.features")
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[peak_shared, ],
    n = 50000
  )

  enriched_motif <- FindMotifs(multiome_seurat_subset,
                               features = peak_shared,
                               background = peaks.matched)

  return (enriched_motif)
}

celltype_array <- c(1, 6, 7, 8)
# celltype_array <- c(1, i)

enriched_motif_ls <- lapply(celltype_array, function(ct){
  enriched_motif <- enriched_motif_celltype(multiome_seurat, ct, peak_shared)
  enriched_motif$celltype <- rep(ct, nrow(enriched_motif))
  return (enriched_motif)
})

names(enriched_motif_ls) <- celltype_array
# enriched_motif_df <- do.call(cbind, enriched_motif_ls)

# enriched_motif_df <- data.frame(motif = c(enriched_motif_ls[[1]]$motif, enriched_motif_ls[[2]]$motif),
#                                 motif_name = c(enriched_motif_ls[[1]]$motif.name, enriched_motif_ls[[2]]$motif.name),
#                                 p_adjust = c(enriched_motif_ls[[1]]$p.adjust, enriched_motif_ls[[2]]$p.adjust),
#                                 enrichment = c(enriched_motif_ls[[1]]$fold.enrichment, enriched_motif_ls[[2]]$fold.enrichment),
#                                 celltype = c(enriched_motif_ls[[1]]$celltype, enriched_motif_ls[[2]]$celltype)
#                       )

# motif_array <- unlist(lapply(seq_along(celltype_array), function(ct_idx){enriched_motif_ls[[ct_idx]]$motif}))
# motif_array <- do.call(cbind, lapply(seq_along(celltype_array), function(ct_idx){enriched_motif_ls[[ct_idx]]$motif}))
enriched_motif_df <- data.frame(motif = unlist(lapply(seq_along(celltype_array), function(ct_idx){enriched_motif_ls[[ct_idx]]$motif})),
                                motif_name = unlist(lapply(seq_along(celltype_array), function(ct_idx){enriched_motif_ls[[ct_idx]]$motif.name})),
                                p_adjust = unlist(lapply(seq_along(celltype_array), function(ct_idx){enriched_motif_ls[[ct_idx]]$p.adjust})),
                                enrichment = unlist(lapply(seq_along(celltype_array), function(ct_idx){enriched_motif_ls[[ct_idx]]$fold.enrichment})),
                                celltype = unlist(lapply(seq_along(celltype_array), function(ct_idx){enriched_motif_ls[[ct_idx]]$celltype}))
                                )

enriched_motif_df$celltype <- as.factor(enriched_motif_df$celltype)

ggplot(enriched_motif_df, aes(x = enrichment, y = -log10(p_adjust), color = celltype)) +
  geom_point()

# Convert to wide format for enrichment comparison
enrichment_wide_df <- enriched_motif_df %>%
  dplyr::select(motif, motif_name, enrichment, celltype) %>%
  tidyr::pivot_wider(
    names_from = celltype,
    values_from = enrichment
  )

colnames(enrichment_wide_df)[3:ncol(enrichment_wide_df)] <- paste0("cell_", celltype_array)

enrichment_wide_df <- as.data.frame(enrichment_wide_df)

# plot motif enrichment score between cell types
# always show linear trend
ggplot(enrichment_wide_df, aes(x = cell_1, y = cell_6)) +
  geom_point()

pairs(enrichment_wide_df[, 3:6])
# motif_sig <- enriched_motif_df %>%
#   filter(p_adjust <= 0.05 & enrichment >= 1) %>%
#   pull(motif)
#
# motif_sig <- unique(enriched_motif_df$motif)

# filter significant motif enrichment
enriched_motif_df <- enriched_motif_df %>%
  filter(p_adjust <= 0.05 & enrichment >= 1)

motif_rank_df <- enriched_motif_df %>%
  arrange(celltype, p_adjust) %>%
  group_by(celltype)

motif_rank_df <- enriched_motif_df %>%
  group_by(celltype) %>%
  slice_min(order_by = -enrichment, n = 20, with_ties = FALSE) %>%
  # slice_min(order_by = p_adjust, n = 20, with_ties = FALSE) %>%
  ungroup()


motif_rank_df <- motif_rank_df %>%
  mutate(label = paste0(signif(enrichment, 2), " (p=", signif(p_adjust, 2), ")"))

motif_label_wide <- motif_rank_df %>%
  dplyr::select(motif, motif_name, celltype, label) %>%
  tidyr::pivot_wider(names_from = celltype, values_from = label)


# found EGR1, EGR2, TFDP1, E2F6, SP1, SP3, KLF5, KLF6, ELK1, KLF2
# motifs relevant for HSC -> myeloid differentiation

important_myeloid_TFs <- c(
  "EGR1", "EGR2", "TFDP1", "E2F6", "SP1", "SP3", "KLF5", "KLF6", "KLF2", "ELK1",
  "RUNX1", "SPI1", "CEBPA", "GATA2", "GFI1", "IRF8", "STAT5A", "STAT5B",
  "MYB", "HOXA9", "NFIL3",
  "KLF4", "MAFB", "ZEB2", "ID2", "BCL11A"
)

# check how many motifs enriched in each cell type

enriched_motif_df %>%
  group_by(celltype) %>%
  summarize(num_motifs = n())

filtered_df <- enriched_motif_df %>%
                filter(motif_name %in% important_myeloid_TFs)

filtered_df %>%
  filter(celltype == 8) %>%
  pull(motif_name)

motif_enriched_ct_list <- lapply(celltype_array, function(ct) {
  motif_enriched_ct <- enriched_motif_df %>%
  filter(celltype == ct & p_adjust <= 0.05) %>%
  pull(motif_name)
})

names(motif_enriched_ct_list) <- celltype_array

# Choose and order specific cell types
celltype_order <- c("celltype_1", "celltype_6", "celltype_7", "celltype_8")

# Pivot to wide format
enrichment_matrix <- filtered_df %>%
  filter(celltype %in% celltype_array) %>%
  dplyr::select(motif_name, celltype, enrichment) %>%
  pivot_wider(names_from = celltype, values_from = enrichment) %>%
  tibble::column_to_rownames("motif_name") %>%
  as.matrix()

# make sure in order
enrichment_matrix <- enrichment_matrix[, intersect(celltype_array, colnames(enrichment_matrix))]

library(ComplexHeatmap)

Heatmap(
  enrichment_matrix,
  name = "Motif Enrichment",
  top_annotation = col_anno,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE
)


filtered_df %>% summarise(.by = celltype)

common_motif <- intersect(enriched_motif_df$motif[enriched_motif_df$celltype == 1],
                          enriched_motif_df$motif[enriched_motif_df$celltype == i])
#unique motif for cell type 1
motif_1 <- setdiff(enriched_motif_df$motif[enriched_motif_df$celltype == 1], common_motif)
#unique motif for cell type 8
motif_8 <- setdiff(enriched_motif_df$motif[enriched_motif_df$celltype == i], common_motif)

plot(enriched_motif_df$enrichment[enriched_motif_df$motif %in% common_motif & enriched_motif_df$celltype == 1],
     enriched_motif_df$enrichment[enriched_motif_df$motif %in% common_motif & enriched_motif_df$celltype == i])


# summarize(
#   motifs = paste0(motif.name, " (", signif(p.adjust, 2), ")") %>%
#             head(10) %>%           # top 10 motifs (optional)
#             paste(collapse = ", ")
#   )

#check all cells
enriched_motif <- enriched_motif_celltype(multiome_seurat, c(1, 6, 7, 8), peak_shared)

enriched_motif %>% filter(p.adjust <= 0.05)

# check differentially accessible peaks in stem cells and myeloid ----

peak_ct_1 <- de_peaks_bw_ct %>%
  filter(p_val_adj <= 0.05 & avg_log2FC >= 1.5 & pct.1 >= pct.2) %>%
  rownames()

peak_ct_8 <- de_peaks_bw_ct %>%
  filter(p_val_adj <= 0.05 & avg_log2FC <= -1.5 & pct.1 <= pct.2) %>%
  rownames()

peaks_de <- c(peak_ct_1, peak_ct_8)


# calculate peak mean for each cell type
peaks_de_mean_ls <- lapply(c(1, 6, 7, 8), function(ct){
  cell_id <- which(Idents(multiome_seurat) == ct)
  #take mean of each peak across a cell type
  peaks_mean <- apply(multiome_seurat[["ATAC"]]$data[peaks_de, cell_id], 1, mean)
  return (peaks_mean)
})

peaks_de_mean_ls <- do.call(cbind, peaks_de_mean_ls)
colnames(peaks_de_mean_ls) <- celltype_array

Heatmap(
  peaks_de_mean_ls,
  name = "Peak accessibility",
  top_annotation = col_anno,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE
)

peaks_de_mean_ls_zscaled <- t(scale(t(peaks_de_mean_ls)))

Heatmap(
  peaks_de_mean_ls_zscaled,
  name = "z score of Peak accessibility",
  top_annotation = col_anno,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE
)

# check enriched motifs

enriched_motif_de_ls <- lapply(celltype_array, function(ct){
  enriched_motif <- enriched_motif_celltype(multiome_seurat, ct, peaks_de)
  enriched_motif$celltype <- rep(ct, nrow(enriched_motif))
  return (enriched_motif)
})

names(enriched_motif_de_ls) <- celltype_array

enriched_motif_de_df <- data.frame(motif = unlist(lapply(seq_along(celltype_array), function(ct_idx){enriched_motif_de_ls[[ct_idx]]$motif})),
                                   motif_name = unlist(lapply(seq_along(celltype_array), function(ct_idx){enriched_motif_de_ls[[ct_idx]]$motif.name})),
                                   p_adjust = unlist(lapply(seq_along(celltype_array), function(ct_idx){enriched_motif_de_ls[[ct_idx]]$p.adjust})),
                                   enrichment = unlist(lapply(seq_along(celltype_array), function(ct_idx){enriched_motif_de_ls[[ct_idx]]$fold.enrichment})),
                                   celltype = unlist(lapply(seq_along(celltype_array), function(ct_idx){enriched_motif_de_ls[[ct_idx]]$celltype}))
)

enriched_motif_de_df$celltype <- as.factor(enriched_motif_de_df$celltype)

# Convert to wide format for enrichment comparison
enrichment_wide_de_df <- enriched_motif_de_df %>%
  dplyr::select(motif, motif_name, enrichment, celltype) %>%
  tidyr::pivot_wider(
    names_from = celltype,
    values_from = enrichment
  )

colnames(enrichment_wide_de_df)[3:ncol(enrichment_wide_de_df)] <- paste0("cell_", celltype_array)

enrichment_wide_de_df <- as.data.frame(enrichment_wide_de_df)

# filter significant motif enrichment
enriched_motif_de_df <- enriched_motif_de_df %>%
  filter(p_adjust <= 0.05 & enrichment >= 1)

#find number of enriched motifs in each cell type
enriched_motif_de_df %>%
  group_by(celltype) %>%
  summarise(num_motifs = n())


motif_enriched_ct_de_list <- lapply(celltype_array, function(ct) {
  motif_enriched_ct <- enriched_motif_de_df %>%
    filter(celltype == ct & p_adjust <= 0.05) %>%
    pull(motif_name)
})

names(motif_enriched_ct_de_list) <- celltype_array

shared_motifs <- Reduce(intersect, motif_enriched_ct_de_list)

# find important myeloid motif
filtered_de_df <- enriched_motif_de_df %>%
  filter(motif_name %in% important_myeloid_TFs)

# find shared motifs
filtered_de_df <- enriched_motif_de_df %>%
  filter(motif_name %in% shared_motifs)


# Pivot to wide format
enrichment_de_matrix <- filtered_de_df %>%
  filter(celltype %in% celltype_array) %>%
  dplyr::select(motif_name, celltype, enrichment) %>%
  pivot_wider(names_from = celltype, values_from = enrichment) %>%
  tibble::column_to_rownames("motif_name") %>%
  as.matrix()

# make sure in order
enrichment_de_matrix <- enrichment_de_matrix[, intersect(celltype_array, colnames(enrichment_matrix))]

#remove rows containing NAs
enrichment_de_matrix <- enrichment_de_matrix[apply(!is.na(enrichment_de_matrix), 1, all), ]

Heatmap(
  enrichment_de_matrix,
  name = "Motif Enrichment",
  top_annotation = col_anno,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE
)


# check TF activity by chromvar ----

motif_sig_names <- c(rownames(enrichment_matrix), rownames(enrichment_de_matrix))
motif_sig_names <- unique(motif_sig_names)

# match motif names
motif_name <- unlist(multiome_seurat@assays[["ATAC"]]@motifs@motif.names)
motif_sig <- names(motif_name)[match(motif_sig_names, motif_name)]
# check
all(unlist(multiome_seurat@assays[["ATAC"]]@motifs@motif.names[motif_sig]) == motif_sig_names)


chromvar_mean_ls <- lapply(celltype_array, function(ct){
  id_cell <- which(Idents(multiome_seurat) == ct)
  chromvar_mean <- apply(multiome_seurat[["chromvar"]]$data[motif_sig, id_cell], 1, mean)
  return (chromvar_mean)
})

chromvar_mean_ls <- do.call(cbind, chromvar_mean_ls)
colnames(chromvar_mean_ls) <- celltype_array
rownames(chromvar_mean_ls) <- motif_sig_names

pheatmap::pheatmap(chromvar_mean_ls,
                   scale = "row",
                   cluster_cols = FALSE,
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   annotation_names_col = FALSE,
                   annotation_col = anno_col)

Heatmap(chromvar_mean_ls,
        name = "TF activity",  # This sets the legend title for the color scale
        top_annotation = col_anno,    # Add column annotation
        # col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        col = colorRampPalette(c("blue", "white", "red"))(100),
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_columns = FALSE)

chromvar_mean_ls_zscaled <- t(scale(t(chromvar_mean_ls)))

Heatmap(chromvar_mean_ls_zscaled,
        name = "z score of TF activity",  # This sets the legend title for the color scale
        top_annotation = col_anno,    # Add column annotation
        # col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        col = colorRampPalette(c("blue", "white", "red"))(100),
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_columns = FALSE)

# check TF gene expression ----

TF_expression_mean_ls <- lapply(celltype_array, function(ct){
  id_cell <- which(Idents(multiome_seurat) == ct)
  TF_expression_mean <- apply(expm1(multiome_seurat[["SCT"]]$data[motif_sig_names, id_cell]), 1, mean)
  return (log1p(TF_expression_mean))
})

TF_expression_mean_ls <- do.call(cbind, TF_expression_mean_ls)
colnames(TF_expression_mean_ls) <- celltype_array
rownames(TF_expression_mean_ls) <- motif_sig_names

pheatmap::pheatmap(TF_expression_mean_ls,
                   scale = "row",
                   cluster_cols = FALSE,
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   annotation_names_col = FALSE,
                   annotation_col = anno_col)

Heatmap(TF_expression_mean_ls,
        cluster_rows = FALSE,
        name = "TF exprsn",  # This sets the legend title for the color scale
        top_annotation = col_anno,    # Add column annotation
        # col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        col = colorRampPalette(c("blue", "white", "red"))(100),
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_columns = FALSE)

TF_expression_mean_ls_zscaled <- t(scale(t(TF_expression_mean_ls)))

Heatmap(TF_expression_mean_ls_zscaled,
        cluster_rows = FALSE,
        name = "z score of TF exprsn",  # This sets the legend title for the color scale
        top_annotation = col_anno,    # Add column annotation
        # col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        col = colorRampPalette(c("blue", "white", "red"))(100),
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_columns = FALSE)

#
# TF_expression_ls <- lapply(celltype_array, function(ct){
#   id_cell <- which(Idents(multiome_seurat) == ct)
#   exp_mat <- multiome_seurat[["SCT"]]$data[motif_sig_names, id_cell]
#   exp_mat <- rbind(rep(ct, ncol(exp_mat)), exp_mat)
#   return (exp_mat)
# })
#
# TF_expression_ls <- do.call(cbind, TF_expression_ls)
#
# TF_expression_ls <- as.data.frame(t(TF_expression_ls))
# colnames(TF_expression_ls)[1] <- "celltype"
TF_expression_ls$celltype <- as.factor(TF_expression_ls$celltype)

TF_expression_ls_long_df <- TF_expression_ls %>%
                              pivot_longer(cols = c(2:15),
                                           names_to = "TF",
                                           values_to = "Expression")
ggplot(TF_expression_ls_long_df, aes(x = TF, y = Expression, fill = celltype)) +
  geom_violin(scale = "width", trim = TRUE) +
  theme_minimal() +
  labs(x = "Transcription Factor", y = "Value", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


multiome_seurat_subset <- subset(multiome_seurat, idents = c(1, 6, 7, 8))

dot_p <- DotPlot(multiome_seurat_subset,
        assay = "SCT",
        features = motif_sig_names,
        scale = FALSE)

dot_p +
  coord_flip() +
  scale_x_discrete(position = "bottom") +    # move celltype (now x) labels to top
  scale_y_discrete(position = "left") +  # move TF (now y) labels to right
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),   # top axis labels (celltypes)
    axis.text.y = element_text(hjust = 0),                # right axis labels (TFs)
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "right"
  )

dot_p_zscaled <- DotPlot(multiome_seurat_subset,
                 assay = "SCT",
                 features = motif_sig_names,
                 scale = TRUE)

dot_p_zscaled +
  coord_flip() +
  scale_x_discrete(position = "bottom") +    # move celltype (now x) labels to top
  scale_y_discrete(position = "left") +  # move TF (now y) labels to right
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),   # top axis labels (celltypes)
    axis.text.y = element_text(hjust = 0),                # right axis labels (TFs)
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "right"
  )

chromvar_dot_p <- DotPlot(multiome_seurat_subset,
                 assay = "chromvar",
                 features = motif_sig,
                 scale = FALSE)

chromvar_df <- as.data.frame(chromvar_dot_p[["data"]])

# need to do this because DotPlot use expm1 for averaging
mean_chromvar_ls <- lapply(celltype_array, function(ct){
  mean_chromvar <- apply(multiome_seurat[["chromvar"]]$data[motif_sig, which(Idents(multiome_seurat) == ct)], 1, mean)
  return (mean_chromvar)
})

# mean_chromvar_ls <- do.call(rbind, mean_chromvar_ls)
mean_chromvar_ls <- unlist(mean_chromvar_ls)

chromvar_dot_p[["data"]]$avg.exp <- as.numeric(mean_chromvar_ls)

chromvar_dot_p$data$features.plot <- factor(motif_sig_names[match(chromvar_dot_p$data$features, motif_sig)], levels = motif_sig_names)


# the plot use scaled average expression
# so previous is only for plotting the scatter plot later
chromvar_dot_p +
  coord_flip() +
  scale_x_discrete(position = "bottom") +    # move celltype (now x) labels to top
  scale_y_discrete(position = "left") +  # move TF (now y) labels to right
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),   # top axis labels (celltypes)
    axis.text.y = element_text(hjust = 0),                # right axis labels (TFs)
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "right"
  )

chromvar_dot_p_zscaled <- DotPlot(multiome_seurat_subset,
                         assay = "chromvar",
                         features = motif_sig,
                         scale = TRUE)

chromvar_dot_p_zscaled$data$features.plot <- factor(motif_sig_names[match(chromvar_dot_p_zscaled$data$features, motif_sig)], levels = motif_sig_names)

chromvar_dot_p_zscaled +
  coord_flip() +
  scale_x_discrete(position = "bottom") +    # move celltype (now x) labels to top
  scale_y_discrete(position = "left") +  # move TF (now y) labels to right
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),   # top axis labels (celltypes)
    axis.text.y = element_text(hjust = 0),                # right axis labels (TFs)
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "right"
  )

#plot mean TF gene expression vs activity for TFs
exprsn_activity_df <- data.frame(expression = dot_p[["data"]][["avg.exp"]],
                                 activity = chromvar_dot_p[["data"]][["avg.exp"]],
                                 motif =  dot_p[["data"]][["features.plot"]],
                                 celltype = dot_p[["data"]][["id"]])

exprsn_activity_df$celltype <- as.factor(exprsn_activity_df$celltype)


ggplot(exprsn_activity_df, aes(x = expression, y = activity, color = celltype)) +
  geom_point() +
  facet_wrap(~ celltype) +
  labs(x = "Average TF expression",
       y = "Average TF activity")

ggplot(exprsn_activity_df, aes(x = expression, y = activity, color = celltype)) +
  geom_point() +
  labs(x = "Average TF expression",
       y = "Average TF activity")


#check correlation

#check motif name matches motif IDs

id_cells <- which(Idents(multiome_seurat) %in% c(1, 6, 7, 8))

# TF_name <- "NFIL3"
# TF_name <- "SP3"

# check correlation
cor_array <- sapply(motif_sig_names, function(TF_name){
  TF_id <- motif_sig[match(TF_name, motif_sig_names)]
  cor_val <- cor(multiome_seurat[["SCT"]]$data[TF_name, id_cells],
      multiome_seurat[["chromvar"]]$data[TF_id, id_cells])
  return (cor_val)
})

cor_array
# multiome_seurat[["ATAC"]]@motifs@motif.names[[motif_sig[match("NFIL3", motif_sig_names)]]]

plot(multiome_seurat[["SCT"]]$data[TF_name, id_cells],
     multiome_seurat[["chromvar"]]$data[TF_id, id_cells])



