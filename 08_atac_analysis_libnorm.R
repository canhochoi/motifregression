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

# find shared peaks with similar accessibility in cell type 1 and 8 ----
# another way
# find shared peaks from motifregression package
peak_shared <- find_shared_peaks(multiome_seurat, 1, 8)

# remove chrM
peak_shared <- peak_shared[!grepl("^chrM-", peak_shared)]



multiome_seurat <- normalize_library_size(multiome_seurat, lib_size = 1e4)

# calculate peak mean for each cell type
peaks_mean_ls <- lapply(c(1, 6, 7, 8), function(ct){
  cell_id <- which(Idents(multiome_seurat) == ct)
  #take mean of each peak across a cell type
  peaks_mean <- apply(multiome_seurat[["ATAC_libnormed"]]$data[peak_shared, cell_id], 1, mean)
  # peaks_mean <- log1p(peaks_mean)
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

library(tidyverse)
peaks_mean_ls_df <- as.data.frame(peaks_mean_ls)
peaks_mean_ls_df$peaks <- rownames(peaks_mean_ls)

peaks_mean_ls_df <- peaks_mean_ls_df %>%
                      pivot_longer(
                        cols = -peaks,
                        names_to = "celltype",
                        values_to = "peak_accessibility"
                      )

ggplot(peaks_mean_ls_df, aes(x = celltype, y = peak_accessibility)) +
  geom_boxplot() +
  # geom_jitter() +
  # geom_box(trim = TRUE, scale = "width", fill = "skyblue") +
  theme_minimal() +
  labs(title = "Peak Accessibility per Cell Type",
       x = "Cell Type",
       y = "Mean Peak Accessibility")

peaks_mean_ls_zscaled <- t(scale(t(peaks_mean_ls)))

Heatmap(peaks_mean_ls_zscaled,
        name = "z score of peak accessibility",  # This sets the legend title for the color scale
        top_annotation = col_anno,    # Add column annotation
        # col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        col = colorRampPalette(c("blue", "white", "red"))(100),
        show_row_names = FALSE,
        show_column_names = FALSE,
        cluster_columns = FALSE)



# find enriched motifs ----



celltype_array <- c(1, 6, 7, 8)
# celltype_array <- c(1, i)

# enriched_motif_celltype is in motifregression package
# atac_utils.R
enriched_motif_ls <- lapply(celltype_array, function(ct){
  enriched_motif <- enriched_motif_celltype(multiome_seurat, ct, peak_shared)
  enriched_motif$celltype <- rep(ct, nrow(enriched_motif))
  return (enriched_motif)
})

names(enriched_motif_ls) <- celltype_array

enriched_motif_df <- data.frame(motif = unlist(lapply(seq_along(celltype_array), function(ct_idx){enriched_motif_ls[[ct_idx]]$motif})),
                                motif_name = unlist(lapply(seq_along(celltype_array), function(ct_idx){enriched_motif_ls[[ct_idx]]$motif.name})),
                                p_adjust = unlist(lapply(seq_along(celltype_array), function(ct_idx){enriched_motif_ls[[ct_idx]]$p.adjust})),
                                enrichment = unlist(lapply(seq_along(celltype_array), function(ct_idx){enriched_motif_ls[[ct_idx]]$fold.enrichment})),
                                celltype = unlist(lapply(seq_along(celltype_array), function(ct_idx){enriched_motif_ls[[ct_idx]]$celltype}))
)

enriched_motif_df$celltype <- as.factor(enriched_motif_df$celltype)

ggplot(enriched_motif_df, aes(x = enrichment, y = -log10(p_adjust), color = celltype)) +
  geom_point()


motif_rank_df <- enriched_motif_df %>%
  arrange(celltype, p_adjust) %>%
  group_by(celltype)

motif_rank_df <- enriched_motif_df %>%
  group_by(celltype) %>%
  slice_min(order_by = -enrichment, n = 20, with_ties = FALSE) %>%
  # slice_min(order_by = p_adjust, n = 20, with_ties = FALSE) %>%
  ungroup()


ggplot(motif_rank_df, aes(
  x = celltype,
  y = motif_name,
  size = enrichment,
  color = -log10(p_adjust)  # or use `color = p_adjust` for raw scale
)) +
  geom_point(alpha = 0.8) +
  scale_color_viridis_c(option = "C", direction = -1, name = "-log10(p-adjust)") +
  scale_size_continuous(name = "Enrichment") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.title = element_blank()
  ) +
  labs(
    title = "Motif Enrichment by Cell Type",
    subtitle = "Dot size = enrichment, color = significance"
  )


# found EGR1, EGR2, TFDP1, E2F6, SP1, SP3, KLF5, KLF6, ELK1, KLF2
# motifs relevant for HSC -> myeloid differentiation

important_myeloid_TFs <- c(
  "EGR1", "EGR2", "TFDP1", "E2F6", "SP1", "SP3", "KLF5", "KLF6", "KLF2", "ELK1",
  "RUNX1", "SPI1", "CEBPA", "GATA2", "GFI1", "IRF8", "STAT5A", "STAT5B",
  "MYB", "HOXA9", "NFIL3",
  "KLF4", "MAFB", "ZEB2", "ID2", "BCL11A"
)

filtered_df <- enriched_motif_df %>%
  filter(motif_name %in% important_myeloid_TFs)

filtered_df <- filtered_df %>%
  arrange(celltype, p_adjust) %>%
  group_by(celltype)


ggplot(filtered_df, aes(
  x = celltype,
  y = motif_name,
  size = enrichment,
  color = -log10(p_adjust)  # or use `color = p_adjust` for raw scale
)) +
  geom_point(alpha = 0.8) +
  scale_color_viridis_c(option = "C", direction = -1, name = "-log10(p-adjust)") +
  scale_size_continuous(name = "Enrichment") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.title = element_blank()
  ) +
  labs(
    title = "Motif Enrichment by Cell Type",
    subtitle = "Dot size = enrichment, color = significance"
  )

# Pivot to wide format: rows = motifs, columns = cell types
enrichment_matrix <- filtered_df %>%
  filter(celltype %in% celltype_array) %>%
  ungroup() %>%
  dplyr::select(motif_name, celltype, enrichment) %>%
  pivot_wider(names_from = celltype, values_from = enrichment)

# Set motif names as rownames
enrichment_mat <- as.data.frame(enrichment_matrix)
rownames(enrichment_mat) <- enrichment_mat$motif_name
enrichment_mat$motif_name <- NULL

enrichment_mat <- enrichment_mat %>%
  filter(if_all(everything(), ~ .x >= 1))

# Convert to matrix
enrichment_mat <- as.matrix(enrichment_mat)

motif_p <- Heatmap(
  enrichment_mat,
  name = "Motif Enrichment",
  top_annotation = col_anno,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE
)

# Draw the heatmap and return the full object
motif_p <- draw(motif_p)
# get the row ordering
row_order_vec <- row_order(motif_p)


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

# TF activity on shared peaks ----
# for TFs important for myeloid

# motif_sig <- unique(filtered_df$motif)
# motif_sig_names <- unlist(multiome_seurat[["ATAC"]]@motifs@motif.names)[motif_sig]

# multiome_seurat[["ATAC_shared"]] <- CreateAssayObject(counts = multiome_seurat@assays[["ATAC"]]$counts[peak_shared, ])
multiome_seurat[["ATAC_shared"]] <- CreateChromatinAssay(counts = multiome_seurat@assays[["ATAC"]]$counts[peak_shared, ])
multiome_seurat[["ATAC_shared"]]$data <- multiome_seurat@assays[["ATAC"]]$data[peak_shared, ]

pwm <- getMatrixSet(x = JASPAR2020,
                    opts = list(species = 9606, all_versions = FALSE)
                    )

# DefaultAssay(multiome_seurat) <- "ATAC"
multiome_seurat <- AddMotifs(object = multiome_seurat,
                       genome = BSgenome.Hsapiens.UCSC.hg38,
                       pfm = pwm,
                       assay = "ATAC_shared")

multiome_seurat <- RunChromVAR(multiome_seurat,
                               new.assay.name = "chromvar_shared",
                               # motif.matrix = multiome_seurat@assays[["ATAC"]]@motifs@data[peak_shared, ],
                               assay = "ATAC_shared",
                               genome = BSgenome.Hsapiens.UCSC.hg38)

motif_sig_names <- enrichment_matrix$motif_name[row_order_vec]
motif_sig <- names(multiome_seurat[["ATAC"]]@motifs@motif.names)[match(motif_sig_names, unlist(multiome_seurat[["ATAC"]]@motifs@motif.names))]

chromvar_mean_ls <- lapply(celltype_array, function(ct){
  id_cell <- which(Idents(multiome_seurat) == ct)
  chromvar_mean <- apply(multiome_seurat[["chromvar_shared"]]$data[motif_sig, id_cell], 1, mean)
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
        cluster_rows = FALSE,
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
  # motif_sig_names[row_order_vec] follows ordering from chromvar ordering above
  TF_expression_mean <- apply(expm1(multiome_seurat[["SCT"]]$data[motif_sig_names, id_cell]), 1, mean)
  return (log1p(TF_expression_mean))
})

TF_expression_mean_ls <- do.call(cbind, TF_expression_mean_ls)
colnames(TF_expression_mean_ls) <- celltype_array
rownames(TF_expression_mean_ls) <- motif_sig_names



pheatmap::pheatmap(TF_expression_mean_ls,
                   scale = "row",
                   cluster_cols = FALSE,
                   cluster_rows = FALSE,
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   annotation_names_col = FALSE,
                   annotation_col = anno_col)

Heatmap(TF_expression_mean_ls,
        cluster_rows = FALSE,
        name = "TF expr",  # This sets the legend title for the color scale
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

# check which TFs have most changes ----

which(abs(log2(TF_expression_mean_ls[, 1]/TF_expression_mean_ls[, 4])) >= 1)
which(abs(log2(abs(chromvar_mean_ls[, 1]/chromvar_mean_ls[, 4]))) >= 1)
which(abs(log2(enrichment_mat[, 1]/enrichment_mat[, 4])) >= 0.5)

# combine motif enrichment, TF activity and TF gene expression ----



library(tidyverse)

# Expression
expr_df <- as.data.frame(TF_expression_mean_ls) %>%
  rownames_to_column("TF") %>%
  pivot_longer(-TF, names_to = "celltype", values_to = "expression")

# chromVAR activity
activity_df <- as.data.frame(chromvar_mean_ls) %>%
  rownames_to_column("TF") %>%
  pivot_longer(-TF, names_to = "celltype", values_to = "activity")

# Motif enrichment
enrich_df <- as.data.frame(enrichment_mat) %>%
  rownames_to_column("TF") %>%
  pivot_longer(-TF, names_to = "celltype", values_to = "enrichment")

merged_df <- expr_df %>%
  left_join(activity_df, by = c("TF", "celltype")) %>%
  left_join(enrich_df, by = c("TF", "celltype"))

# Make sure celltype is treated as a factor
merged_df$celltype <- factor(merged_df$celltype, levels = c("1", "6", "7", "8"))

# TF expression vs motif
ggplot(merged_df, aes(x = expression, y = enrichment, color = celltype)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~celltype) +
  theme_minimal() +
  labs(title = "TF Expression vs Motif Enrichment", x = "Expression", y = "Motif Enrichment")

ggplot(merged_df, aes(x = expression, y = activity, color = celltype)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~celltype) +
  theme_minimal() +
  labs(title = "TF Activity (chromVAR) vs Expression", x = "Expression", y = "chromVAR Activity")

ggplot(merged_df, aes(x = enrichment, y = activity, color = celltype)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~celltype) +
  theme_minimal() +
  labs(title = "Motif Enrichment vs TF Activity", x = "Motif Enrichment", y = "chromVAR Activity")

cor(merged_df$expression, merged_df$activity)
cor(merged_df$expression, merged_df$enrichment)
cor(merged_df$activity, merged_df$enrichment)

cor.test(merged_df$expression, merged_df$enrichment)
cor.test(merged_df$expression, merged_df$activity)
cor.test(merged_df$activity, merged_df$enrichment)

# check gene functions in these shared peaks ----
library(GenomicRanges)

# Convert to GRanges object
peak_gr <- StringToGRanges(peak_shared, sep = c("-", "-"))


library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peak_annot <- annotatePeak(
  peak = peak_gr,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db",
  tssRegion = c(-3000, 3000),
  verbose = FALSE
)

plotAnnoPie(peak_annot)
plotAnnoBar(peak_annot)
plotDistToTSS(peak_annot)


# get nearest genes
nearest_genes <- as.data.frame(peak_annot)$geneId
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = nearest_genes,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# check cell cycle genes ----
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

cc_genes <- c(s.genes, g2m.genes)

intersect(gene_symbols, cc_genes)

# for shared static peaks
peak_annot_df <- as.data.frame(peak_annot)
peak_annot_df$peak_name <- paste0(peak_annot_df$seqnames, "-", peak_annot_df$start, "-", peak_annot_df$end)


peak_cc_genes <- peak_annot_df %>%
  filter(SYMBOL %in% cc_genes) %>%
  pull(peak_name)

peak_cc_genes_pair <- peak_annot_df %>%
  filter(SYMBOL %in% cc_genes) %>%
  pull(peak_name, SYMBOL)

peak_cc_genes_annotated <- peak_annot_df %>%
  filter(SYMBOL %in% cc_genes) %>%
  mutate(phase = case_when(
    SYMBOL %in% s.genes ~ "S",
    SYMBOL %in% g2m.genes ~ "G2M",
    TRUE ~ NA_character_  # Optional fallback
  ))

# check average of these peaks

peaks_mean_ls_cc_df <- peaks_mean_ls_df %>%
                          filter(peaks %in% peak_cc_genes_annotated$peak_name)

peaks_mean_ls_cc_df <- peaks_mean_ls_cc_df %>%
                        left_join(
                          peak_cc_genes_annotated %>% select(peak_name, phase),
                          by = c("peaks" = "peak_name")
                          )

ggplot(peaks_mean_ls_cc_df, aes(x = celltype, y = peak_accessibility, fill = phase)) +
  geom_boxplot() +
  geom_jitter() +
  theme_minimal() +
  labs(title = "Peaks near cell cycle genes",
       x = "Cell Type",
       y = "Mean Peak Accessibility")

cc_genes_shared <- intersect(cc_genes, rownames(multiome_seurat[["SCT"]]))

# S genes index
id_s <- match(s.genes, cc_genes_shared)
id_s <- id_s[!is.na(id_s)]

# G2M genes index
id_g2m <- match(g2m.genes, cc_genes_shared)
id_g2m <- id_g2m[!is.na(id_g2m)]

cc_genes_mean_list <- lapply(c(1, 6, 7, 8), function(ct){
  ct_id <- which(Idents(multiome_seurat) == ct)
  mean_gene_expr <- apply(expm1(multiome_seurat[["SCT"]]$data[cc_genes_shared, ct_id]), 1, mean)
  return (log1p(mean_gene_expr))
  })

cc_genes_mean_list <- do.call(cbind, cc_genes_mean_list)
names(cc_genes_mean_list) <- c(1, 6, 7, 8)

cc_genes_mean_list_df <- as.data.frame(cc_genes_mean_list)
colnames(cc_genes_mean_list_df) <- c(1, 6, 7, 8)
cc_genes_mean_list_df$phase <- rep("NA", nrow(cc_genes_mean_list_df))

cc_genes_mean_list_df$phase[id_s] <- "S"
cc_genes_mean_list_df$phase[id_g2m] <- "G2M"
cc_genes_mean_list_df$phase <- as.factor(cc_genes_mean_list_df$phase)

cc_genes_mean_list_df <- cc_genes_mean_list_df %>%
                          pivot_longer(cols = -phase,
                                       names_to = "cell_type",
                                       values_to = "gene_expression")

ggplot(cc_genes_mean_list_df, aes(x = cell_type, y = gene_expression, fill = phase)) +
  geom_boxplot() +
  geom_jitter() +
  theme_minimal() +
  labs(title = "Cell cycle genes for each phase",
       x = "Cell Type",
       y = "Mean Gene Expression")

# check correlation between mean peak accessibility and gene expression

for (ct in c(1, 6, 7, 8)){
  p_g_ct <- sapply(seq_along(peak_cc_genes_pair), function(id){
    gene_name <- names(peak_cc_genes_pair)[id]
    peak_name <- peak_cc_genes_pair[id]
    ct_id <- which(Idents(multiome_seurat) == ct)

    # cor_val <- cor(multiome_seurat[["ATAC_libnormed"]]$data[peak_name, ct_id],
    #                multiome_seurat[["SCT"]]$data[gene_name, ct_id])
    return (cor_val)
  })
}
# cause memory crash
# p_g_cor <- cor(t(as.numeric(multiome_seurat[["ATAC_libnormed"]]$data[peak_cc_genes_pair, ct_id])),
#                t(as.numeric(multiome_seurat[["SCT"]]$data[names(peak_cc_genes_pair), ct_id])))
# p_g_cor <- diag(p_g_cor)

# p_g_cor <- cor(t(peaks_mean_ls[unique(peaks_mean_ls_cc_df$peaks), ]), t(cc_genes_mean_list))

# pheatmap(p_g_cor)

# GO enrichment ----
library(clusterProfiler)

go_enrich <- enrichGO(
  gene = nearest_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # biological process
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

dotplot(go_enrich, showCategory = 10)

# filter for cell cycle process
cell_cycle_go <- go_enrich@result %>%
  dplyr::filter(grepl("cell cycle", Description, ignore.case = TRUE))  %>%

# Convert back to enrichResult-like object (optional but needed for dotplot)
go_enrich_subset <- go_enrich
go_enrich_subset@result <- cell_cycle_go

# Plot
dotplot(go_enrich_subset, showCategory = 10) +
  ggtitle("Selected GO Biological Processes")

autophagy_go <- go_enrich@result %>%
  dplyr::filter(grepl("autophagy", Description, ignore.case = TRUE))  # Case-insensitive match

# Get unique ENTREZ IDs involved in autophagy-related GO terms
autophagy_entrez <- unlist(strsplit(autophagy_go$geneID, split = "/"))
autophagy_entrez <- unique(autophagy_entrez)

# autophagy_entrez <- as.character(autophagy_entrez)
# autophagy_entrez <- autophagy_entrez[autophagy_entrez != "" & !is.na(autophagy_entrez)]
#
# valid_entrez <- keys(org.Hs.eg.db, keytype = "ENTREZID")

# autophagy_entrez_valid <- intersect(autophagy_entrez, valid_entrez)

# Convert to gene symbols
autophagy_symbols <- mapIds(
  org.Hs.eg.db,
  keys = autophagy_entrez,
  keytype = "SYMBOL",
  column = "ENTREZID",
  multiVals = "first"
)

autophagy_gene_shared <- intersect(autophagy_entrez, rownames(multiome_seurat[["SCT"]]$data))

celltype_array <- c(1, 6, 7, 8)
autophagy_expression_mean_ls <- lapply(celltype_array, function(ct){
  id_cell <- which(Idents(multiome_seurat) == ct)
  # motif_sig_names[row_order_vec] follows ordering from chromvar ordering above
  TF_expression_mean <- apply(expm1(multiome_seurat[["SCT"]]$data[autophagy_gene_shared, id_cell]), 1, mean)
  return (log1p(TF_expression_mean))
})

autophagy_expression_mean_ls <- do.call(cbind, autophagy_expression_mean_ls)
colnames(autophagy_expression_mean_ls) <- celltype_array
rownames(autophagy_expression_mean_ls) <- autophagy_gene_shared

autophagy_expression_mean_ls <- autophagy_expression_mean_ls[
  apply(autophagy_expression_mean_ls, 1, function(row) any(row >= 0.5)),
]

pheatmap::pheatmap(autophagy_expression_mean_ls,
                   scale = "none",
                   cluster_cols = FALSE,
                   cluster_rows = TRUE,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   annotation_names_col = FALSE,
                   annotation_col = anno_col)

Heatmap(autophagy_expression_mean_ls,
        cluster_rows = TRUE,
        name = "autophagy gene expr",  # This sets the legend title for the color scale
        top_annotation = col_anno,    # Add column annotation
        # col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        col = colorRampPalette(c("blue", "white", "red"))(100),
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_columns = FALSE)


id_autophagy <- which(abs(log2(autophagy_expression_mean_ls[, 1]/autophagy_expression_mean_ls[, 4])) >= 0.5)

VlnPlot(multiome_seurat,
        assay = "SCT",
        features = names(id_autophagy),
        idents = c(1, 6, 7, 8))

peak_anno_df <- as.data.frame(peak_anno)
peak_anno_df$peak_name <- paste0(peak_anno_df$seqnames, "-", peak_anno_df$start, "-", peak_anno_df$end)

peak_anno_df %>%
  filter(SYMBOL %in% c("MTDH", "BCL2", "ATG7", "LRBA", "ARHGAP26"))

peak_array <- peak_anno_df %>%
  # filter(SYMBOL %in% c("BCL2")) %>%
  # filter(SYMBOL %in% c("LRBA")) %>%
  dplyr::filter(SYMBOL %in% c("ATG7")) %>%
  dplyr::pull(peak_name)

library(EnsDb.Hsapiens.v86)

bcl2_range <- genes(EnsDb.Hsapiens.v86)[which(genes(EnsDb.Hsapiens.v86)$gene_name == "BCL2")]

bcl2_region <- resize(bcl2_range, width = width(bcl2_range) + 50000, fix = "center")

sample_array <- c("C1", "C2", "D", "T")
# for (i in seq_along(sample_array)){
#   multiome_seurat@assays[["ATAC"]]@fragments[[i]]@path <- paste0("~/luan/Juno/multiome_Rusland/multiome/data/raw/atac_fragments_", sample_array[i], ".tsv.gz")
# }

fragment_paths <- paste0("~/luan/Juno/multiome_Rusland/multiome/data/raw/atac_fragments_", sample_array, ".tsv.gz")

fragment_list <- lapply(seq_along(fragment_paths), function(path_id) {
  # remove C1 suffix to match the barcode in fragments file
  sample_cells <- colnames(multiome_seurat)[multiome_seurat@meta.data[["dataset"]] == sample_array[path_id]]
  sample_cells <- sub("^[^_]+_", "", sample_cells)

  CreateFragmentObject(
    path = fragment_paths[path_id],
    cells = sample_cells
  )
})

multiome_seurat[["ATAC"]]@fragments <- fragment_list

Fragments(multiome_seurat)

multiome_seurat[["ATAC"]]@fragments[[1]]@cells

Cells(multiome_seurat)

# can show normalized signals
# because the cell barcodes in Seurat differ by C1_ etc .. from fragments files
CoveragePlot(
  object = multiome_seurat,
  region = peak_array[2],
  idents = c(1, 6, 7, 8),         # only show selected cell types
  annotation = TRUE,
  peaks = TRUE,
  tile = FALSE
)

# plot regions in peak_array
peak_min_max <- peak_anno_df %>%
  # filter(SYMBOL %in% c("BCL2")) %>%
  # dplyr::filter(SYMBOL %in% c("LRBA")) %>%
  dplyr::filter(SYMBOL %in% c("ATG7")) %>%
  summarise(
    chr = seqnames,
    min_start = min(start),
    max_end = max(end)
  )

# for BCL2 gene
PeakPlot(
  object = multiome_seurat,
  region = paste0("chr18-", peak_min_max$min_start, "-", peak_min_max$max_end)
)

FindRegion <- function(
    object,
    region,
    sep = c("-", "-"),
    assay = NULL,
    extend.upstream = 0,
    extend.downstream = 0
) {
  if (!is(object = region, class2 = "GRanges")) {
    # first try to convert to coordinates, if not lookup gene
    region <- tryCatch(
      expr = suppressWarnings(
        expr = StringToGRanges(regions = region, sep = sep)
      ),
      error = function(x) {
        region <- LookupGeneCoords(
          object = object,
          assay = assay,
          gene = region
        )
        return(region)
      }
    )
    if (is.null(x = region)) {
      stop("Gene not found")
    }
  }
  region <- suppressWarnings(expr = Extend(
    x = region,
    upstream = extend.upstream,
    downstream = extend.downstream
  )
  )
  return(region)
}

# for BCL2 gene
region <- paste0("chr18-", peak_min_max$min_start, "-", peak_min_max$max_end)
# for LRBA gene and ATG7 gene
region <- paste0(peak_min_max$chr, "-", peak_min_max$min_start, "-", peak_min_max$max_end)

region <- FindRegion(object = multiome_seurat,
                     region = region,
                     sep = c("-", "-"),
                     assay = "ATAC",
                     extend.upstream = 0,
                     extend.downstream = 0)

library(GenomicRanges)
library(stringr)
# Split into chromosome, start, end
peak_split <- str_split(peak_array, "-", simplify = TRUE)

# Convert to GRanges
peaks <- GRanges(
  seqnames = peak_split[, 1],
  ranges = IRanges(
    start = as.numeric(peak_split[, 2]),
    end = as.numeric(peak_split[, 3])
  )
)

peak.intersect <- subsetByOverlaps(x = peaks, ranges = region)
peak.df <- as.data.frame(x = peak.intersect)
start.pos <- start(x = region)
end.pos <- end(x = region)
chromosome <- seqnames(x = region)

peak.df$start[peak.df$start < start.pos] <- start.pos
peak.df$end[peak.df$end > end.pos] <- end.pos
peak.plot <- ggplot(data = peak.df) +
                    geom_segment(aes(x = start, y = 0, xend = end, yend = 0), size = 5, data = peak.df)

peak.plot <- peak.plot + theme_classic() + ylab(label = "Peaks") +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  xlab(label = paste0(chromosome, " position (bp)")) +
  xlim(c(start.pos, end.pos))

peak.plot <- peak.plot + scale_color_manual(values = "dimgrey") +
  theme(legend.position = "none")

peak.plot

gene_plot <- AnnotationPlot(
  object = multiome_seurat,
  # region = "BCL2",
  # region = "LRBA"
  region = "ATG7"
)

gene_plot

CombineTracks(
  plotlist = list(peak.plot, gene_plot),
  heights = c(1, 1)
)

# for LRBA
# Define a common region: e.g., union of your peak and gene
common_start <- min(peak.df$start, 150264435)
common_end <- max(peak.df$end, 151015727)

# for ATG6
# based on peak.df
# common_start <- min(peak.df$start, 11277219)
# common_end <- max(peak.df$end, 11309028)

# based on AnnotationPlot
common_start <- min(peak.df$start, 11272309)
common_end <- max(peak.df$end, 11557665)



peak.plot <- peak.plot + xlim(common_start, common_end)

peak.plot


CombineTracks(
  plotlist = list(peak.plot, gene_plot),
  heights = c(1, 1)
)

VlnPlot(multiome_seurat,
        features = peak_array,
        idents = c(1, 6, 7, 8))


