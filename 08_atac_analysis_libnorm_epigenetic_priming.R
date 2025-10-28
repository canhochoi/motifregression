

# https://pmc.ncbi.nlm.nih.gov/articles/PMC7789351/#notes3
# epigentic priming suggests some peaks that are exclusively shared between
# blood stem cells and terminal cells remain accessible in both cell types


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

# result from shared static peaks
# save.image(file = "/data1/soldatr/luan/projects/motifregression/data/08_shared_static_peaks.RData")
load(file = "/data1/soldatr/luan/projects/motifregression/data/08_shared_static_peaks.RData")


# has multiome_seurat_noNA
load(file='/data1/soldatr/luan/projects/motifregression/result/CH/multiome_metacell.RData')

multiome_seurat <- readRDS("~/luan/Juno/multiome_Rusland/multiome/data/processed/multiome_seurat.rds")

DefaultAssay(multiome_seurat) <- "ATAC"

# find shared peaks with similar accessibility in cell type 1 and 8 ----
# another way
multiome_seurat <- normalize_library_size(multiome_seurat, lib_size = 1e4)

# find shared peaks from motifregression package
peak_shared <- find_shared_peaks(multiome_seurat, 1, 8)

# remove chrM
peak_shared <- peak_shared[!grepl("^chrM-", peak_shared)]

# check if these peaks are statistically insignificant

peak_shared_de <- FindMarkers(multiome_seurat,
                              assay = "ATAC",
                              features = peak_shared,
                              ident.1 = 1,
                              ident.2 = 8)

ggplot(peak_shared_de, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point() +
  xlim(c(-1, 1)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(x = "log2 FC", y = "-log(p_adj)", title = "shared_peak_1_8")



# number of DA peaks that statistically insignificant
peak_shared_de %>% filter(p_val_adj <= 0.05) %>% nrow()

DA_peaks_shared <- peak_shared_de %>% filter(p_val_adj <= 0.05) %>% rownames()

# check gene functions in these shared peaks ----
library(GenomicRanges)

# Convert to GRanges object
peak_gr <- StringToGRanges(peak_shared, sep = c("-", "-"))

length(peak_gr) == length(peak_shared)

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

# convert to data frame
peak_annot_df <- as.data.frame(peak_annot)
peak_annot_df$peak_name <- paste0(peak_annot_df$seqnames, "-", peak_annot_df$start, "-", peak_annot_df$end)

nrow(peak_annot_df) == length(peak_shared)

peak_annot_df %>% filter(peak_name %in% DA_peaks_shared) %>% pull(annotation) %>% table()

peak_annot_df %>% filter(peak_name %in% DA_peaks_shared) %>% pull(SYMBOL)

# check if these peaks are accessible in other terminal cell types

# 10 is erythroid
peak_shared_1_10 <- find_shared_peaks(multiome_seurat, 1, 10)
# 4 is lymphoid
peak_shared_1_4 <- find_shared_peaks(multiome_seurat, 1, 4)
# 5 is DC
peak_shared_1_5 <- find_shared_peaks(multiome_seurat, 1, 5)

combined_peaks_all <- unique(c(peak_shared_1_10, peak_shared_1_4, peak_shared_1_5))

# plot Venn diagram
library(ggvenn)

venn_list <- list(
  '1_10' = peak_shared_1_10,
  '1_5' = peak_shared_1_5,
  '1_4' = peak_shared_1_4,
  '1_8' = peak_shared
)

ggvenn(venn_list, fill_color = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"))

# find shared peaks and corresponding genes
shared_peaks_across_ct <- Reduce(intersect, venn_list)

diff_shared_peaks_da <- FindAllMarkers(multiome_seurat,
                                assay = "ATAC",
                                features = shared_peaks_across_ct)

# Filter for selected clusters and significance
filtered_df <- diff_shared_peaks_da %>%
  filter(cluster %in% c(4, 5, 8, 10)) %>%
  mutate(
    neg_log10_padj = -log10(p_val_adj),
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 1 ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -1 ~ "Down",
      TRUE ~ "Not Significant"
    )
  )

ggplot(filtered_df, aes(x = avg_log2FC, y = neg_log10_padj, color = significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Significant" = "gray")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  facet_wrap(~ cluster, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Volcano Plot of Shared Peaks Across Clusters",
    x = "log2 Fold Change",
    y = "-log10 Adjusted p-value"
  )


diff_peaks_da <- FindAllMarkers(multiome_seurat,
                                assay = "ATAC")



# find genes whose promoters contain shared peaks
genes_of_shared_peaks_across_ct <- peak_annot_df %>%
  dplyr::filter(!(annotation == "Distal Intergenic")) %>%
  dplyr::filter(peak_name %in% shared_peaks_across_ct) %>%
  pull(SYMBOL)
# pull(geneId)

genes_of_shared_peaks_across_ct <- intersect(genes_of_shared_peaks_across_ct, rownames(multiome_seurat[["SCT"]]$data))

# find DE genes
diff_genes_de <- FindAllMarkers(multiome_seurat,
                                assay = "SCT",
                                features = genes_of_shared_peaks_across_ct)

# Define cell types of interest
celltypes_to_plot <- c(4, 5, 8, 10)

# Filter the data for selected cell types
filtered_df <- diff_genes_de %>%
  filter(cluster %in% celltypes_to_plot)

# Create a new column for significance coloring
filtered_df <- filtered_df %>%
  mutate(sig = p_val_adj < 0.05 & abs(avg_log2FC) >= 1)

# Plot with faceting
ggplot(filtered_df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ cluster, scales = "free") +
  scale_color_manual(values = c("grey", "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  labs(
    title = "Volcano Plots for Cell Types 4, 5, 8, and 10 vs the rest",
    x = "avg_log2FC",
    y = "-log10(adj. p-value)"
  ) +
  theme_minimal()


# Get unique genes per cluster
de_genes_list <- filtered_df %>%
  filter(cluster %in% c(4, 5, 8, 10), p_val_adj <= 0.05, abs(avg_log2FC) >= 1) %>%
  group_by(cluster) %>%
  summarise(genes = list(unique(gene))) %>%
  tibble::deframe()  # Named list: names = cluster, values = vector of genes

library(ggvenn)
ggvenn(de_genes_list, fill_color = c("#8da0cb",  "#e78ac3", "#fc8d62", "#66c2a5"))

# plot heatmap of all DE genes
de_genes <- unique(unlist(de_genes_list))

de_genes_mean_expr <- get_mean_expression_for_heatmap(multiome_seurat, celltype_array, de_genes, assay = "SCT")

anno_col <- data.frame(celltype = celltype_array)

anno_col$celltype <- as.factor(anno_col$celltype)
#important to get correct annotation
rownames(anno_col) <- colnames(de_genes_mean_expr)

highly_expr_genes_8 <- rownames(de_genes_mean_expr)[de_genes_mean_expr[, 4] >= 1]

highly_expr_de_genes_8 <- intersect(highly_expr_genes_8, de_genes_list[["8"]])

de_genes_mean_expr[highly_expr_de_genes_8, ]

presence_ls <- lapply(c(4, 8, 5, 10), function(ct){
  highly_expr_de_genes_8 %in% de_genes_list[[as.character(ct)]]
})

presence_ls <- do.call(rbind, presence_ls)
rownames(presence_ls) <- c(4, 8, 5, 10)
colnames(presence_ls) <- highly_expr_de_genes_8
presence_ls

# Create custom row labels
custom_labels <- rownames(de_genes_mean_expr)
custom_labels[!custom_labels %in% highly_expr_de_genes_8] <- ""

# Plot the full matrix with custom row labels
pheatmap::pheatmap(de_genes_mean_expr,
                   scale = "none",
                   cluster_cols = FALSE,
                   cluster_rows = TRUE,
                   show_rownames = TRUE,
                   labels_row = custom_labels,  # use custom labels
                   show_colnames = FALSE,
                   annotation_names_col = FALSE,
                   annotation_col = anno_col)

Idents(multiome_seurat) <- factor(Idents(multiome_seurat), levels = sort(as.numeric(levels(Idents(multiome_seurat)))))
VlnPlot(multiome_seurat,
        assay = "SCT",
        features = highly_expr_de_genes_8,
        idents = c(1, 6, 7, 8, 4, 5, 10))

# find peaks of these genes CCSER1 and LYST

peak_of_highly_expr_de_genes_8 <- peak_annot_df %>%
  dplyr::filter(SYMBOL %in% highly_expr_de_genes_8) %>%
  mutate(order = match(SYMBOL, highly_expr_de_genes_8)) %>%
  arrange(order) %>%
  dplyr::pull(peak_name)

VlnPlot(multiome_seurat,
        assay = "ATAC_libnormed",
        features = peak_of_highly_expr_de_genes_8,
        idents = c(1, 6, 7, 8, 4, 5, 10))


# find genes uniquely in 8

unique_de_genes_8 <- setdiff(de_genes_list[["8"]], c(de_genes_list[["4"]], de_genes_list[["5"]], de_genes_list[["10"]]))

unique_de_genes_8

# plot heatmap of these gene expression

celltype_array <- c(1, 6, 7, 8)
unique_de_genes_8_mean_expr <- get_mean_expression_for_heatmap(multiome_seurat, celltype_array, unique_de_genes_8, assay = "SCT")


anno_col <- data.frame(celltype = celltype_array)

anno_col$celltype <- as.factor(anno_col$celltype)
#important to get correct annotation
rownames(anno_col) <- colnames(unique_de_genes_8_mean_expr)

# Create custom row labels: only label MDK and PCSK5, others are blank
custom_labels <- rownames(unique_de_genes_8_mean_expr)
custom_labels[!custom_labels %in% c("FOS", "FOSB", "MAFG", "ATG7", "CALR", "CLEC11A", "PIK3CB", "CAT", "MAP3K1")] <- ""

# Plot the full matrix with custom row labels
pheatmap::pheatmap(unique_de_genes_8_mean_expr,
                   scale = "none",
                   cluster_cols = FALSE,
                   cluster_rows = TRUE,
                   show_rownames = TRUE,
                   labels_row = custom_labels,  # use custom labels
                   show_colnames = FALSE,
                   annotation_names_col = FALSE,
                   annotation_col = anno_col)

non_distal_peaks_df <- peak_annot_df %>%
  dplyr::filter(!(annotation == "Distal Intergenic"))

shared_nearest_genes <- non_distal_peaks_df %>%
  dplyr::filter(SYMBOL %in% unique_de_genes_8) %>%
  # pull(SYMBOL)
  pull(geneId)

shared_nearest_genes

# find how many DA peaks between 1 and 8 is in shared_peaks_across_ct
# there are 80 DA peaks
intersect(DA_peaks_shared, shared_peaks_across_ct)

# find static peaks uniquely in myeloid ----

myeloid_unique_peaks <- setdiff(peak_shared, combined_peaks_all)

length(myeloid_unique_peaks)

myeloid_unique_peaks_gene <- peak_annot_df %>%
  dplyr::filter(!(annotation == "Distal Intergenic")) %>%
  dplyr::filter((peak_name %in% myeloid_unique_peaks)) %>%
  pull(geneId)

length(myeloid_unique_peaks_gene)

myeloid_gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = myeloid_unique_peaks_gene,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

myeloid_gene_symbols <- intersect(myeloid_gene_symbols, rownames(multiome_seurat[["SCT"]]$data))

myeloid_gene_de <- FindMarkers(multiome_seurat,
                               assay = "SCT",
                               features = myeloid_gene_symbols,
                               ident.1 = 8,
                               ident.2 = 1)

myeloid_gene_de <- myeloid_gene_de %>%
  tibble::rownames_to_column(var = "gene")

myeloid_gene_de <- myeloid_gene_de %>%
  mutate(sig = ifelse(abs(avg_log2FC) >= 1 & p_val_adj <= 0.05, "Significant", "Not Significant"))

sig_genes <- myeloid_gene_de %>%
  filter(sig == "Significant")

library(ggrepel)

ggplot(myeloid_gene_de, aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal() +
  geom_text_repel(
    data = sig_genes,
    aes(label = gene),
    size = 3,
    max.overlaps = 100,
    box.padding = 0.4
  ) +
  labs(title = "8 vs 1",
       x = "Avg_log2_FC",
       y = "-log10(Adj_P_value)",
       color = "Significance")


gene_expression_mean_ls <- lapply(celltype_array, function(ct){
  id_cell <- which(Idents(multiome_seurat) == ct)
  # motif_sig_names[row_order_vec] follows ordering from chromvar ordering above
  gene_expression_mean <- apply(expm1(multiome_seurat[["SCT"]]$data[myeloid_gene_symbols, id_cell]), 1, mean)
  return (log1p(gene_expression_mean))
})

gene_expression_mean_ls <- do.call(cbind, gene_expression_mean_ls)
colnames(gene_expression_mean_ls) <- celltype_array
rownames(gene_expression_mean_ls) <- myeloid_gene_symbols


anno_col <- data.frame(celltype = celltype_array)

anno_col$celltype <- as.factor(anno_col$celltype)
#important to get correct annotation
rownames(anno_col) <- colnames(gene_expression_mean_ls)

# Create custom row labels: only label MDK and PCSK5, others are blank
custom_labels <- rownames(gene_expression_mean_ls)
custom_labels[!custom_labels %in% c("MDK", "PCSK5")] <- ""

# Plot the full matrix with custom row labels
pheatmap::pheatmap(gene_expression_mean_ls,
                   scale = "none",
                   cluster_cols = FALSE,
                   cluster_rows = TRUE,
                   show_rownames = TRUE,
                   labels_row = custom_labels,  # use custom labels
                   show_colnames = FALSE,
                   annotation_names_col = FALSE,
                   annotation_col = anno_col)


# check shared peaks among terminal cells and stem cells ----
non_distal_peaks_df <- peak_annot_df %>%
  dplyr::filter(!(annotation == "Distal Intergenic"))

non_distal_peaks_shared <- non_distal_peaks_df %>%
  dplyr::filter((peak_name %in% shared_peaks_across_ct)) %>%
  pull(peak_name)

shared_nearest_genes <- non_distal_peaks_df %>%
  dplyr::filter(peak_name %in% non_distal_peaks_shared) %>%
  # pull(SYMBOL)
  pull(geneId)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = shared_nearest_genes,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# do GO enrichment
library(clusterProfiler)

go_enrich <- enrichGO(
  gene = shared_nearest_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # biological process
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

dotplot(go_enrich, showCategory = 15)

library(enrichplot)
library(dplyr)

# Sort GO terms by decreasing GeneRatio
go_enrich_sorted <- go_enrich@result %>%
  mutate(GeneRatio_num = sapply(GeneRatio, function(x) eval(parse(text = x)))) %>%
  arrange(desc(GeneRatio_num))

# Select top N with highest GeneRatio
top_go <- go_enrich_sorted %>%
  arrange(desc(GeneRatio_num)) %>%
  head(20)

# Plot manually using ggplot2
ggplot(top_go, aes(x = GeneRatio_num, y = reorder(Description, GeneRatio_num))) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_color_viridis_c(option = "C", direction = -1) +
  theme_minimal() +
  labs(
    x = "Gene Ratio",
    y = "GO Biological Process",
    title = "Top 20 Enriched GO Terms (by Gene Ratio)"
  )



differentiation_go <- go_enrich@result %>%
  dplyr::filter(grepl("differentiation", Description, ignore.case = TRUE) & p.adjust <= 0.05)

differentiation_go <- differentiation_go %>%
                        mutate(GeneRatio_num = sapply(GeneRatio, function(x) eval(parse(text = x))))

# Plot manually using ggplot2
ggplot(differentiation_go, aes(x = GeneRatio_num, y = reorder(Description, GeneRatio_num))) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_color_viridis_c(option = "C", direction = -1) +
  theme_minimal() +
  labs(
    x = "Gene Ratio",
    y = "GO Biological Process",
    title = "GO for differentiation"
  )


cell_cycle_go <- go_enrich@result %>%
  dplyr::filter(grepl("cell cycle", Description, ignore.case = TRUE) & p.adjust <= 0.05)

# Find GO terms with "autophagy" in their description
autophagy_terms <- go_enrich@result %>%
  filter(grepl("autophagy", Description, ignore.case = TRUE))

# Extract genes from each autophagy-related GO term
autophagy_genes <- unique(unlist(strsplit(autophagy_terms$geneID, "/")))

autophagy_genes_shared <- intersect(autophagy_genes, rownames(multiome_seurat[["SCT"]]$data))


# calculate peaks accessibility of shared peaks in the promoters ----
celltype_array <- c(1, 6, 7, 8)
celltype_array <- c(1, 6, 7, 8, 4, 5, 10)
shared_peak_expression_mean_ls <- lapply(celltype_array, function(ct){
  id_cell <- which(Idents(multiome_seurat) == ct)
  # motif_sig_names[row_order_vec] follows ordering from chromvar ordering above
  # only normalize by lib size so no need to expm1 and log1p
  peak_expression_mean <- apply((multiome_seurat[["ATAC_libnormed"]]$data[non_distal_peaks_shared, id_cell]), 1, mean)
  return ((peak_expression_mean))
})

shared_peak_expression_mean_ls <- do.call(cbind, shared_peak_expression_mean_ls)
colnames(shared_peak_expression_mean_ls) <- celltype_array
rownames(shared_peak_expression_mean_ls) <- non_distal_peaks_shared


# anno_col <- data.frame(celltype = c(1, 6, 7, 8))

anno_col <- data.frame(celltype = celltype_array)

anno_col$celltype <- as.factor(anno_col$celltype)
#important to get correct annotation
rownames(anno_col) <- colnames(shared_peak_expression_mean_ls)

pheatmap::pheatmap(shared_peak_expression_mean_ls,
                   # scale = "row",
                   scale = "none",
                   cluster_cols = FALSE,
                   cluster_rows = TRUE,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   annotation_names_col = FALSE,
                   annotation_col = anno_col)


pheatmap::pheatmap(shared_peak_expression_mean_ls,
                   # scale = "row",
                   scale = "row",
                   cluster_cols = TRUE,
                   cluster_rows = TRUE,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   annotation_names_col = FALSE,
                   annotation_col = anno_col)




# calculate how many lymphoid genes are differentially expressed in lymphoid ----
# filter out gene sets involved with myeloid
myeloid_go <- go_enrich@result %>%
  dplyr::filter(grepl("myeloid", Description, ignore.case = TRUE) & p.adjust <= 0.05)

myeloid_diff_genes <- myeloid_go %>%
  dplyr::filter(grepl("differentiation", Description, ignore.case = TRUE) & p.adjust <= 0.05) %>%
  pull(geneID) %>%
  strsplit(split = "/") %>%
  unlist()

myeloid_diff_genes_in_multiome <- intersect(myeloid_diff_genes, rownames(multiome_seurat[["SCT"]]$data))

Idents(multiome_seurat) <- multiome_seurat$projection
multiome_seurat <- PrepSCTFindMarkers(multiome_seurat, assay = "SCT", verbose = TRUE)

myeloid_diff_genes_de_ls <- lapply(c(4, 5, 8, 10), function(ct){
  myeloid_diff_genes_de <- FindMarkers(multiome_seurat,
                                       assay = "SCT",
                                       features = myeloid_diff_genes_in_multiome,
                                       ident.1 = ct,
                                       ident.2 = 1)
  # myeloid_diff_genes_de_genes <- myeloid_diff_genes_de %>%
  #                                  dplyr::filter(pct.1 >= 0.05 & pct.2 >= 0.05 & p_val_adj <= 0.05 & abs(avg_log2FC) >= 1) %>%
  #                                  rownames()

  # return (myeloid_diff_genes_de_genes)

  return (myeloid_diff_genes_de)

})

names(myeloid_diff_genes_de_ls) <- c(4, 5, 8, 10)

for (ct in c(4, 5, 8, 10)){
  myeloid_diff_genes_de <- myeloid_diff_genes_de_ls[[as.character(ct)]]
  ggplot(myeloid_diff_genes_de, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point() +
    # xlim(c(-1, 1)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    labs(x = "log2 FC", y = "-log(p_adj)", title = paste0(ct, " vs 1"))

}

Reduce(intersect, myeloid_diff_genes_de_ls)

ggvenn(myeloid_diff_genes_de_ls, fill_color = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"))


# calculate gene expression of autophagy genes whose promoters contain shared peaks ----


# gene_symbols_shared <- intersect(gene_symbols, rownames(multiome_seurat[["SCT"]]$data))

autophagy_gene_expression_mean_ls <- lapply(celltype_array, function(ct){
  id_cell <- which(Idents(multiome_seurat) == ct)
  # motif_sig_names[row_order_vec] follows ordering from chromvar ordering above
  gene_expression_mean <- apply(expm1(multiome_seurat[["SCT"]]$data[autophagy_genes_shared, id_cell]), 1, mean)
  return (log1p(gene_expression_mean))
})

autophagy_gene_expression_mean_ls <- do.call(cbind, autophagy_gene_expression_mean_ls)
colnames(autophagy_gene_expression_mean_ls) <- celltype_array
rownames(autophagy_gene_expression_mean_ls) <- autophagy_genes_shared


autophagy_1_8 <- FindMarkers(multiome_seurat,
            assay = "SCT",
            features = autophagy_genes_shared,
            ident.1 = 1,
            ident.2 = 8)

autophagy_1_8 <- autophagy_1_8 %>%
                    filter(pct.1 > 0.1 & pct.2 > 0.1 & p_val_adj <= 0.05)

autophagy_1_4 <- FindMarkers(multiome_seurat,
                             assay = "SCT",
                             features = autophagy_genes_shared,
                             ident.1 = 1,
                             ident.2 = 4)

autophagy_1_4 <- autophagy_1_4 %>%
  filter(pct.1 > 0.1 & pct.2 > 0.1 & p_val_adj <= 0.05)


autophagy_1_5 <- FindMarkers(multiome_seurat,
                              assay = "SCT",
                              features = autophagy_genes_shared,
                              ident.1 = 1,
                              ident.2 = 5)

autophagy_1_5 <- autophagy_1_5 %>%
  filter(pct.1 > 0.1 & pct.2 > 0.1 & p_val_adj <= 0.05)

# try to see if can use autophagy genes to cluster cell types
combined_diff_genes <- unique(c(rownames(autophagy_1_10),
                        rownames(autophagy_1_4),
                        rownames(autophagy_1_5),
                        rownames(autophagy_1_8)))

celltype_array_combined <- levels(Idents(multiome_seurat))
autophagy_gene_expression_mean_ls <- lapply(celltype_array_combined, function(ct){
  id_cell <- which(Idents(multiome_seurat) == ct)
  # motif_sig_names[row_order_vec] follows ordering from chromvar ordering above
  gene_expression_mean <- apply(expm1(multiome_seurat[["SCT"]]$data[combined_diff_genes, id_cell]), 1, mean)
  return (log1p(gene_expression_mean))
})

autophagy_gene_expression_mean_ls <- do.call(cbind, autophagy_gene_expression_mean_ls)
colnames(autophagy_gene_expression_mean_ls) <- celltype_array_combined
rownames(autophagy_gene_expression_mean_ls) <- combined_diff_genes

anno_col_shared <- data.frame(celltype = celltype_array_combined)
anno_col_shared$celltype <- as.factor(anno_col_shared$celltype)
#important to get correct annotation
rownames(anno_col_shared) <- colnames(autophagy_gene_expression_mean_ls)

# see if cell types can be clustered
pheatmap::pheatmap(autophagy_gene_expression_mean_ls,
                   # scale = "row",
                   scale = "row",
                   cluster_cols = TRUE,
                   cluster_rows = TRUE,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   annotation_names_col = FALSE,
                   annotation_col = anno_col_shared)


shared_diff_genes <- Reduce(intersect,
                     list(rownames(autophagy_1_10),
                          rownames(autophagy_1_4),
                          rownames(autophagy_1_5),
                          rownames(autophagy_1_8))
              )

autophagy_de_genes_1_8 <- setdiff(rownames(autophagy_1_8), shared_diff_genes)
autophagy_de_genes_1_10 <- setdiff(rownames(autophagy_1_10), shared_diff_genes)
autophagy_de_genes_1_4 <- setdiff(rownames(autophagy_1_4), shared_diff_genes)
autophagy_de_genes_1_5 <- setdiff(rownames(autophagy_1_5), shared_diff_genes)



pheatmap::pheatmap(autophagy_gene_expression_mean_ls[autophagy_de_genes_1_8, ],
                   # scale = "row",
                   scale = "none",
                   cluster_cols = FALSE,
                   cluster_rows = FALSE,
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   annotation_names_col = FALSE,
                   annotation_col = anno_col)

pheatmap::pheatmap(autophagy_gene_expression_mean_ls[autophagy_de_genes_1_8, ],
                   scale = "row",
                   # scale = "none",
                   cluster_cols = FALSE,
                   cluster_rows = TRUE,
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   annotation_names_col = FALSE,
                   annotation_col = anno_col)

VlnPlot(multiome_seurat,
        assay = "SCT",
        features = shared_diff_genes,
        idents = c(1, 6, 7, 8, 4, 5, 10))


# among non-distal peaks, find exclusively peaks shared by stem cell and myeloid ----
non_distal_peaks_df <- peak_annot_df %>%
                    dplyr::filter(!(annotation == "Distal Intergenic"))

non_distal_peaks_1_8 <- non_distal_peaks_df %>%
  dplyr::filter(!(peak_name %in% combined_peaks_all)) %>%
  pull(peak_name)

length(non_distal_peaks_1_8)

# get nearest genes of uniquely shared peaks between 1 and 8 ----

nearest_genes <- non_distal_peaks_df %>%
  dplyr::filter(peak_name %in% non_distal_peaks_1_8) %>%
  # pull(SYMBOL)
  pull(geneId)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = nearest_genes,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# do GO enrichment
# library(clusterProfiler)
#
# go_enrich <- enrichGO(
#   gene = nearest_genes,
#   OrgDb = org.Hs.eg.db,
#   keyType = "ENTREZID",
#   ont = "BP",  # biological process
#   pAdjustMethod = "BH",
#   qvalueCutoff = 0.05,
#   readable = TRUE
# )

# dotplot(go_enrich, showCategory = 10)

# filter out gene sets involved with myeloid
myeloid_go <- go_enrich@result %>%
  dplyr::filter(grepl("myeloid", Description, ignore.case = TRUE) & p.adjust <= 0.05)

# extract out the genes

myeloid_gene <- sapply(myeloid_go$geneID, function(gene_str) {strsplit(gene_str, "/")})
myeloid_gene <- unique(unlist(myeloid_gene))

# calculate peaks accessibility of the peaks in the promoters ----
celltype_array <- c(1, 6, 7, 8)
peak_expression_mean_ls <- lapply(celltype_array, function(ct){
  id_cell <- which(Idents(multiome_seurat) == ct)
  # motif_sig_names[row_order_vec] follows ordering from chromvar ordering above
  # only normalize by lib size so no need to expm1 and log1p
  peak_expression_mean <- apply((multiome_seurat[["ATAC_libnormed"]]$data[non_distal_peaks_1_8, id_cell]), 1, mean)
  return ((peak_expression_mean))
})

peak_expression_mean_ls <- do.call(cbind, peak_expression_mean_ls)
colnames(peak_expression_mean_ls) <- celltype_array
rownames(peak_expression_mean_ls) <- non_distal_peaks_1_8

anno_col <- data.frame(celltype = c(1, 6, 7, 8))
anno_col$celltype <- as.factor(anno_col$celltype)
#important to get correct annotation
rownames(anno_col) <- colnames(peak_expression_mean_ls)

pheatmap::pheatmap(peak_expression_mean_ls,
                   # scale = "row",
                   scale = "none",
                   cluster_cols = FALSE,
                   cluster_rows = TRUE,
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   annotation_names_col = FALSE,
                   annotation_col = anno_col)


# calculate gene expression of genes of the peaks in the promoters ----

gene_symbols_shared <- intersect(gene_symbols, rownames(multiome_seurat[["SCT"]]$data))

gene_symbols_shared <- intersect(myeloid_gene, rownames(multiome_seurat[["SCT"]]$data))


myeloid_gene_de <- FindMarkers(multiome_seurat,
                               assay = "SCT",
                               features = gene_symbols_shared,
                               ident.1 = 1,
                               ident.2 = 8)

gene_expression_mean_ls <- lapply(celltype_array, function(ct){
  id_cell <- which(Idents(multiome_seurat) == ct)
  # motif_sig_names[row_order_vec] follows ordering from chromvar ordering above
  gene_expression_mean <- apply(expm1(multiome_seurat[["SCT"]]$data[gene_symbols_shared, id_cell]), 1, mean)
  return (log1p(gene_expression_mean))
})

gene_expression_mean_ls <- do.call(cbind, gene_expression_mean_ls)
colnames(gene_expression_mean_ls) <- celltype_array
rownames(gene_expression_mean_ls) <- gene_symbols_shared




pheatmap::pheatmap(gene_expression_mean_ls,
                   scale = "row",
                   # scale = "none",
                   cluster_cols = FALSE,
                   cluster_rows = TRUE,
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   annotation_names_col = FALSE,
                   annotation_col = anno_col)

VlnPlot(multiome_seurat,
        assay = "SCT",
        idents = c(1, 6, 7, 8),
        features = c("PHF19", "NPRL3", "CORO2A", "RALGDS", "SH2B2"),
        # features = c("KLF2", "ZFP36", "JAG1", "HK2")
        )


# do motif enrichment ----
celltype_array <- c(1, 6, 7, 8)

# enriched_motif_celltype is in motifregression package
# atac_utils.R
enriched_motif_ls <- lapply(celltype_array, function(ct){
  enriched_motif <- enriched_motif_celltype(multiome_seurat, ct, non_distal_peaks_1_8)
  enriched_motif$celltype <- rep(ct, nrow(enriched_motif))
  return (enriched_motif)
})

names(enriched_motif_ls) <- celltype_array


# among distal peaks, find exclusively peaks shared by stem cell and myeloid ----
# find distal intergenic
distal_peaks_df <- peak_annot_df %>%
  dplyr::filter(annotation == "Distal Intergenic")

distal_peaks_1_8 <- distal_peaks_df %>%
  filter(!(peak_name %in% combined_peaks_all)) %>%
  pull(peak_name)

celltype_array <- c(1, 6, 7, 8)

# enriched_motif_celltype is in motifregression package
# atac_utils.R
enriched_motif_ls <- lapply(celltype_array, function(ct){
  enriched_motif <- enriched_motif_celltype(multiome_seurat, ct, distal_peaks_1_8)
  enriched_motif$celltype <- rep(ct, nrow(enriched_motif))
  return (enriched_motif)
})

names(enriched_motif_ls) <- celltype_array


