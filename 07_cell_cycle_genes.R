devtools::load_all("~/luan/projects/motifregression/")
devtools::document()

library(Signac)
library(Seurat)
library(RCurl)
library(AnnotationHub)
library(dplyr)

# has multiome_seurat_noNA
load(file='/data1/soldatr/luan/projects/motifregression/result/CH/multiome_metacell.RData')

multiome_seurat <- readRDS("~/luan/Juno/multiome_Rusland/multiome/data/processed/multiome_seurat.rds")


# get cell cycle genes list ----
#https://hbctraining.github.io/scRNA-seq/lessons/cell_cycle_scoring.html

# more genes with G1 phase
# table S1 of https://celldiv.biomedcentral.com/articles/10.1186/s13008-019-0058-4#MOESM1
# and this one
# but have only G1/S 
# table S3 of https://academic.oup.com/jmcb/article/11/8/703/5188008#supplementary-data

cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv")
cell_cycle_genes <- read.csv(text = cc_file)

# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah,
              pattern = c("Homo sapiens", "EnsDb"),
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb,
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)


# Get gene names for Ensembl IDs for each gene
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")

# Perform cell cycle scoring ----

DefaultAssay(multiome_seurat) <- "SCT"
multiome_seurat <- CellCycleScoring(multiome_seurat,
                                   g2m.features = g2m_genes,
                                   s.features = s_genes)

# Perform PCA and color by cell cycle phase ----
# https://satijalab.org/seurat/archive/v3.0/cell_cycle_vignette.html

DefaultAssay(multiome_seurat) <- "SCT"
multiome_seurat <- ScaleData(multiome_seurat,
                             features = rownames(multiome_seurat))

multiome_seurat <- RunPCA(multiome_seurat,
                          assay = "SCT",
                          features = c(g2m_genes, s_genes))

# Visualize the PCA, grouping by cell cycle phase
p1 <- DimPlot(multiome_seurat,
        reduction = "pca",
        group.by= "Phase")

p1 + labs(title = "Cell cycle genes")

#try peaks linked to cell cycle genes to separate cells ----

# multiome_seurat <- LinkPeaks_modified(multiome_seurat,
#                                      peak.assay = "ATAC",
#                                      expression.assay = "SCT",
#                                      genes.use = c(g2m_genes, s_genes))

multiome_seurat <- LinkPeaks(multiome_seurat,
                            peak.assay = "ATAC",
                            expression.assay = "SCT",
                            genes.use = c(g2m_genes, s_genes))


DefaultAssay(multiome_seurat) <- "ATAC"

link <- Links(multiome_seurat)
link <- as.data.frame(link)
link$adj_pval <- p.adjust(link$pvalue, method = 'BH')


head(link)

multiome_seurat <- ScaleData(multiome_seurat,
                             assay = "ATAC",
                             features = link$peak)

multiome_seurat <- RunPCA(multiome_seurat,
                          assay = "ATAC",
                          features = link$peak)

p2 <- DimPlot(multiome_seurat,
        reduction = "pca",
        group.by= "Phase")

p2 + labs(title = "Peaks linked to genes")

cellcycle_common <- intersect(c(g2m_genes, s_genes), rownames(multiome_seurat[["Promoter"]]$data))


# try promoter activity to separate cell cycle genes ----
multiome_seurat <- NormalizeData(multiome_seurat,
                                 assay = "Promoter",
                                 normalization.method = 'LogNormalize',
                                 scale.factor = median(multiome_seurat$nCount_Promoter))

multiome_seurat <- ScaleData(multiome_seurat,
                             assay = "Promoter",
                             features = cellcycle_common)

multiome_seurat <- RunPCA(multiome_seurat,
                          assay = "Promoter",
                          features = cellcycle_common)

p3 <- DimPlot(multiome_seurat,
        reduction = "pca",
        group.by= "Phase")

p3 + labs(title = "Gene activity via promoter")


#



library(scMEGA)
library(SummarizedExperiment)

p2g_df <- PeakToGene(
  as.matrix(multiome_seurat[["ATAC"]]$data),
  as.matrix(multiome_seurat[["SCT"]]$data[c(g2m_genes, s_genes), ]),
  genome = "hg38",
  max.dist = 1e6,
  method = "correlation"
)


#try differential peaks for each cell cycles

Idents(multiome_seurat) <- multiome_seurat$Phase
levels(multiome_seurat)
#find differentially accessible peaks for each cell cycle
de_peaks_cc <- FindAllMarkers(multiome_seurat,
                assay = "ATAC",
                group.by = "Phase",
                only.pos = TRUE)

head(de_peaks_cc)

de_peaks_cc <- de_peaks_cc %>%
                  filter(pct.1 > pct.2 & p_val_adj <= 0.05)

de_peaks_cc <- de_peaks_cc %>%
                  distinct(gene, .keep_all = TRUE)
#count how many peaks in each phases
de_peaks_cc %>% count(cluster)



top_peaks <- de_peaks_cc %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC)

head(top_peaks)

unique(top_peaks$gene)

#plot heatmap of these peaks based on average values for each phase
peaks_chosen <- top_peaks$gene

peaks_chosen <- de_peaks_cc$gene

length(unique(de_peaks_cc$gene))

peaks_mean_list <- lapply(unique(multiome_seurat$Phase), function(phase){
  cell_phase <- which(multiome_seurat$Phase == phase)
  peaks_mean <- apply(multiome_seurat[["ATAC"]]$data[peaks_chosen, cell_phase], 1, mean)
  return (peaks_mean)
})

peaks_mean <- do.call(cbind, peaks_mean_list)
colnames(peaks_mean) <- unique(multiome_seurat$Phase)

#if use top DA peaks
annol_row <- data.frame(Phase = top_peaks$cluster)
rownames(annol_row) <- top_peaks$gene

#if use all DA peaks
annol_row <- data.frame(Phase = de_peaks_cc$cluster)
rownames(annol_row) <- de_peaks_cc$gene

head(annol_row)

anno_col <- data.frame(
  Phase = colnames(peaks_mean)  # or the grouping variable, like G1, G2M, S
)
rownames(anno_col) <- colnames(peaks_mean)



anno_colors <- list(
  cluster = c("G1" = "#E41A1C", "G2M" = "#377EB8", "S" = "#4DAF4A")
)


library(pheatmap)

pheatmap(peaks_mean,
         scale = "row",
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_names_col = FALSE,
         annotation_names_row = FALSE,
         annotation_col = anno_col,
         annotation_colors = anno_colors,
         annotation_row = annol_row)

#try to cluster all
cor_mat <- cor(multiome_seurat[["ATAC"]]$scale.data[top_peaks$gene, ], method = "spearman")
dist_mat <- as.dist(1 - cor_mat)  # convert correlation to distance
hc <- hclust(dist_mat, method = "ward.D2")
# cell_order <- hc$order

cell_clusters <- cutree(hc, 3)
table(cell_clusters)

table(multiome_seurat$Phase)

anno_col_cell <- data.frame(Phase = multiome_seurat$Phase)
rownames(anno_col_cell) <- colnames(multiome_seurat)

anno_colors <- list(
  cluster = c("G1" = "#E41A1C", "G2M" = "#377EB8", "S" = "#4DAF4A")
)

pheatmap(multiome_seurat[["ATAC"]]$scale.data[peaks_chosen, ],
         cluster_rows = TRUE,          # cluster genes
         # cluster_cols = hc,            # use your custom cell clustering
         cluster_cols = TRUE,
         scale = "row",               # already scaled
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_col = anno_col_cell,
         annotation_colors = anno_colors)



#use these peaks for PCA
#https://rpubs.com/crazyhottommy/pca-in-action

library(irlba)

DefaultAssay(multiome_seurat) <- "ATAC"

multiome_seurat <- ScaleData(multiome_seurat,
                             features = rownames(multiome_seurat))

multiome_seurat <- RunPCA(multiome_seurat,
                          assay = "ATAC",
                          features = de_peaks_cc$gene)
                          # features = top_peaks$gene)

p2 <- DimPlot(multiome_seurat,
              reduction = "pca",
              group.by= "Phase")

p2 + labs(title = "DA peaks from each phase")

#try to do manually
# https://divingintogeneticsandgenomics.com/post/pca-projection/

peaks_chosen <- top_peaks$gene
peaks_chosen <- de_peaks_cc$gene

# peak_mat <- multiome_seurat[["ATAC"]]$data[peaks_chosen, ]

peak_mat_scaled <- multiome_seurat[["ATAC"]]$scale.data[peaks_chosen, ]

peak_mat_scaled <- t(peak_mat_scaled)

#try log count

peak_mat <- multiome_seurat[["ATAC"]]$counts[peaks_chosen, ]

peak_mat <- apply(peak_mat, 2, function(col){log(col/sum(col)*1e4+1)})

peak_mat_scaled <- scale(t(peak_mat))

# Keep 100 PCs. The orginal seurat object kept 100 PCs

pca_peak_mat_scaled <- irlba(peak_mat_scaled, nv = min(100, length(peaks_chosen)-1))

#peak loadings (V matrix)
peak_loadings <- pca_peak_mat_scaled$v
# rownames(peak_loadings) <- top_peaks$gene
rownames(peak_loadings) <- peaks_chosen

colnames(peak_loadings) <- paste0("PC", 1:ncol(peak_loadings))

# Get PCA embeddings/cell embeddings (U matrix * D matrix)
cell_embeddings <- pca_peak_mat_scaled$u %*% diag(pca_peak_mat_scaled$d)  # Cell embeddings (10k cells in rows)
dim(cell_embeddings)

rownames(cell_embeddings) <- colnames(multiome_seurat)
colnames(cell_embeddings) <- colnames(peak_loadings)

head(cell_embeddings)

library(ggplot2)

cell_embeddings_df <- data.frame(PC1 = as.numeric(cell_embeddings[, 1]),
                                 PC2 = as.numeric(cell_embeddings[, 2]),
                                 Phase = multiome_seurat$Phase)

ggplot(cell_embeddings_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Phase)) +
  theme_classic(base_size = 14) +
  scale_x_continuous(breaks = pretty(cell_embeddings_df$PC1, n = 5)) +
  scale_y_continuous(breaks = pretty(cell_embeddings_df$PC2, n = 5))

#how about averaging then do PCA
#use


#try SCT counts

# rna_mat <- multiome_seurat[["SCT"]]$counts[c(g2m_genes, s_genes), ]

# rna_mat <- apply(rna_mat, 2, function(col){log(col/sum(col)*1e4+1)})

# rna_mat_scaled <- scale(t(rna_mat))

# common_genes <- intersect(rownames(multiome_seurat[["SCT"]]$scale.data), c(g2m_genes, s_genes))

# length(c(g2m_genes, s_genes))

rna_mat_scaled <- multiome_seurat[["SCT"]]$scale.data[c(g2m_genes, s_genes), ]

rna_mat_scaled <- t(rna_mat_scaled)
# Keep 20 PCs. The orginal seurat object kept 100 PCs

pca_rna_mat_scaled <- irlba(rna_mat_scaled, nv = 20)

#peak loadings (V matrix)
rna_loadings <- pca_rna_mat_scaled$v
# rownames(peak_loadings) <- top_peaks$gene
rownames(rna_loadings) <- c(g2m_genes, s_genes)

colnames(rna_loadings) <- paste0("PC", 1:20)

# Get PCA embeddings/cell embeddings (U matrix * D matrix)
cell_embeddings <- pca_rna_mat_scaled$u %*% diag(pca_rna_mat_scaled$d)  # Cell embeddings (10k cells in rows)
dim(cell_embeddings)

rownames(cell_embeddings) <- colnames(multiome_seurat[["SCT"]]$scale.data)
colnames(cell_embeddings) <- colnames(rna_loadings)

head(cell_embeddings)

library(ggplot2)

cell_embeddings_df <- data.frame(PC1 = as.numeric(cell_embeddings[, 1]),
                                 PC2 = as.numeric(cell_embeddings[, 2]),
                                 Phase = multiome_seurat$Phase)

ggplot(cell_embeddings_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Phase)) +
  theme_classic(base_size = 14) +
  scale_x_continuous(breaks = pretty(cell_embeddings_df$PC1, n = 5)) +
  scale_y_continuous(breaks = pretty(cell_embeddings_df$PC2, n = 5))

