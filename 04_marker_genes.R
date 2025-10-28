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

devtools::load_all("~/luan/projects/motifregression/")
devtools::document()


# has multiome_seurat_noNA
load(file='/data1/soldatr/luan/projects/motifregression/result/CH/multiome_metacell.RData')

multiome_seurat <- readRDS("~/luan/Juno/multiome_Rusland/multiome/data/processed/multiome_seurat.rds")

#follow this
#https://www.nature.com/articles/s41467-024-45199-x#Sec2
#https://github.com/fangfang0906/Single_cell_multiome_palate/blob/main/Reproducibility/Figure_1.R

DefaultAssay(multiome_seurat_noNA) <- "SCT"
Idents(multiome_seurat_noNA) <- multiome_seurat_noNA$projection

gene_markers_all <- FindAllMarkers(multiome_seurat_noNA,
                                   only.pos = TRUE,
                                   min.pct = 0.1,
                                   min.diff.pct = 0.1,
                                   logfc.threshold = 0.1)
head(gene_markers_all)

hist(gene_markers_all$avg_log2FC)
summary(gene_markers_all$avg_log2FC)

gene_markers_all <- gene_markers_all[gene_markers_all$pct.1 > gene_markers_all$pct.2 & gene_markers_all$p_val_adj<0.05, ]

head(gene_markers_all)

top_genes <- gene_markers_all %>%
        group_by(cluster) %>%
        top_n(5, avg_log2FC)

head(top_genes)

#check that some genes do not have gene coordinates

setdiff(top_genes$gene, rownames(multiome_seurat_noNA[["Promoter"]]))

DotPlot(multiome_seurat_noNA, features = unique(top_genes$gene)) +
  RotatedAxis() +
  theme(axis.title.y = element_blank(), axis.title.x=element_blank(),
        axis.text.x = element_text(size = 10, angle = 30, vjust = 0.8, hjust=0.5),
        axis.text.y = element_text(size = 10))+
  scale_color_viridis_c() + ggtitle("Gene expression")

DefaultAssay(multiome_seurat_noNA) <- "Promoter"

DotPlot(multiome_seurat_noNA, features = unique(top_genes$gene)) +
  RotatedAxis() +
  theme(axis.title.y = element_blank(), axis.title.x=element_blank(),
        axis.text.x = element_text(size = 10,angle = 45, vjust = 0.8, hjust=0.5),
        axis.text.y = element_text(size = 10))+
  scale_color_viridis_c() + ggtitle("Gene activity")

#core TFs
#https://www.nature.com/articles/ncb2709

TFs_core <- c("ERG", "TEL", "SCL", "RUNX1", "PU.1", "NFE2", "MITF",
              "MEIS1", "LYL1", "LMO2", "LDB1", "HHEX", "GFI1B", "GFI1",
              "GATA2", "GATA1", "FLI1", "ETO2")

TFs_shared <- intersect(TFs_core, rownames(multiome_seurat_noNA[["SCT"]]$data))



cor_mat <- cor(as.matrix(multiome_seurat_noNA[["SCT"]]$data[TFs_shared, ]), method = "spearman")
cor_mat[is.na(cor_mat)] <- 0
dist_mat <- as.dist((1 - cor_mat)/2)
hc <- hclust(dist_mat, method = "average")  # or "ward.D2", "complete", etc.
plot(hc)
#cut into 10 groups
clusters <- cutree(hc, k = 10)

library(pheatmap)

# Make sure column names in expression match names in clusters
expr_mat_tf <- multiome_seurat_noNA[["SCT"]]$data[TFs_shared, ]

# Create annotation for columns (cells)
anno_col <- data.frame(Cluster = factor(clusters))
rownames(anno_col) <- names(clusters)

# Plot heatmap
pheatmap(expr_mat_tf,
         scale = "row",  # normalize each gene
         clustering_method = "complete",  # optional
         annotation_col = anno_col,
         show_colnames = FALSE,
         main = "TF expression by cell cluster")
