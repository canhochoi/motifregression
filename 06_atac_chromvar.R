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

library(tibble)
library(tidyr)

devtools::load_all("~/luan/projects/motifregression/")
devtools::document()


# has multiome_seurat_noNA
load(file='/data1/soldatr/luan/projects/motifregression/result/CH/multiome_metacell.RData')

multiome_seurat <- readRDS("~/luan/Juno/multiome_Rusland/multiome/data/processed/multiome_seurat.rds")


#from running 05_link_peaks_genes.R
multiome_seurat_noNA <- readRDS("/data1/soldatr/luan/projects/motifregression/result/CH/multiome_seurat_noNA_linkedpeaks.rds")

DefaultAssay(multiome_seurat_noNA) <- "ATAC"

#get list of motif PWM matrix for human
pwm <- getMatrixSet(x = JASPAR2020,
                    opts = list(species = 9606, all_versions = FALSE)
                    )
#add motif object to assay
#https://stuartlab.org/signac/articles/motif_vignette

multiome_seurat_noNA <- AddMotifs(object = multiome_seurat_noNA,
                                  genome = BSgenome.Hsapiens.UCSC.hg38,
                                  pfm = pwm)

#run chromvar

multiome_seurat_noNA <- RunChromVAR(
  object = multiome_seurat_noNA,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(multiome_seurat_noNA) <- "chromvar"
Idents(multiome_seurat_noNA) <- factor(Idents(multiome_seurat_noNA), levels = as.factor(sort(as.numeric(levels(multiome_seurat_noNA)))))
Idents(multiome_seurat_noNA)

#find cell type specific TF activity
TF_celltype <- FindAllMarkers(multiome_seurat_noNA,
                              only.pos = TRUE,
                              assay = "chromvar")
head(TF_celltype)

ggplot(TF_celltype, aes(x = avg_log2FC, y = -log(p_val_adj))) +
  geom_point(aes(color = p_val_adj < 0.05), alpha = 0.6) +
  facet_wrap(~ cluster, scales = "free_y") +
  theme_bw() +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot by Cluster",
       x = "log2 Fold Change",
       y = "-log10 Adjusted p-value",
       color = "FDR < 0.05")

#filter
TF_celltype <- TF_celltype %>%
                  filter(pct.1 > pct.2 & p_val_adj <= 0.05)
#find unique TFs
DA_TFs <- unique(TF_celltype$gene)

length(DA_TFs)
#follow this
#https://www.nature.com/articles/s41467-021-22368-w

average_chromvar <- Seurat::AverageExpression(multiome_seurat_noNA,
                                      assays = "chromvar",
                                      features = DA_TFs)
average_chromvar <- as.data.frame(average_chromvar)


#sort columns (i.e cell types) for maximum TF activity
#order command sorts indeces wiht values from small to high
average_chromvar <- average_chromvar[order(max.col(average_chromvar)), ]
colnames(average_chromvar) <- levels(Idents(multiome_seurat_noNA))

library(pheatmap)

pheatmap(average_chromvar, scale = "row",
         cluster_cols = FALSE, cluster_rows = FALSE,
         show_rownames = FALSE)

#check TF gene expression and chromvar activity across cell types

motifs <- multiome_seurat_noNA[["ATAC"]]@motifs@motif.names

motif_names <- data.frame(genes = unlist(motifs)) %>%
                rownames_to_column(var = "motif") %>%
                arrange(genes)

#make each motif corresponds to a gene
motif_names <- mutate(motif_names, gene = stringr::str_split(genes, pattern = "::")) %>%
  unnest(gene) %>%
  select(motif, gene) %>%
  distinct() %>%
  arrange() %>%
  filter(!stringr::str_detect(gene, pattern = "var.")) %>%
  as.data.frame()

#fix EWSR1-FLI1
motif_names$gene <- stringr::str_replace(motif_names$gene, pattern = "EWSR1-FLI1", replacement = "EWSR1")

# motif_gene_presence <- intersect(motif_names$gene, rownames(multiome_seurat_noNA[["SCT"]]$data))

# multiome_seurat_noNA[["SCT"]]$data[motif_gene_presence, ]

df_motif_list <- lapply(seq_along(motif_names$motif), function(motif_idx){
  motif <- motif_names$motif[motif_idx]
  print(motif)
  gene_motif <- motif_names[motif_idx, "gene"]
  print(gene_motif)
  DefaultAssay(multiome_seurat_noNA) <- "chromvar"
  motif_chromvar <- FindAllMarkers(multiome_seurat_noNA,
                                   logfc.threshold = 0,
                                   min.pct = 0,
                                   features = motif)
  DefaultAssay(multiome_seurat_noNA) <- "SCT"
  motif_expression <- FindAllMarkers(multiome_seurat_noNA,
                                     logfc.threshold = 0,
                                     min.pct = 0,
                                     features = gene_motif)
  motif_combined <- tryCatch(full_join(motif_chromvar, motif_expression, by = "cluster"),
           error = function(e) NULL)


  if (!is.null(motif_combined)){
    #remove NAs
    motif_combined <- motif_combined[complete.cases(motif_combined), ]

    if (nrow(motif_combined) > 0){
      df_motif <- data.frame(chromvar = motif_combined$avg_log2FC.x,
                             rna = motif_combined$avg_log2FC.y,
                             celltype = motif_combined$cluster,
                             motif = motif,
                             gene = gene_motif,
                             chromvar_pval = motif_combined$p_val_adj.x,
                             gene_pval = motif_combined$p_val_adj.y
      )
      return (df_motif)
    } else {
      return (NULL)
    }
  } else {
    return (NULL)
  }
})

df_motif <- bind_rows(df_motif_list)
df_motif <- dplyr::arrange(df_motif, gene, motif)

saveRDS(df_motif, file = "~/luan/projects/motifregression/result/CH/df_corr_chromvar_TF_exp.rds")

df_motif <- readRDS(file = "~/luan/projects/motifregression/result/CH/df_corr_chromvar_TF_exp.rds")

df_motif <- df_motif %>%
              filter(chromvar_pval <= 0.05 & gene_pval <= 0.05)
df_motif <- df_motif %>%
              mutate(df_motif, gene_motif = paste0(gene,"_",motif))

TF_motif_unique <- unique(df_motif$gene_motif)

#calculate pearson correlation between chromvar and motif gene expression across cell types
df_pearson <- lapply(TF_motif_unique, function(combo){
  require(ggrepel)

  df <- dplyr::filter(df_motif, gene_motif == combo)
  cor_res <- tryCatch(cor.test(df$chromvar, df$rna, method = "pearson", conf.level = 0.95),
                      error = function(e) NULL)
  df$cor <- tryCatch(signif(cor_res$estimate, 2), error = function(e) NULL)
  df$pval <- tryCatch(signif(cor_res$p.value, 2), error = function(e) NULL)
  df$max_exp <- round(max(abs(df$rna)), 2)
  df$max_chrom <- round(max(abs(df$chromvar)), 2)
  df$num_celltypes <- length(unique(df$celltype))
  return (df)
})





