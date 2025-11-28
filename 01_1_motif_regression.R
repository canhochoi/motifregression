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


load(file='/data1/soldatr/luan/projects/motifregression/result/CH/multiome_metacell.RData')


library(randomForestSRC)
library(parallel)
library(scattermore)
library(data.table)
library(propagate)
library(tidyr)

library(rrpack)
library(Matrix)

plot_dir <- "~/luan/projects/motifregression/result/CH/plots/"

Idents(MC_seurat) <- factor(MC_seurat$cluster_cell, levels = c(1:10))


umap_p <- DimPlot(MC_seurat, label = FALSE, reduction = "wnn.umap",
        cols = MC_seurat@misc$celltype_colors)
ggsave(paste0(plot_dir, "umap_metacell_CH.pdf"), umap_p)

logcount_data <- apply(MC_seurat[["ATAC"]]$counts, 2, function(col_val){log(col_val/sum(col_val)*10000+1)})
logcount_data <- Matrix(logcount_data, sparse = TRUE)
peak_TF_mat <- MC_seurat[["ATAC"]]@motifs@data

#split to training and testing
nfolds <- 5
foldid_cell <- sample(rep(seq(nfolds), length = dim(logcount_data)[1]))
table(foldid_cell)
fold_idx <- 1
idx <- which(foldid_cell == fold_idx)


peak_TF_mat_train <- Matrix::Matrix(peak_TF_mat[-c(idx), ], sparse = TRUE)
peak_TF_mat_test <- Matrix::Matrix(peak_TF_mat[c(idx), ], sparse = TRUE)


Y_cell_train <- Matrix::Matrix(logcount_data[-c(idx), , drop = FALSE], sparse = TRUE)
Y_cell_test <- Matrix::Matrix(logcount_data[idx, , drop = FALSE], sparse = TRUE)


#do low rank ridge regression
rfit <- cv.rrr(Y = as.matrix(Y_cell_train), X = as.matrix(peak_TF_mat_train), nfold = 5, maxrank = 60)

rfit

#get y_predict from reduced rank regression
TF_activity_mat <- rfit$coef
y_predict_rr <- peak_TF_mat_test %*% TF_activity_mat
colnames(TF_activity_mat) <- colnames(MC_seurat[["ATAC"]]$data)
rownames(TF_activity_mat) <- names(MC_seurat[["ATAC"]]@motifs@motif.names)
regression_quality_list_rr <- check_regression(y_predict_rr, as.matrix(Y_cell_test))

df_rr <- data.table(prediction = as.vector(y_predict_rr),
                    observation = as.vector(Y_cell_test))

scatter_p_rr <- ggplot(df_rr, aes(x = prediction, y = observation)) +
  geom_scattermore(pointsize = 3, alpha = 0.1, pixels = c(1000, 1000), interpolate = TRUE) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Predicted peak accessibility", y = "Observed peak accessibility",
       title = paste0("Low rank, Pearson = ", round(regression_quality_list_rr$cor_val, 2), ", Rsq = ", round(regression_quality_list_rr$r_squared, 2))) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"))
scatter_p_rr
ggsave(paste0(plot_dir, "rrr_motif_regression.pdf"), scatter_p_rr)

# celltype_list <- unique(MC_seurat$cluster_cell)
celltype_list <- levels(Idents(MC_seurat))

#check correlation
cor_ct_rr <- mclapply(celltype_list, function(ct){
  idx <- which(MC_seurat$cluster_cell == ct)
  # cor_val <- cor(as.matrix(y_predict_rr[, idx, drop = FALSE]), as.matrix(Y_cell_test[, idx, drop = FALSE]))
  return (cor(as.vector(y_predict_rr[, idx, drop = FALSE]), as.vector(Y_cell_test[, idx, drop = FALSE])))
}, mc.cores = 4)

DefaultAssay(MC_seurat) <- "ATAC"

#get list of motif PWM matrix for human
pwm <- getMatrixSet(x = JASPAR2020,
                    opts = list(species = 9606, all_versions = FALSE)
)
#add motif object to assay
#https://stuartlab.org/signac/articles/motif_vignette

MC_seurat <- AddMotifs(object = MC_seurat,
                       genome = BSgenome.Hsapiens.UCSC.hg38,
                       pfm = pwm)

#run chromvar

MC_seurat <- RunChromVAR(
  object = MC_seurat,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# save(MC_seurat, file = "~/luan/projects/motifregression/result/CH/MC_seurat.rds")
saveRDS(MC_seurat, file = "~/luan/projects/motifregression/result/CH/MC_seurat.rds")

MC_seurat <- readRDS("~/luan/projects/motifregression/result/CH/MC_seurat.rds")

DefaultAssay(MC_seurat) <- "chromvar"
Idents(MC_seurat) <- factor(Idents(MC_seurat), levels = as.factor(sort(as.numeric(levels(MC_seurat)))))
Idents(MC_seurat)

#filter motifs to correspond to a gene
motifs <- MC_seurat[["ATAC"]]@motifs@motif.names

library(tibble)
library(tidyr)
library(stringr)
library(dplyr)

motif_names <- data.frame(genes = unlist(motifs)) %>%
  rownames_to_column(var = "motif") %>%
  arrange(genes)


#make each motif corresponds to a gene
motif_names <- mutate(motif_names, gene = stringr::str_split(genes, pattern = "::")) %>%
  unnest(gene) %>%
  dplyr::select(motif, gene) %>%
  dplyr::distinct() %>%
  arrange(gene) %>%
  dplyr::filter(!stringr::str_detect(gene, pattern = "var.")) %>%
  as.data.frame()

#fix EWSR1-FLI1
motif_names$gene <- stringr::str_replace(motif_names$gene, pattern = "EWSR1-FLI1", replacement = "EWSR1")

#remove genes without presence in data
motif_names <- motif_names %>%
            filter(gene %in% rownames(MC_seurat[["SCT"]]$data))

#some genes have multiple motifs
length(intersect(unique(motif_names$gene), rownames(MC_seurat[["SCT"]]$data)))

#plot correlation with TF gene expression----
#compare with chromvar
#use match to have one-to-one correspondence

#find which motifs in motifs name match the filtered motifs
matching_rows <- names(motifs)[match(motif_names$motif, names(motifs))]

all(matching_rows == motif_names$motif)


#check TF activity correlation with TF gene expression
cor_val <- cor(as.matrix(t(TF_activity_mat[matching_rows, , drop = FALSE])), as.matrix(t(MC_seurat[['SCT']]$data[motif_names$gene, , drop = FALSE])))
# cor_val_TF <- sapply(seq(dim(MC_seurat[['SCT']]$data)[2]), function(i){cor_val[i,i]})

cor_val_TF <- diag(cor_val)[is.na(diag(cor_val)) == FALSE]

#for chromvar
TF_activity_chromvar <- MC_seurat[["chromvar"]]$data
cor_val_chromvar <- cor(as.matrix(t(TF_activity_chromvar[matching_rows, , drop = FALSE])), as.matrix(t(MC_seurat[['SCT']]$data[motif_names$gene, , drop = FALSE])))

cor_val_chromvar <- diag(cor_val_chromvar)[is.na(diag(cor_val_chromvar)) == FALSE]

#plot
cor_val_TF_df <- data.frame(cor_val_TF)

histogram_plot_TF <- ggplot(cor_val_TF_df, aes(x = cor_val_TF)) +
  geom_histogram(binwidth = 0.05, color = "black", fill = "#F8766D", alpha = 0.5) +
  labs(title = paste0("Motif TF"), x = "TF activity - gene correlation Value", y = "Frequency") +
  theme_minimal(base_size = 22)

cor_val_chromvar_df <- data.frame(cor_val_chromvar)

histogram_plot_chromvar <- ggplot(cor_val_chromvar_df, aes(x = cor_val_chromvar)) +
  geom_histogram(binwidth = 0.05, color = "black", fill = "skyblue", alpha = 0.5) +
  labs(title = paste0("Chromvar"), x = "TF activity - gene correlation Value", y = "Frequency") +
  theme_minimal(base_size = 22)

histogram_plot_TF
histogram_plot_chromvar

cor_val_TF_df$method <- "Motif"
cor_val_chromvar_df$method <- "ChromVAR"

combined_df <- rbind(
  data.frame(cor_val = cor_val_TF_df$cor_val_TF, method = "TF"),
  data.frame(cor_val = cor_val_chromvar_df$cor_val_chromvar, method = "ChromVAR")
)

hist_p <- ggplot(combined_df, aes(x = cor_val, fill = method)) +
  geom_histogram(binwidth = 0.05, position = "identity", alpha = 0.25, color = "black") +
  scale_fill_manual(values = c("TF" = "#F8766D", "ChromVAR" = "skyblue")) +
  labs(title = "TF activity - gene correlation", x = "Correlation Value", y = "Frequency") +
  theme_minimal(base_size = 22)

ggsave(paste0(plot_dir, "histogram_chromvar_motif.pdf"),
       hist_p,
       width = 161,
       height = 161,
       units = "mm")

#check correlation and its significance
cor_list <- lapply(seq_along(matching_rows), function(idx_motif){
  cor_result <- cor.test(as.vector(TF_activity_mat[matching_rows[idx_motif], , drop = FALSE]), as.vector(MC_seurat[['SCT']]$data[motif_names$gene[idx_motif], , drop = FALSE]))
  return (list(p_val = cor_result$p.value, cor_val = cor_result$estimate))
})


cor_df <- as.data.frame(do.call(rbind, cor_list))
cor_df$cor_val <- as.numeric(cor_df$cor_val)
cor_df$p_val <- as.numeric(cor_df$p_val)
cor_df$motif_id <- matching_rows
# rownames(cor_df) <- matching_rows


cor_list_chromvar <- lapply(seq_along(matching_rows), function(idx_motif){
  cor_result <- cor.test(as.vector(TF_activity_chromvar[matching_rows[idx_motif], , drop = FALSE]), as.vector(MC_seurat[['SCT']]$data[motif_names$gene[idx_motif], , drop = FALSE]))
  return (list(p_val = cor_result$p.value, cor_val = cor_result$estimate))
})

cor_df_chromvar <- as.data.frame(do.call(rbind, cor_list_chromvar))
cor_df_chromvar$cor_val <- as.numeric(cor_df_chromvar$cor_val)
cor_df_chromvar$p_val <- as.numeric(cor_df_chromvar$p_val)
cor_df_chromvar$motif_id <- matching_rows
# rownames(cor_df_chromvar) <- single_motif_common

motif_sig_df <- cor_df %>%
  dplyr::filter(p_val <= 0.05)

motif_sig_chromvar_df <- cor_df_chromvar %>%
  dplyr::filter(p_val <= 0.05)

# common_motif_sig <- intersect(rownames(motif_sig_chromvar_df), rownames(motif_sig_df))

common_motif_sig <- intersect(motif_sig_chromvar_df$motif_id, motif_sig_df$motif_id)


#plot common motif correlation between chromvar and our method

motif_sig_df_common <- motif_sig_df %>%
            filter(motif_id %in% common_motif_sig) %>%
            slice(match(common_motif_sig, motif_id))

motif_sig_chromvar_df_common <- motif_sig_chromvar_df %>%
  filter(motif_id %in% common_motif_sig) %>%
  slice(match(common_motif_sig, motif_id))


cor_motif_chromvar <- data.frame(motif = motif_sig_df_common$cor_val,
                                 chromvar = motif_sig_chromvar_df_common$cor_val)

rownames(cor_motif_chromvar) <- common_motif_sig

cor.test(cor_motif_chromvar$motif, cor_motif_chromvar$chromvar)

cor_sc_p <- ggplot(cor_motif_chromvar, aes(x = chromvar, y = motif)) +
  geom_point(alpha = 0.5, color = "black", size = 5) +
  labs(title = "TF activity and TF expression, correlation", x = "Chromvar", y = "Motif Regression") +
  theme_minimal(base_size = 30) +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
  ) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")
cor_sc_p

ggsave(paste0(plot_dir, "TF_activity_gene_motif_vs_chromvar.pdf"),
       cor_sc_p,
       width = 161,
       height = 161,
       units = "mm")

#check correlation between chromvar and motif regression

library(scMEGA)

MC_seurat <- AddTrajectory(object = MC_seurat,
                      trajectory = c("1", "6", "7", "8"),
                      group.by = "cluster_cell",
                      reduction = "wnn.umap",
                      dims = 1:2,
                      use.all = FALSE)

MC_seurat_traj <- MC_seurat[, !is.na(MC_seurat$Trajectory)]

#continuousSet is from names(ArchRPalettes)

traj_p <- TrajectoryPlot(object = MC_seurat_traj,
               reduction = "wnn.umap",
               # continuousSet = "blueYellow",
               continuousSet = "grove",
               size = 1,
               addArrow = FALSE)

traj_p

ggsave(paste0(plot_dir, "traj_1_678.pdf"),
       traj_p,
       width = 161,
       height = 161,
       units = "mm")

sel.tfs <- SelectTFs(object = MC_seurat_traj,
                     return.heatmap = TRUE,
                     cor.cutoff = 0.4)

df.cor <- sel.tfs$tfs
ht <- sel.tfs$heatmap

pdf(file = paste0(plot_dir, "TF_activity_expression.pdf"),
    width = 10,
    height = 20)
draw(ht)  # or heatmapt_TF_p if that's your heatmap object
dev.off()


sel.genes <- SelectGenes(object = MC_seurat_traj,
                         labelTop1 = 0,
                         labelTop2 = 0)

df.p2g <- sel.genes$p2g
ht <- sel.genes$heatmap

draw(ht)


tf.gene.cor <- GetTFGeneCorrelation(object = MC_seurat_traj,
                                    tf.use = df.cor$tfs,
                                    gene.use = unique(df.p2g$gene),
                                    tf.assay = "chromvar",
                                    gene.assay = "SCT",
                                    trajectory.name = "Trajectory")

mat.cor <- tf.gene.cor %>%
  as.data.frame() %>%
  select(c(tf, gene, correlation)) %>%
  tidyr::pivot_wider(names_from = tf, values_from = correlation) %>%
  textshape::column_to_rownames("gene")


ht <- GRNHeatmap(tf.gene.cor,
                 tf.timepoint = df.cor$time_point)

draw(ht)

#associate gene to TFs
motif.matching <- MC_seurat[["ATAC"]]@motifs@data
colnames(motif.matching) <- MC_seurat_traj[["ATAC"]]@motifs@motif.names
motif.matching <- motif.matching[unique(df.p2g$peak), unique(tf.gene.cor$tf)]


df.grn <- GetGRN(motif.matching = motif.matching,
                 df.cor = tf.gene.cor,
                 df.p2g = df.p2g)

# define colors for nodes representing TFs (i.e., regulators)
df.cor <- df.cor[order(df.cor$time_point), ]
tfs.timepoint <- df.cor$time_point
names(tfs.timepoint) <- df.cor$tfs

# plot the graph, here we can highlight some genes
df.grn2 <- df.grn %>%
  subset(correlation > 0.5) %>%
  dplyr::select(c(tf, gene, correlation)) %>%
  dplyr::rename(weights = correlation)

library(igraph)
library(ggraph)

p <- GRNPlot(df.grn2,
             tfs.timepoint = tfs.timepoint,
             show.tf.labels = TRUE,
             seed = 42,
             plot.importance = FALSE,
             min.importance = 2,
             remove.isolated = FALSE)

print(p)

#visualize individual TFs
MC_seurat_traj <- AddTargetAssay(object = MC_seurat_traj, df.grn = df.grn2)

p1 <- PseudotimePlot(object = MC_seurat_traj, tf.use = "GATA3")

p1

p2 <- PseudotimePlot(object = MC_seurat_traj, tf.use = "FOXA3")

p2

p3 <- PseudotimePlot(object = MC_seurat_traj, tf.use = "CEBPD")

ggsave(paste0(plot_dir, "CEBPD.pdf"),
       p3,
       width = 161,
       height = 161,
       units = "mm")


#do MRF - too large so causing OOM
cores_max <- detectCores()
print(cores_max)

options(rf.cores = 12, mc.cores = 12)

y_df <- data.frame(Y_cell_train[, , drop = FALSE])
x_df <- data.frame(peak_TF_mat_train[, , drop = FALSE]*1)

ntree_val <- 5
obj <- rfsrc(get.mv.formula(colnames(y_df)),
             data.frame(y_df, x_df),
             ntree = ntree_val,
             importance=TRUE, nsplit = 10, splitrule = "mahalanobis")
y_predict_obj <- predict(obj, data.frame(peak_TF_mat_test[, , drop = FALSE]))
y_predict_mrf <- as.matrix(get.mv.predicted(y_predict_obj))
regression_quality_list <- check_regression(y_predict = y_predict_mrf, y_observed = as.matrix(Y_cell_test[, , drop = FALSE]))

