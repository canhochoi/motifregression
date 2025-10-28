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
#https://github.com/fangfang0906/Single_cell_multiome_palate/blob/main/Reproducibility/Figure_2.R
#https://www.nature.com/articles/s41467-024-45199-x#Sec2


gene_sig <- gene_markers_all %>%
  filter(p_val_adj <= 0.05 & avg_log2FC > 0.1) %>%
  pull(gene)

gene_sig <- unique(gene_sig)

DefaultAssay(multiome_seurat_noNA) <- "ATAC"

multiome_seurat_noNA <- RegionStats(multiome_seurat_noNA, genome = BSgenome.Hsapiens.UCSC.hg38)

#https://www.nature.com/articles/s41467-025-59306-z#code-availability
#fix p-value function which corresponds to 0 - 0.5 only

multiome_seurat_noNA <- LinkPeaks_modified(multiome_seurat_noNA,
                          peak.assay = "ATAC",
                          expression.assay = "SCT",
                          genes.use = gene_sig)

saveRDS(multiome_seurat_noNA, file = "/data1/soldatr/luan/projects/motifregression/result/CH/multiome_seurat_noNA_linkedpeaks.rds")

# note that the null distribution can generate bimodal data
# https://www.nature.com/articles/s41598-023-31040-w#Sec7

# alternatively can shuffle cell type
# https://www.sciencedirect.com/science/article/pii/S2666979X23001179#bib5

multiome_seurat_noNA <- readRDS("/data1/soldatr/luan/projects/motifregression/result/CH/multiome_seurat_noNA_linkedpeaks.rds")

DefaultAssay(multiome_seurat_noNA) <- "ATAC"
link <- Links(multiome_seurat_noNA)
link <- as.data.frame(link)
link$adj_pval <- p.adjust(link$pvalue, method = 'BH')

#identify which cluster each gene is uniquely expressed in
link$gene_cluster <- gene_markers_all[link$gene, 'cluster']

head(link)

# same as this
# all(gene_markers_all[link$gene, 'cluster'] == gene_markers_all[match(link$gene, gene_markers_all$gene), 'cluster'])
# link$gene_cluster <- gene_markers_all[match(link$gene, gene_markers_all$gene), 'cluster']


lvl <- levels(multiome_seurat_noNA)
lvl <- as.factor(sort(as.numeric(lvl)))

link <- arrange(link, factor(gene_cluster,levels=lvl))

link <- link[link$adj_pval<0.05 & link$score>0,]

#summarize number of peaks linked to each gene
t <- link %>%
  group_by(gene, gene_cluster) %>%
  summarise(n=n(), .groups = "drop")

table(t$gene_cluster)

t <- t %>%
  mutate(gene_cluster = factor(gene_cluster, levels = sort(as.numeric(levels(gene_cluster))))) %>%
  arrange(gene_cluster)

table(t$gene_cluster)


#plot number of linked peaks per marker genes in each cell cluster
ggplot(t, aes(x = gene_cluster, y = n)) +
  geom_violin(drop = FALSE) +
  geom_jitter() +
  theme_minimal() +
  labs(x = "Cell type", y = "Number of linked peaks per genes")

# Heatmap ####
# link_1 <- link[link$score>0.1,]
# link_2 <- link[link$score>0 & link$gene_cluster %in% c("Neuronal","Myocytes"),]
# link_combined <- rbind(link_1,link_2)

link <- arrange(link,factor(gene_cluster,levels=lvl))

# gene expression
gene.use <- unique(link$gene)
length(gene.use)

DefaultAssay(multiome_seurat_noNA) <- "SCT"
# multiome_seurat_noNA <- ScaleData(multiome_seurat_noNA, vars.to.regress = c("S.Score", "G2M.Score"), features=gene.use)

multiome_seurat_noNA <- ScaleData(multiome_seurat_noNA, features=gene.use)

rna <- multiome_seurat_noNA[["SCT"]]$scale.data[gene.use,]

#split cells into different clusters
asplit_cells <- split(colnames(multiome_seurat_noNA), multiome_seurat_noNA$projection)

#check
table(multiome_seurat_noNA$projection[asplit_cells[[2]]])

library(zoo)

n <- 10
#calculate average for every n cells in a cell type

means <- do.call(cbind, lapply(lvl, function(x){
  df <- rna[,asplit_cells[[x]]]
  t(apply(df,1,function(x){rollapply(x,n,mean,by=n)}))
}))

#do floor division
celltype_major <- unlist(lapply(lvl, function(x) rep(x, length(asplit_cells[[x]])%/%n)))

anno_col <- data.frame(celltype_major)
rownames(anno_col) <- colnames(means) <- paste0(seq(1:ncol(means)),colnames(means))
anno_col <- arrange(anno_col,factor(celltype_major,levels=lvl))
anno_col$celltype_major <- factor(anno_col$celltype_major,levels=lvl)

annotation_color <- list(celltype_major = multiome_seurat_noNA@misc$celltype_colors)

library(pheatmap)
pheatmap(means[link$gene, rownames(anno_col)],
         cluster_rows = F, cluster_cols = F, scale = "row",
         # breaks = seq(-2, 2, length = length(yellow2blue) + 1), col = yellow2blue,
         breaks = seq(-2, 2, length=101),
         annotation_col = anno_col,
         annotation_colors = annotation_color,
         show_colnames = F, show_rownames = F,
         annotation_names_col = F,
         main = 'Gene expression')

# chromatin accessibility
peak.use <- unique(link$peak)
length(peak.use)
DefaultAssay(multiome_seurat_noNA) <- "ATAC"

atac <- RunTFIDF(multiome_seurat_noNA[["ATAC"]]$counts[peak.use,])

asplit_cells <- split(colnames(multiome_seurat_noNA), multiome_seurat_noNA$projection)
means_peak <- do.call(cbind, lapply(lvl, function(x){
  df <- atac[,asplit_cells[[x]]]
  t(apply(df,1,function(x){rollapply(x,n,mean,by=n)}))
}))

celltype_major <- unlist(lapply(lvl, function(x) rep(x, length(asplit_cells[[x]])%/%n)))

anno_col2 <- data.frame(celltype_major)
rownames(anno_col2) <- colnames(means_peak) <- paste0(seq(1:ncol(means_peak)),colnames(means_peak))
anno_col2 <- arrange(anno_col2,factor(celltype_major,levels=lvl))

pheatmap(means_peak[link$peak,rownames(anno_col2)],
         cluster_rows = F, cluster_cols = F, scale = "row",
         # breaks = seq(-2, 2, length = 101), color = hcl.colors(100, "BluYl"),
         breaks = seq(-2, 2, length=101),
         annotation_col = anno_col2,
         annotation_colors = annotation_color,
         show_colnames = F, show_rownames = F,
         annotation_names_col = F,
         main = 'Chromatin accessibility')


# Correlation Plot ####
library(ggpubr)
#can try CRHBP, MLLT3 also for 1
gene <- "HLF"
gene <- "CRHBP"
gene <- "MLLT3"

tmp <- link[link$gene_cluster=="1",]
tmp <- na.omit(tmp[tmp$gene==gene,])
tmp <- arrange(tmp,-score)
peak1 <- tmp$peak[1]
region <- unlist(lapply(tmp[,'peak'],function(x){as.numeric(strsplit(x,"-",fixed=T)[[1]][c(2,3)])}))
peak2 <- paste(strsplit(tmp[1,'peak'],"-",fixed=T)[[1]][1],
               min(region,na.rm = T),
               max(region,na.rm=T),sep="-")

df <- data.frame(gene_expr=means[gene,],celltype=anno_col$celltype_major)
# p1 <- ggplot(df,aes(celltype,gene_expr,col=celltype))+geom_boxplot()+
#   theme_bw()+ggtitle(paste0(gene, " expression"))+
#   xlab("Cell type")+ylab("Relative expression")+
#   theme(axis.text.x = element_blank(),axis.ticks.x=element_blank())
p1 <- ggplot(df,aes(celltype,gene_expr,col=celltype))+geom_boxplot()+
  theme_bw() + ggtitle(paste0(gene, " expression"))+
  xlab("Cell type")+ylab("Relative expression")+
  theme(legend.position = "none",axis.ticks.x=element_blank())

df <- data.frame(atac=means_peak[peak1,],celltype=anno_col2$celltype_major)
p2 <- ggplot(df,aes(celltype,atac,col=celltype))+geom_boxplot()+
  theme_bw()+ggtitle(paste0(peak1, " accessibility"))+
  xlab("Cell type")+ylab("Chromatin accessibility")+
  # theme(axis.text.x = element_blank(),axis.ticks.x=element_blank())
  theme(legend.position = "none",axis.ticks.x=element_blank())

df <- data.frame(gene_expr=means[gene,],atac=means_peak[peak1,])
p3 <- ggplot(df,aes(gene_expr,atac))+geom_point()+geom_smooth(color="grey")+
  theme_bw()+ggtitle(paste0(gene, " & ", peak1))+
  xlab("gene expression")+ylab("chromatin accessibility")+
  stat_cor(method = "spearman", label.x = -0.5, label.y = 3.5)

p1 + p2 + p3


#plot coverage plot

Idents(multiome_seurat_noNA) <- factor(Idents(multiome_seurat_noNA), levels = as.factor(sort(as.numeric(levels(Idents(multiome_seurat_noNA))))))

#for HLF
CoveragePlot(
  object = multiome_seurat_noNA,
  region = gene,
  features = gene,
  expression.assay = "SCT",
  idents = lvl,
  extend.upstream = 250000,
  extend.downstream = 5000
)

#for CRHBP
CoveragePlot(
  object = multiome_seurat_noNA,
  region = gene,
  features = gene,
  expression.assay = "SCT",
  idents = lvl,
  extend.upstream = 15000,
  extend.downstream = 1000
)

#for MLLT3

CoveragePlot(
  object = multiome_seurat_noNA,
  region = gene,
  features = gene,
  expression.assay = "SCT",
  idents = lvl,
  extend.upstream = 1000,
  extend.downstream = 150000
)


#classify the peaks into exon, intron, intergenic, etc ....

