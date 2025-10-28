#for IRIS
multiome_seurat <- readRDS("~/luan/Juno/multiome_Rusland/multiome/data/processed/multiome_seurat.rds")
load(file = '~/luan/Juno/multiome_Rusland/multiome/data/processed/multiome_input.RData')
annotation.clusters <- readRDS(file = "~/luan/Juno/multiome_Rusland/multiome/data/processed/clusters_annotation.rds")
old.scrna <- new.env()
load(file = '~/luan/Juno/multiome_Rusland/RNA_analysis/data/processed/scrna.RData', envir = old.scrna)
save.image(file = "TF_gene_Rusland_multiome.RData")

#in Juno
multiome_seurat <- readRDS("/home/soldatr/CH/pipeline/multiome/data/processed/multiome_seurat.rds")
load(file = '/home/soldatr/CH/pipeline/multiome/data/processed/multiome_input.RData')
annotation.clusters <- readRDS(file = "/home/soldatr/CH/pipeline/multiome/data/processed/clusters_annotation.rds")
old.scrna <- new.env()
load(file = '/home/soldatr/CH/pipeline/RNA_analysis/data/processed/scrna.RData', envir = old.scrna)
save.image(file = "TF_gene_Rusland_multiome.RData")
