setwd("new/result/case2")

library(Seurat)
library(harmony)
library(dplyr)
library(tibble)
library(ggplot2)
library(readr)

seurat_standard_normalize_and_scale <- function(object){
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)
  object <- ScaleData(object, features = rownames(object))
  object <- RunPCA(object, features = VariableFeatures(object = object))
  return(object)
}

object.data <- read_csv("data.csv.gz") %>% column_to_rownames(var = "...1")
object.data[1:4,1:4]
dim(object.data)
set.seed(2024)

object.data[object.data < 0] <- 0
object <- CreateSeuratObject(counts = na.omit(object.data), min.cells = 3, min.features = 40)

str(object)

object[["orig.ident"]] <- unlist(lapply(strsplit(object@assays$RNA@data@Dimnames[[2]],"\\."),function(x)x[1]))
object[["CellType"]] <- unlist(lapply(strsplit(object@assays$RNA@data@Dimnames[[2]],"_"),function(x)x[1]))

object <- seurat_standard_normalize_and_scale(object)
object <- RunHarmony(object,group.by.vars = "orig.ident")
object <- seurat_standard_normalize_and_scale(object)
object <- RunUMAP(object, dims = 1:20, reduction = "pca")

library(RColorBrewer)
cell_type_cols <- c(brewer.pal(6, "Set1"), "#FF34B3", "#BC8F8F", "#20B2AA", "#00F5FF", "#FFA500", "#ADFF2F", "#FF6A6A", "#7FFFD4", "#AB82FF", "#90EE90", "#00CD00", "#008B8B", "#6495ED", "#FFC1C1", "#CD5C5C", "#8B008B", "#FF3030", "#7CFC00", "#000000", "#708090")
svg(paste0("cluster-umap.svg"))
DimPlot(object, reduction = "umap", group.by = "CellType", label = T, label.size = 3, cols= cell_type_cols, pt.size = 1,repel = T)
dev.off()

svg(paste0("cluster-pca.svg"))
DimPlot(object, reduction = "pca", group.by = "CellType", label = T, cols= cell_type_cols, pt.size = 1.5,repel = T)
dev.off()

unique(object[["CellType"]]$CellType)

diffs = list(
  c("HSC-MPPs","LMPs"),
  c("LMPs","MEMPs"),
  c("MEMPs","Megakaryocytes"),
  c("MEMPs","Erythroid cells"),
  c("MEMPs","Mast cells"),
  c("LMPs","Monocytes"),
  c("LMPs","GPs"),
  c("GPs","Granulocytes"),
  c("LMPs","pDCs"),
  c("LMPs","NK cells"),
  c("LMPs","Pro-B cells"),
  c("Pro-B cells","Pre-B cells"),
  c("Pre-B cells","Mature B cells")
)

output <- "output"
if(!dir.exists(output)) {
  dir.create(output)
}

for (item in diffs){
  print(item)
  h <- item[2]
  c <- item[1]
  markers <- FindMarkers(object, group.by="CellType", ident.1=h, ident.2=c, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25, max.cells.per.ident = 250)
  write.csv(markers,paste0("output/markers_top200-[",h,"-",c,"].csv"))
  write.table(data.frame(rownames(head(markers, 200))),paste0("output/markers_top200-[",h,"-",c,"].txt"),row.names=F,col.names=F, quote=F)
}