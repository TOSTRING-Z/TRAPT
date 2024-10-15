library(data.table)  # 用于高效处理大数据集
library(dplyr)       # 用于数据操作和转换
library(edgeR)
library(limma)
library(svglite)
library(ggplot2)


############################ 表达矩阵和分组信息整理 ##############################

# 读取表达矩阵
setwd("new/result/case1")
exp_brca <- fread("TCGA.BRCA.sampleMap_HiSeqV2.gz", header = T, sep = '\t', data.table = F)

# 查看表达矩阵
head(exp_brca)[1:5, 1:5]
# 去重
exp_brca <- avereps(exp_brca, exp_brca$sample) #可以对重复基因名取平均表达量再去重,可以用limma包中的avereps函数


# 把基因名转换为行名
exp_brca<-as.data.frame(exp_brca)
rownames(exp_brca) <- exp_brca$sample
exp_brca <- exp_brca[ , -1]

# 选取样本名14和15位置元素，因为它们可以代表样本类型
# 01 到 09 表示不同类型的肿瘤样本，10 到 19 表示不同类型的正常样本
# 01（原发性实体瘤）和 11（实体正常组织）是最常见的，06 则表示转移

gp <- substring(colnames(exp_brca), 14, 15)

# 分组
brca_tumor <- exp_brca[, as.numeric(gp) < 10]
brca_normal <- exp_brca[, as.numeric(gp) >= 10]

# 按顺序存储在一个矩阵中
exp_brca <- cbind(brca_tumor, brca_normal)

group <- c(rep('tumor', ncol(brca_tumor)), rep('normal', ncol(brca_normal)))
group <- factor(group, levels = c("normal", "tumor"))

########################### 使用 limma 进行差异分析 ############################

exprSet <- exp_brca

# 创建设计矩阵，指定组别信息
design <- model.matrix(~0 + factor(group))
colnames(design) = levels(factor(group))
rownames(design) = colnames(exprSet)

str(exprSet)

# 创建 DGEList 对象
exprSet[] <- lapply(exprSet, as.numeric)
dge <- DGEList(counts = exprSet, group = group)

# 使用filterByExpr() 进行自动过滤，去除低表达基因
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# 归一化，得到的归一化系数被用作文库大小的缩放系数
dge <- calcNormFactors(dge)

# 使用 voom 方法进行标准化
v <- voom(dge, design, plot = TRUE, normalize = "quantile")

# 使用线性模型进行拟合
fit <- lmFit(v, design)

con <- paste(rev(levels(group)), collapse = "-")

# 创建对比矩阵
cont.matrix <- makeContrasts(contrasts = c(con), levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# 获取差异表达基因结果
tempOutput <- topTable(fit2, coef = con, n = Inf)
DEG_limma_voom <- na.omit(tempOutput)


re<-cbind(rownames(DEG_limma_voom),DEG_limma_voom)
colnames(re)[1]<-"Gene_Symbol"

write.table(re,"Breast_cancer_limma.txt",row.names=F,sep="\t",quote=F)

########################### 找出前500差异基因 ############################

significant_genes <- re[re$P.Value < 0.05, ]
# 根据 logFC 排序
significant_genes <- significant_genes[order(abs(significant_genes$logFC), decreasing = TRUE), ]
top_significant_genes <- as.matrix(significant_genes[1:500,1])

colnames(top_significant_genes)[1]<-"Gene"	
write.table(top_significant_genes,"genes.txt",row.names=F,sep="\t",quote=F)

########################### 画火山图 ############################

res = read.table("./Breast_cancer_limma.txt",header=TRUE,row.names=1)
head(DEG_limma_voom)

library(tidyverse)
res <- res |> rownames_to_column("SYMBOL") |> arrange(desc(logFC))

# 给差异基因打标签，logFC > 1且 padj < 0.05认为是上调基因，logFC < -1且 padj < 0.05认为是下调基因
df <- res |>  
  mutate(significant = case_when(logFC > 1 & adj.P.Val < 0.05 ~ "Up",
                                 abs(logFC) < 1 | adj.P.Val > 0.05 ~ "None",
                                 logFC < -0.5 & adj.P.Val < 0.05 ~ "Down"))
df$significant |> as.factor()

head(df)

library(ggrepel)
library(ggfun)
library(grid)

# 指定显示上调前5名的基因，和下调前5名的基因
p <- ggplot(data = df) + 
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), 
                 color = logFC,
                #  size = -log10(adj.P.Val)
                 )) + 
  geom_text_repel(data =  df %>%
                    tidyr::drop_na() %>%
                    dplyr::filter(significant != "None") %>%
                    dplyr::filter(significant != "Down") %>%
                    dplyr::arrange(desc(-log10(adj.P.Val))) %>%
                    dplyr::slice(1:5) %>%
                    dplyr::filter(significant == "Up"),
                  aes(x = logFC, y = -log10(adj.P.Val), label = SYMBOL),
                  box.padding = 0.5,
                  nudge_x = 0.5,
                  nudge_y = 0.2,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 10,
                  direction = "y",
                  hjust = "left",
                  max.overlaps = 100
  )+
  geom_text_repel(data =  df %>%
                    tidyr::drop_na() %>%
                    dplyr::filter(significant != "None") %>%
                    dplyr::filter(significant != "Up") %>%
                    dplyr::arrange(desc(-log10(adj.P.Val))) %>%
                    dplyr::slice(1:5) %>%
                    dplyr::filter(significant == "Down"),
                  aes(x = logFC, y = -log10(adj.P.Val), label = SYMBOL),
                  box.padding = 0.5,
                  nudge_x = -0.2,
                  nudge_y = 0.2,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 10,
                  direction = "y", 
                  hjust = "left",
                  max.overlaps = 100
  ) + 
  scale_color_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                        values = seq(0, 1, 0.2)) +
  scale_fill_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                       values = seq(0, 1, 0.2)) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 4) + 
  scale_size(range = c(1,7)) + 
  ggtitle(label = "Volcano Plot") + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.background = element_roundrect(color = "#808080", linetype = 1),
        axis.text = element_text(size = 13, color = "#000000"),
        axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5)
  )

p

ggsave(filename = "./volcano_plot_limma_voom.svg", plot = p, device = "svg", width = 6, height = 5)
dev.off()