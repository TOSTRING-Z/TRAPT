install.packages("svglite")
library(data.table)  # 用于高效处理大数据集
library(dplyr)       # 用于数据操作和转换
library(edgeR)       # 差异分析二号选手
library(limma)       # 差异分析三号选手
library(svglite)

############################ 表达矩阵和分组信息整理 ##############################

# 读取表达矩阵
setwd("E:\\SS硕士文件\\研一下学期_2024-2-27开始\\sc\\zgr改稿\\跑差异基因\\tcga3")
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

library(limma)
exprSet <- exp_brca

# 创建设计矩阵，指定组别信息
design <- model.matrix(~0 + factor(group))
colnames(design) = levels(factor(group))
rownames(design) = colnames(exprSet)

str(exprSet)

# 创建 DGEList 对象
exprSet[] <- lapply(exprSet, as.numeric)
dge <- DGEList(counts = exprSet, group = group)

#使用filterByExpr() 进行自动过滤，去除低表达基因
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

###找出前500差异基因
significant_genes <- re[re$P.Value < 0.05, ]
# 根据 logFC 排序
significant_genes <- significant_genes[order(abs(significant_genes$logFC), decreasing = TRUE), ]
top_significant_genes <- as.matrix(significant_genes[1:500,1])

colnames(top_significant_genes)[1]<-"Gene"	
write.table(top_significant_genes,"genes.txt",row.names=F,sep="\t",quote=F)

# 将差异表达基因结果保存到 Rdata 文件中
save(DEG_limma_voom, file = 'DEG_limma_voom.Rdata')

# 画火山图
load("./DEG_limma_voom.Rdata")

# 加change列,标记上下调基因
logFC = 1
P.Value = 0.05
k1 <- (DEG_limma_voom$P.Value < P.Value) & (DEG_limma_voom$logFC < -logFC)
k2 <- (DEG_limma_voom$P.Value < P.Value) & (DEG_limma_voom$logFC > logFC)
DEG_limma_voom <- mutate(DEG_limma_voom, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(DEG_limma_voom$change)


# 火山图
p <- ggplot(data = DEG_limma_voom, 
            aes(x = logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha = 0.4, size = 3.5, 
             aes(color = change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values = c("blue4", "grey", "red3"))+
  geom_vline(xintercept = c(-logFC, logFC), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(P.Value), lty = 4, col = "black", lwd = 0.8) +
  theme_bw()
p


ggsave(filename = "./volcano_plot_limma_voom.svg", plot = p, device = "svg", width = 6, height = 5)
dev.off()






