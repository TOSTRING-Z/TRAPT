import os
import pandas as pd
from rpy2.robjects import r
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
# 设置数据路径
path = "new/result/3.11/Fig. 6/GTEx/"
r.setwd(path)
# 读取处理好的基因表达数据
expr = pd.read_csv(f"{path}GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz", sep="\t",header=2,index_col=0).dropna()
expr = expr.set_index("Description")
phenotype = pd.read_csv(f"{path}GTEX_phenotype.gz", sep="\t")
sample = pd.DataFrame(expr.columns,columns=["Sample"])
sample = sample.merge(phenotype,on="Sample")
expr = expr.loc[:,sample["Sample"].values]
_primary_site = sample._primary_site.unique()
## 读取表达矩阵
with localconverter(ro.default_converter + pandas2ri.converter):
    exprSet_R = ro.conversion.py2rpy(expr)

r('''suppressMessages(library(limma))''')
################# 基因差异表达分析 #################
# 分组信息
for tissue in _primary_site:
    print(tissue)
    output = f"{path}Differentially_significant_genes/{tissue.replace(' ','_').strip()}"
    output_up = f"{path}Differentially_significant_genes_up/{tissue.replace(' ','_').strip()}"
    output_down = f"{path}Differentially_significant_genes_down/{tissue.replace(' ','_').strip()}"
    os.system(f"mkdir -p {output} {output_up} {output_down}")
    sample["Group"] = "B"
    sample.loc[sample._primary_site==tissue,"Group"] = "A"
    sample["Group_A"] = (sample["Group"]=="A").astype(int)
    sample["Group_B"] = (sample["Group"]=="B").astype(int)
    group_list = "-".join(list(sample["Group"].astype("str").unique())[::-1])
    ## 制作分组矩阵
    design = sample[["Group_A","Group_B"]].set_index(sample["Sample"])
    with localconverter(ro.default_converter + pandas2ri.converter):
        design_R = ro.conversion.py2rpy(design.rename(columns={"Group_A": "A","Group_B": "B"}))
    ## 制作差异比较矩阵
    contrast_matrix = r['makeContrasts'](group_list, levels=design_R)
    ## 使用limma包来进行差异分析
    fit = r["lmFit"](exprSet_R, design_R)
    fit2 = r["contrasts.fit"](fit, contrast_matrix)
    fit2 = r["eBayes"](fit2)
    tempOutput = r["topTable"](fit2, coef=1, number=float('inf'), lfc=0.5849625007211562)
    with localconverter(ro.default_converter + pandas2ri.converter):
        nrDEG = ro.conversion.rpy2py(tempOutput)
        nrDEG.dropna(axis=0, how='any', inplace=True)
        nrDEG.index = nrDEG.index.map(lambda i:expr.index[int(i)-1].split(".")[0])
    ## 筛选差异显著基因abs(log2FoldChange)>1.5,P.Value<0.05
    Differentially_significant_genes = nrDEG[nrDEG["adj.P.Val"] < 0.01].query('logFC>0.5849625007211562 or logFC<-0.5849625007211562')
    Differentially_significant_genes.head(500).index.to_frame().to_csv(output, sep="\t",header=False,index=False)
    up = nrDEG[nrDEG["adj.P.Val"] < 0.01].query('logFC>0.5849625007211562')
    up.head(500).index.to_frame().to_csv(output_up, sep="\t",header=False,index=False)
    down = nrDEG[nrDEG["adj.P.Val"] < 0.01].query('logFC<-0.5849625007211562')
    down.head(500).index.to_frame().to_csv(output_down, sep="\t",header=False,index=False)

