import pandas as pd
import numpy as np

# Chen, CH., Zheng, R., Tokheim, C. et al. Determinants of transcription factor regulatory range. Nat Commun 11, 2472 (2020). https://doi.org/10.1038/s41467-020-16106-x
"""shell
wget -c https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-16106-x/MediaObjects/41467_2020_16106_MOESM4_ESM.csv -O new/result/2.5/41467_2020_16106_MOESM4_ESM.csv
"""

tr_type = pd.read_csv("new/result/2.5/41467_2020_16106_MOESM4_ESM.csv",index_col=0)[["assignment"]]
tr_type.index.name = "TR"


# TRAPT
data = pd.read_csv("new/result/2.5/decay_rate.txt",sep="\t")

data_g = data.groupby("TR").agg({
    "decay rate": np.mean,
    "sample": len
})
# data_g = data_g.query("sample > 3")
mean = data_g["decay rate"].mean()
alpha = lambda m: np.log(2/m - 1)/90000
l_m = 0.22
h_m = 0.22
data_g["group"] = data_g["decay rate"].map(lambda s:"short-range TF" if s > alpha(l_m) else ("long-range TF" if s < alpha(h_m) else "other"))

data_g_inter = data_g.reset_index().merge(tr_type.reset_index(),on="TR")
same = data_g_inter.query("group==assignment").reset_index(drop=True)

g_long = data_g_inter.query("group=='long-range TF'").reset_index(drop=True)
g_short = data_g_inter.query("group=='short-range TF'").reset_index(drop=True)

same_long = same.query("assignment=='long-range TF'").reset_index(drop=True)
same_short = same.query("assignment=='short-range TF'").reset_index(drop=True)

pri_long = data_g_inter.query("assignment=='long-range TF'").reset_index(drop=True)
pri_short = data_g_inter.query("assignment=='short-range TF'").reset_index(drop=True)
print(f"g_long: {len(g_long)}, g_short: {len(g_short)}\n\
same_long: {len(same_long)}, same_short: {len(same_short)}\n\
pri_long: {len(pri_long)}, pri_short: {len(pri_short)}")

data_g.to_csv("new/result/2.5/decay_rate[aggregation].txt",sep="\t")

# 韦恩图
from matplotlib_venn import venn2, venn2_circles
import venn
import matplotlib.pyplot as plt
from scipy.stats import hypergeom

def draw_venn(type_,name,g_,pri_):
    # 定义集合
    set1 = set(g_.TR)
    set2 = set(pri_.TR)
    # 超几何检验的参数
    x = len(set1 & set2)
    M = len(data_g)
    n = len(data_g_inter)
    N = len(set1)
    # 执行超几何检验 (in set - 1, background, in background, set)
    p_value = hypergeom.sf(x - 1, M, n, N)
    print(x - 1, M, n, N)
    print(f"交集大小：{x}")
    print(f"超几何检验的P值：{p_value}")

    mycolor = [
        [0.10588235294117647, 0.6196078431372549, 0.4666666666666667, 0.6],
        [0.9058823529411765, 0.1607843137254902, 0.5411764705882353, 0.6],
    ]
    labels = venn.get_labels([set1, set2], fill=["number", "logic", "percent"])

    def labels_formatter(c):
        for label in labels.values():
            if f": {c} " in label:
                return label.split(": ")[1]

    plt.figure(figsize=(10, 8))
    ax = venn2(
        [set1, set2],
        set_labels=("TRAPT", "Chen. et al."),
        set_colors=mycolor,
        subset_label_formatter=labels_formatter,
    )
    c = venn2_circles(subsets=[set1, set2], linestyle="dashed")
    # 在图像的右下角添加文本，并设置字体大小为 12
    plt.text(1, 0, 'p = {:.2e}'.format(p_value), horizontalalignment='right', verticalalignment='bottom', transform=plt.gca().transAxes, fontsize=12)
    plt.title("# of %s-range TF" % type_)
    plt.savefig(f"new/result/2.5/{name}.svg")

draw_venn("long","long",g_long,pri_long)
draw_venn("short","short",g_short,pri_short)

