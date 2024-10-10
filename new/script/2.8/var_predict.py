from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--rank_x", type=str, default="new/result/2.8/files/rank_H3K27ac.csv")
parser.add_argument("--rank_y", type=str, default="new/result/2.8/files/rank_ATAC.csv")
parser.add_argument("--output_path", type=str, default="new/result/2.8/files/rank_scatterplot.pdf")
parser.add_argument("--title", type=str, default="H3K27ac/ATAC rank")
args = parser.parse_args()

rank_x = pd.read_csv(args.rank_x)
rank_y = pd.read_csv(args.rank_y)
data = pd.merge(rank_x,rank_y,on=["id","tr"])
data["Type"] = data.id.str.extract("_(down|up)500$", expand=False)
data = data.groupby(["Type","tr"]).agg({
    'rank_x': np.mean,
    'rank_y': np.mean
}).reset_index()
data["x"] = 1/np.power(data.rank_x,0.5)
data["y"] = 1/np.power(data.rank_y,0.5)
data["score"] = data[["x","y"]].mean(axis=1)
data["top"] = data[["rank_x","rank_y"]].min(axis=1) <= 10
# 散点图
f, ax = plt.subplots(figsize=(6, 6))
sns.despine(f, ax, left=False, bottom=False)
sns.scatterplot(
    x="x",
    y="y",
    hue="Type",
    hue_order=["down", "up"],
    size="score",
    palette=["#4dbbd5","#e64b35"],
    data=data,
    ax=ax,
)
# 添加虚线，从坐标原点到右上角45度
ax.plot([0, 1], [0, 1], linestyle='--', color='black')

# 移除顶部和右侧的轴脊
sns.despine(f, ax, left=False, bottom=False)
quadrant = data[data["top"]].reset_index(drop=True)
for i in range(len(quadrant)):
    plt.annotate(
        quadrant.loc[i, "tr"],
        (
            quadrant.loc[i, "x"],
            quadrant.loc[i, "y"],
        ),
    )
plt.title(args.title)
plt.savefig(args.output_path)
plt.clf()