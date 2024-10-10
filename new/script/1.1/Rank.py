import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.preprocessing import minmax_scale
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--library", type=str, default="library")
parser.add_argument("--name", type=str, default=None)
parser.add_argument("--output_path", type=str, default=None)
parser.add_argument("--rank_path", type=str, default=None)
parser.add_argument("--type", type=str, default="ALL",help="ALL/TF/CR/TcoF/CR_TcoF")
parser.add_argument("--columns", type=str, default="TRAPT")
parser.add_argument("--source", type=str, default="KnockTF")
args = parser.parse_args()

name = args.name
title = args.type
output_path = args.output_path
rank_path = args.rank_path

color_palette = ["#3c5488", "#f39b7f", "#8491b4", "#91d1c2", "#fccde5"]

columns = args.columns.split(",")
data = pd.read_csv(f"{args.library}/TRs_info.txt", sep="\t")

ALL = set(data["tr_base"].drop_duplicates())

TF = set(
    data.loc[data["tr"].apply(lambda x: "Sample" in x), "tr_base"].drop_duplicates()
)
TcoF = set(
    data.loc[data["tr"].apply(lambda x: "TcoF" in x), "tr_base"].drop_duplicates()
)
CR = set(data.loc[data["tr"].apply(lambda x: "CR" in x), "tr_base"].drop_duplicates())


for method in columns:
    if method == "ChEA3":
        d = pd.read_csv(
            "other/chea3/down/ABL1@DataSet_03_001_down500$ARCHS4--Coexpression.tsv",
            sep="\t",
        )
        names = set(d["TF"])
    if method == "BART":
        d = pd.read_csv(
            "other/bart/down/ABL1@DataSet_03_001_down500_bart_results.txt", sep="\t"
        )
        names = set(d["TF"])
    if "TRAPT" in method:
        d = pd.read_csv(
            "output/KnockTFv1/AGO1@DataSet_02_95_down500/TR_detail.txt", sep="\t"
        )
        names = set(d["tr_base"])
    if method == "Lisa":
        d = pd.read_csv(
            "other/lisa/down/ABL1@DataSet_03_001_down500.txt.lisa.tsv", sep="\t"
        )
        names = set(d["factor"])
    if method == "i-cisTarget":
        d = pd.read_csv(
            "other/icistarget/down/ABL1@DataSet_03_001_down500/icistarget/statistics.tbl",
            sep="\t",
        )
        names = set(
            d["FeatureDescription"].apply(
                lambda x: x.split(" ")[-1] if x.startswith("ChIP") else x.split(" ")[0]
            )
        )
    inter_tcof = names.intersection(TcoF)
    inter_cr = names.intersection(CR)
    print(
        f"{method}, inter_tcof: {len(inter_tcof)}/{len(TcoF)}, inter_cr: {len(inter_cr)}/{len(CR)}"
    )

CR_TcoF = CR | TcoF

TR = {
    "ALL": ALL,
    "TF": TF,
    "CR": CR,
    "TcoF": TcoF,
    "CR_TcoF": CR_TcoF
}[args.type]

color = {
    "TF": "#e64b35",
    "TcoF": "#4dbbd5",
    "CR": "#00a087",
    "CR_TcoF": "#d476e2",
    "ALL": "#d776c2",
}

def get_rank_dict(file):
    software = pd.read_csv(file)
    software = software.loc[software["tr"].apply(lambda x: x in TR)]
    software_g = software.groupby("rank").agg({"tr": len}).reset_index()
    software_dict = dict(software_g.values)
    return software_dict


def get_auc_(x, y, s=1):
    direction = 1
    dx = np.diff(x)
    if np.any(dx < 0):
        direction = -1
    area = direction * np.trapz(y, x)
    auc_score = area / s
    return round(auc_score, 3)


rank_dict = {}
score_dict = {}
true_dict = {}
for method in columns:
    rank_dict[method] = get_rank_dict(f"{rank_path}/rank_{method}.csv")
    score_dict[method] = 0
    true_dict[method] = 0

# # # # # # # # # # # # # # # # # # # # # # # #
# # # # # #  MRR 柱状图  # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # #
print("MMR score:")
top = 1500
data = pd.DataFrame()
mrrs = []
for method in columns:
    add_data = pd.read_csv(f"{rank_path}/rank_%s.csv" % method)
    add_data = add_data.loc[add_data["tr"].apply(lambda x: x in TR)]
    add_data["Method"] = method
    add_data = add_data[add_data["rank"] <= top]
    mmr = np.round(np.mean(1 / pd.to_numeric(add_data["rank"])), 4)
    add_data["Mean Reciprocal Rank"] = mmr
    print(f"{method}: {mmr}")
    mrrs.append(mmr)
    data = pd.concat([data, add_data], ignore_index=True)

data["-Rank"] = -data["rank"]
data["Reciprocal Rank"] = 1 / data["rank"]

data = (
    data[["Method", "Mean Reciprocal Rank"]]
    .groupby("Method")
    .mean()
    .loc[columns]
    .reset_index()
)
g = sns.barplot(x="Mean Reciprocal Rank", y="Method", data=data, color=color[title])
for index, row in data.iterrows():
    print(row.name, row["Mean Reciprocal Rank"])
    g.text(
        x=row["Mean Reciprocal Rank"],
        y=row.name,
        s=round(row["Mean Reciprocal Rank"], 3),
        color="black",
        ha="left",
    )
sns.despine(bottom=False, left=False)
plt.title(f"Overall {title} recovery performance of differential genes in {args.source}")
plt.xlabel("Mean Reciprocal Rank")
plt.ylabel("Algorithms")
plt.savefig(f"{output_path}/rank_{name}@mmr_bar.svg")
plt.close()


# # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # mean rank 箱线图 # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # #

data = pd.DataFrame()
for method in columns:
    add_data = pd.read_csv(f"{rank_path}/rank_{method}.csv")
    add_data = add_data[add_data.tr.map(lambda x: x in TR)]
    # add_data = add_data[add_data['rank'] != add_data['rank'].max()]
    add_data["method"] = method
    # add_data = add_data[add_data['rank']<=200]
    add_data["rank"] = add_data["rank"] / add_data["rank"].max()
    data = pd.concat([data, add_data], ignore_index=True)

sn = sns.boxplot(
    x="method", y="rank", data=data, palette=color_palette, fliersize=1, width=0.8
)
# sn.set_ylim(top, 0 - 10)
plt.savefig(f"{output_path}/rank_{name}@boxplot.svg")
plt.close()