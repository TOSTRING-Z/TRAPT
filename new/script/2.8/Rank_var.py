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


# # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # mean rank 箱线图 # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # #

data = pd.DataFrame()
for method in columns:
    add_data = pd.read_csv(f"{rank_path}/rank_{method}.csv")
    add_up = add_data[add_data["id"].str.contains("up")]
    add_up = add_up[add_up.tr.map(lambda x: x in TR)]
    add_up["type"] = "up"
    add_up["method"] = method
    add_up["rank"] = add_up["rank"] / add_up["rank"].max()
    add_down = add_data[add_data["id"].str.contains("down")]
    add_down = add_down[add_down.tr.map(lambda x: x in TR)]
    add_down["type"] = "down"
    add_down["method"] = method
    add_down["rank"] = add_down["rank"] / add_down["rank"].max()
    data = pd.concat([data, add_up, add_down], ignore_index=True)

sn = sns.boxplot(
    x="method", y="rank",hue="type", data=data, palette=color_palette, fliersize=1, width=0.8
)
plt.title(f"Ranking of target TRs using knockdown/knockout datasets")
plt.savefig(f"{output_path}/rank_{name}@boxplot.svg")
plt.close()
