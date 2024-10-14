import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--library", type=str, default="library")
parser.add_argument("--name", type=str, default=None)
parser.add_argument("--output_path", type=str, default=None)
parser.add_argument("--rank_path", type=str, default=None)
parser.add_argument("--type", type=str, default="ALL",help="ALL/TF/CR/TcoF/CR_TcoF")
parser.add_argument("--columns", type=str, default="TRAPT")
parser.add_argument("--col_names", type=str, default="TRAPT")
parser.add_argument("--source", type=str, default="KnockTF")
args = parser.parse_args()

name = args.name
title = args.type
output_path = args.output_path
rank_path = args.rank_path

color_palette = ["#3c5488", "#f39b7f", "#8491b4", "#91d1c2", "#fccde5"]

columns = args.columns.split(",")
col_names = args.col_names.split(",")
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
# # # # # # MRR 柱状图  # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # #
print("MMR score:")
data = pd.DataFrame()
for i,method in enumerate(columns):
    add_data = pd.read_csv(f"{rank_path}/rank_%s.csv" % method)
    add_data = add_data.loc[add_data["tr"].apply(lambda x: x in TR)]
    add_up = add_data[add_data["id"].str.contains("up")].copy()
    add_up["Type"] = "up"
    add_up["Method"] = col_names[i]
    mmr = np.round(np.mean(1 / pd.to_numeric(add_up["rank"])),4)
    add_up["Mean Reciprocal Rank"] = mmr
    add_down = add_data[add_data["id"].str.contains("down")].copy()
    add_down["Type"] = "down"
    add_down["Method"] = col_names[i]
    mmr = np.round(np.mean(1 / pd.to_numeric(add_down["rank"])),4)
    add_down["Mean Reciprocal Rank"] = mmr
    data = pd.concat([data, add_up, add_down], ignore_index=True)

data["-Rank"] = -data["rank"]
data["Reciprocal Rank"] = 1 / data["rank"]

data = (
    data[["Method", "Mean Reciprocal Rank", "Type"]]
    .groupby(["Method","Type"])
    .mean()
    .loc[col_names]
    .reset_index()
)
g = sns.barplot(x="Mean Reciprocal Rank", y="Method", hue="Type", data=data)
sns.despine(bottom=False, left=False)
plt.title(f"Overall {title} recovery performance of differential genes in {args.source}")
plt.xlabel("Mean Reciprocal Rank")
plt.ylabel("Algorithms")
plt.savefig(f"{output_path}/rank_{name}@mmr_bar.svg")
plt.close()
