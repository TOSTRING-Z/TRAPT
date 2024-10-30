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
parser.add_argument("--col_names", type=str, default=None)
parser.add_argument("--source", type=str, default="KnockTF")
args = parser.parse_args()

name = args.name
title = args.type
output_path = args.output_path
rank_path = args.rank_path

color_palette = ["#d4738b", "#e8ccbb", "#edf4f7", "#9fdadb", "#648d9c"]

columns = args.columns.split(",")
if args.col_names == None:
    col_names = columns
else:
    col_names = args.col_names.split(",")
col_dict = dict(zip(columns,col_names))
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
    elif method == "BART":
        d = pd.read_csv(
            "other/bart/down/ABL1@DataSet_03_001_down500_bart_results.txt", sep="\t"
        )
        names = set(d["TF"])
    elif "TRAPT" in method:
        d = pd.read_csv(
            "output/KnockTFv1/AGO1@DataSet_02_95_down500/TR_detail.txt", sep="\t"
        )
        names = set(d["tr_base"])
    elif method == "Lisa":
        d = pd.read_csv(
            "other/lisa/down/ABL1@DataSet_03_001_down500.txt.lisa.tsv", sep="\t"
        )
        names = set(d["factor"])
    elif method == "i-cisTarget":
        d = pd.read_csv(
            "other/icistarget/down/ABL1@DataSet_03_001_down500/icistarget/statistics.tbl",
            sep="\t",
        )
        names = set(
            d["FeatureDescription"].apply(
                lambda x: x.split(" ")[-1] if x.startswith("ChIP") else x.split(" ")[0]
            )
        )
    else:
        break
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

# columns = {
#     "TF": columns,
#     "TcoF": list(filter(lambda c:c not in ["BART", "ChEA3"],columns)),
#     "CR": list(filter(lambda c:c not in ["BART", "ChEA3"],columns)),
#     "CR_TcoF": list(filter(lambda c:c not in ["BART", "ChEA3"],columns)),
#     "ALL": columns,
# }[title]


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
# # # # # # # # #  分组柱状图  # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # #
top = 10
data = []
scores = [0] * len(columns)
for i in range(1, top + 1):
    for method in columns:
        if i in rank_dict[method]:
            score_dict[method] += rank_dict[method].get(i)
    data.append(list(score_dict.values()))

data = pd.DataFrame(data, columns=columns)

w = 1 / 2 / data.shape[1]
for i in range(data.shape[1]):
    plt.bar(
        data.index - 0.25 + i * w + 0.5 * w,
        data.iloc[:, i],
        width=w,
        label=data.columns[i],
        color=color_palette[i],
        edgecolor="k",
        zorder=2,
    )

plt.xlabel("Top N", fontsize=11)
plt.ylabel("Number of sets with correct predicted in top N", fontsize=11)
plt.legend(frameon=False)
plt.grid(ls="--", alpha=0.8)
plt.xticks(data.index, data.index + 1, fontsize=10)
plt.tick_params(axis="x", length=0)

plt.tight_layout()
plt.savefig(f"{output_path}/rank_{name}@groupbar.svg")
plt.close()

# # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # #  rank auc 折线图   # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # #
print("Rank AUC:")
top = 10
data = []
score_dict = {}
for method in columns:
    score_dict[method] = 0

for i in range(1, top + 1):
    for method in columns:
        if i in rank_dict[method]:
            score_dict[method] += rank_dict[method].get(i)
    data.append(list(score_dict.values()))

data = pd.DataFrame(data, columns=columns)

x = np.array(range(top))
s = top * data.iloc[:top].values.max()
for i in range(len(columns)):
    col = columns[i]
    col_name = col_dict.get(col)
    auc = get_auc_(x, data[col].values, s)
    print(f"{col}: {auc}")
    plt.fill_between(x, data[col], 0, alpha=0.1, color=color_palette[i])
    plt.plot(
        x,
        data[col],
        "--",
        label=f"{col_name}: {auc}",
        color=color_palette[i],
        linewidth=1.5,
    )
    plt.plot(
        x,
        data[col],
        "o",
        color=color_palette[i],
        markersize=6,
        linewidth=2,
        markeredgecolor="black",
        markeredgewidth=0.5,
    )

plt.title(f"Top N {title} recovery number of differential genes in {args.source}")
plt.ylabel("Number of sets with correct predicted in top N")
plt.xlabel("Top N")
plt.legend(loc="lower right")
plt.savefig(f"{output_path}/rank_{name}@line.svg")
plt.close()

# # # # # # # # # # # # # # # # # # # # # # # #
# # # # # #  boxplot and MRR 柱状图  # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # #
print("MMR score:")
top = 1500
data = pd.DataFrame()
mrrs = []
for method in columns:
    col_name = col_dict.get(method)
    add_data = pd.read_csv(f"{rank_path}/rank_%s.csv" % method)
    add_data = add_data.loc[add_data["tr"].apply(lambda x: x in TR)]
    add_data["Method"] = col_name
    add_data = add_data[add_data["rank"] <= top]
    mmr = np.round(np.mean(1 / pd.to_numeric(add_data["rank"])),4)
    add_data["Mean Reciprocal Rank"] = mmr
    print(f"{col_name}: {mmr}")
    mrrs.append(mmr)
    data = pd.concat([data, add_data], ignore_index=True)

data["-Rank"] = -data["rank"]
data["Reciprocal Rank"] = 1 / data["rank"]

sn = sns.boxplot(x="Method", y="-Rank", data=data, palette=color_palette, fliersize=1)

plt.savefig(f"{output_path}/rank_{name}@boxplot.svg")
plt.close()

data = (
    data[["Method", "Mean Reciprocal Rank"]]
    .groupby("Method")
    .mean()
    .loc[col_names]
    .reset_index()
)

g = sns.barplot(x="Mean Reciprocal Rank", y="Method", data=data, palette=color_palette, edgecolor='#3b3b3b')
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
# # # # # # # # # #  柱状图  # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # #
print("Correct number:")
top = 10
true_dict = {}
for method in columns:
    true_dict[method] = 0

for i in range(1, top + 1):
    for method in columns:
        true_dict[method] += rank_dict[method].get(i) if i in rank_dict[method] else 0

data_bar = pd.DataFrame([true_dict]).T
data_bar = data_bar.reset_index()
data_bar.columns = ["Method", "Correct number"]
data_bar["Method"] = data_bar["Method"].map(col_dict.get)

g = sns.barplot(y="Method", x="Correct number", data=data_bar, palette=color_palette, edgecolor='#3b3b3b')
sns.despine(bottom=False, left=False)
for index, row in data_bar.iterrows():
    print(row.name, row["Correct number"])
    g.text(
        x=row["Correct number"],
        y=row.name,
        s=row["Correct number"],
        color="black",
        ha="left",
    )
sns.despine(bottom=False, left=False)  # 设置是否显示边界线
plt.title(f"Top 10 {title} recovery number of differential genes in {args.source}")
plt.ylabel("Algorithms")
plt.xlabel("Number of sets with correct predicted in top 10")
plt.savefig(f"{output_path}/rank_{name}@bar.svg")
plt.close()

# # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # mean rank 箱线图 # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # #

data = pd.DataFrame()
for method in columns:
    col_name = col_dict.get(method)
    add_data = pd.read_csv(f"{rank_path}/rank_{method}.csv")
    add_data = add_data[add_data.tr.map(lambda x: x in TR)]
    # add_data = add_data[add_data['rank'] != add_data['rank'].max()]
    add_data["method"] = col_name
    # add_data = add_data[add_data['rank']<=200]
    add_data["rank"] = add_data["rank"] / add_data["rank"].max()
    data = pd.concat([data, add_data], ignore_index=True)

sn = sns.boxplot(
    x="method", y="rank", data=data, palette=color_palette, fliersize=1, width=0.8
)
# sn.set_ylim(top, 0 - 10)
plt.savefig(f"{output_path}/rank_{name}@boxplot.svg")
plt.close()
