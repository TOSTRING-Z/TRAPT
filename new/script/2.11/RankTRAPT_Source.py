import os
import re
from glob import glob
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--output_dir", type=str, default=None)
parser.add_argument("--type", type=str, default="")
parser.add_argument("--name", type=str, default=None)
parser.add_argument("--rank_path", type=str, default=None)
parser.add_argument("--model", type=str, default="TRAPT", help="TRAPT/H3K27ac/ATAC")
parser.add_argument("--source", type=str, default=None)
parser.add_argument("--source_exclude", type=bool, default=False)
args = parser.parse_args()

rank_path = args.rank_path
output_dir = args.output_dir
types = args.type.split(",")
name = args.name
sort_col = {
    "TRAPT":"TR activity",
    "H3K27ac":"RP_TR_H3K27ac_auc",
    "ATAC":"RP_TR_ATAC_auc",
    "TR":"RP_TR_auc"
}[args.model]

my_tr = os.listdir(output_dir)
rank = []
top_sum = 0
for t in types:
    dirs = glob(f"{output_dir}/*{t}")
    for dir in sorted(dirs, reverse=False):
        if os.path.exists("%s/TR_detail.txt" % dir):
            TR = re.findall(rf"{output_dir}/(.*?)@.*$", dir)[0]
            tr = re.findall(rf"{output_dir}/(.*?@.*$)", dir)[0]
            if tr not in my_tr or "@" not in dir:
                continue
            output = re.findall(r"/(.*)", dir)[0]
            try:
                summary_data = pd.read_csv("%s/TR_detail.txt" % dir, sep="\t")
                if args.source != None:
                    summary_data = summary_data[
                        summary_data["Source"].str.upper().map(lambda s: s != args.source.upper() if args.source_exclude else s == args.source.upper())
                    ]
                    summary_data = summary_data.reset_index(drop=True)
                summary_data = summary_data.sort_values(
                    sort_col, ascending=False, ignore_index=True
                )
                pass
            except Exception as e:
                print(dir, e)
                continue
            summary_data = summary_data[["tr_base"]]
            summary_data = summary_data.drop_duplicates(ignore_index=True).reset_index()
            summary_data.columns = ["Rank", "TR"]
            summary_data["Rank"] += 1
            if TR not in summary_data["TR"].values.tolist():
                rank_score = 2000
            else:
                rank_score = summary_data.loc[summary_data["TR"] == TR, "Rank"].values[
                    0
                ]
            if rank_score > 1363:
                rank_score = 1364
            rank.append([tr, TR, rank_score])
            if rank_score <= 10:
                top_sum += 1
                print(output, rank_score)

print("Total:", len(rank), "top_sum: ", top_sum)
rank = pd.DataFrame(rank, columns=["id", "tr", "rank"])
rank = rank.sort_values(["rank","tr"])
rank.to_csv(f"{rank_path}/rank_{name}.csv", index=False)

print(rank.query("rank <= 10").groupby("tr").agg({"rank": len}))
