# nohup rsync -P --rsh=ssh zhangguorui@172.23.234.192:/linuxdata3/sc/bin/python_infertf/knockTF/up/*.lisa.tsv .

import os
import re
from glob import glob
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--match_dir", type=str, default="new/result/2.11/output-Lisa")
parser.add_argument("--type", type=str, default="down,up")
parser.add_argument("--input_path", type=str, default="new/result/2.11/bart")
parser.add_argument("--output_path", type=str, default="new/result/2.11/files")
args = parser.parse_args()

my_tr = os.listdir(args.match_dir)
types = args.type.split(",")

rank = []
for t in types:
    files = glob(f"{args.input_path}/{t}/*_bart_results.txt")
    for file in sorted(files, reverse=True):
        TR = re.findall(rf"/{t}/(.*?)@.*?\.txt", file)[0]
        tr = re.findall(rf"/{t}/(.*?@.*?)_bart_results\.txt", file)[0]
        if tr not in my_tr:
            continue
        # print(TR,file)
        summary_data = pd.read_csv(file, sep="\t").iloc[:, 0].reset_index()
        summary_data.columns = ["Rank", "TR"]
        summary_data["Rank"] += 1
        if TR not in summary_data["TR"].values.tolist():
            rank_score = 2000
        else:
            rank_score = summary_data.loc[summary_data["TR"] == TR, "Rank"].values[0]
        if rank_score > 915:
            rank_score = 916
        rank.append([tr, TR, rank_score])

rank = pd.DataFrame(rank, columns=["id", "tr", "rank"])
rank.to_csv("%s/rank_BART.csv" % args.output_path, index=False)
top_sum = len(rank.query('rank <= 10'))
print("Total:", len(rank), "top_sum: ", top_sum)
print(rank.query("rank <= 10").groupby("tr").agg({"rank": len}))
