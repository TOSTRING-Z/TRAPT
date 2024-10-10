# nohup rsync -P --rsh=ssh zhangguorui@10.147.17.156:/linuxdata3/sc/bin/python_infertf/knockTF/up/*.lisa.tsv .

import os
import re
from glob import glob
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--match_dir", type=str, default="new/result/2.11/output-Lisa")
parser.add_argument("--type", type=str, default="down,up")
parser.add_argument("--input_path", type=str, default="new/result/2.11/lisa")
parser.add_argument("--output_path", type=str, default="new/result/2.11/files")
args = parser.parse_args()

my_tr = os.listdir(args.match_dir)
types = args.type.split(",")

rank = []
for t in types:
    files = glob(f'{args.input_path}/{t}/*.lisa.tsv')
    for file in sorted(files, reverse=True):
        TR = re.findall(rf'/{t}/(.*?)@.*?\.lisa\.tsv', file)[0]
        tr = re.findall(rf'/{t}/(.*?@.*?)\.txt\.lisa\.tsv', file)[0]
        if tr not in my_tr:
            continue
        summary_data = pd.read_csv(file, sep='\t').sort_values('summary_p_value')[
            ['factor']
        ]
        summary_data = summary_data.drop_duplicates(ignore_index=True).reset_index()
        summary_data.columns = ['rank', 'factor']
        summary_data['rank'] += 1
        if TR not in summary_data['factor'].values.tolist():
            rank_score = 2000
        else:
            rank_score = summary_data.loc[summary_data['factor'] == TR, 'rank'].values[0]
        if rank_score > 1073:
            rank_score = 1074
        rank.append([TR, rank_score])

rank = pd.DataFrame(rank, columns=['tr', 'rank'])
rank.to_csv('%s/rank_Lisa.csv' % args.output_path, index=False)
top_sum = len(rank.query('rank <= 10'))
print("Total:", len(rank), "top_sum: ", top_sum)
print(rank.query('rank <= 10').groupby('tr').agg({'rank': len}))
