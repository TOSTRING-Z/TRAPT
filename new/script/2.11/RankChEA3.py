import os
import re
from glob import glob
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--match_dir", type=str, default="new/result/2.11/output-Lisa")
parser.add_argument("--input_path", type=str, default="new/result/2.11/chea3")
parser.add_argument("--type", type=str, default="down,up")
parser.add_argument("--output_path", type=str, default="new/result/2.11/files")
args = parser.parse_args()

my_tr = os.listdir(args.match_dir)
types = args.type.split(",")

rank = []
for t in types:
    files = glob(f'{args.input_path}/{t}/*{t}*$Integrated--meanRank.tsv')
    for file in files:
        TR = re.findall(rf'/{t}/(.*?)@', file)[0]
        tr = re.findall(rf'/{t}/(.*?@.*?)\$Integrated--meanRank.tsv', file)[0]
        if tr not in my_tr:
            continue
        summary_data = pd.read_csv(file, sep='\t')[['TF']]
        summary_data = summary_data.drop_duplicates(
            ['TF'], ignore_index=True
        ).reset_index()
        summary_data.columns = ['Rank', 'TR']
        summary_data['Rank'] += 1
        if TR not in summary_data['TR'].values.tolist():
            rank_score = 2000
        else:
            rank_score = summary_data.loc[summary_data['TR'] == TR, 'Rank'].values[0]
        if rank_score > 1601:
            rank_score = 1602
        rank.append([TR, rank_score])
print('TR Total', len(rank))
rank = pd.DataFrame(rank, columns=['tr', 'rank'])
rank.to_csv('%s/rank_ChEA3.csv' % args.output_path, index=False)

print(rank.query('rank <= 10').groupby('tr').agg({'rank': len}))
