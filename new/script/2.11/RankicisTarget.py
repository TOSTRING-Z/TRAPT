import os
import re
from glob import glob
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--match_dir", type=str, default="new/result/2.11/output-Lisa")
parser.add_argument("--type", type=str, default="down,up")
parser.add_argument("--input_path", type=str, default="new/result/2.11/chea3")
parser.add_argument("--output_path", type=str, default="new/result/2.11/files")
args = parser.parse_args()

my_tr = os.listdir(args.match_dir)
types = args.type.split(",")

rank = []
for t in types:
    dirs = glob(f'{args.input_path}/{t}/*')
    for dir in dirs:
        file = f'{dir}/icistarget/statistics.tbl'
        TR = re.findall(rf'/{t}/(.*?)@', file)[0]
        tr = re.findall(rf'/{t}/(.*?@.*?)/icistarget', file)[0]
        if tr not in my_tr:
            continue
        summary_data = pd.read_csv(file, sep='\t')[['FeatureDescription']]
        summary_data.columns = ['factor']
        summary_data['factor'] = (
            summary_data['factor']
            .apply(
                lambda x: x.split(' ')[-1] if x.startswith('ChIP') else x.split(' ')[0]
            )
            .str.strip()
        )
        summary_data = summary_data.drop_duplicates(
            ['factor'], ignore_index=True
        ).reset_index()
        summary_data.columns = ['Rank', 'TR']
        summary_data['Rank'] += 1
        if TR not in summary_data['TR'].values.tolist():
            rank_score = 2000
        else:
            rank_score = summary_data.loc[summary_data['TR'] == TR, 'Rank'].values[0]
        if rank_score > 262:
            rank_score = 263
        rank.append([tr, TR, rank_score])
print('TR Total', len(rank))
rank = pd.DataFrame(rank, columns=["id", "tr", "rank"])
rank.to_csv('%s/rank_i-cisTarget.csv' % args.output_path, index=False)

print(rank.query('rank <= 10').groupby('tr').agg({'rank': len}))
