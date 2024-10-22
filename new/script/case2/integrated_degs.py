import pandas as pd

input = "new/result/case2"
output = "new/result/case2/data.csv"

diffs = (
  ("HSC-MPPs","LMPs"),
  ("LMPs","MEMPs"),
  ("MEMPs","Megakaryocytes"),
  ("MEMPs","Erythroid cells"),
  ("MEMPs","Mast cells"),
  ("LMPs","Monocytes"),
  ("LMPs","GPs"),
  ("GPs","Granulocytes"),
  ("LMPs","pDCs"),
  ("LMPs","NK cells"),
  ("LMPs","Pro-B cells"),
  ("Pro-B cells","Pre-B cells"),
  ("Pre-B cells","Mature B cells")
)

data = []
for item in diffs:
    print(item)
    h = item[1]
    c = item[0]
    markers = pd.read_csv(f"{input}/output/markers_top200-[{h}-{c}].csv",index_col=0)
    markers["cluster"] = f"{h}\n-{c}"
    markers["gene"] = markers.index
    data.append(markers)

data = pd.concat(data)
data['label'] = data.apply(lambda row: 'Up' if (row['avg_log2FC'] > 1 and row['p_val_adj'] < 0.05) else
                                           'None' if (abs(row['avg_log2FC']) < 1 or row['p_val_adj'] > 0.05) else
                                           'Down' if (row['avg_log2FC'] < -1 and row['p_val_adj'] < 0.05) else
                                           'None', axis=1)
data["gene"] = data["gene"].map(lambda x:x.split("(")[0])
data.to_csv(output)