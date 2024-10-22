import pandas as pd
#  os.chdir("/data/zgr/data/TRAPT/tool/")
input = "new/result/case2"
output = "new/result/case2/data.csv"
library = "library"

tr_info = pd.read_csv(f"{library}/TRs_info.txt",sep="\t",index_col=0)
tr_set = set(tr_info.tr_base)

diffs = (
  ("HSC-MPPs","LMPs"),
  ("LMPs","Monocytes"),
  ("LMPs","GPs"),
  ("LMPs","pDCs"),
  ("LMPs","NK cells"),
  ("LMPs","Pro-B cells"),
  ("LMPs","MEMPs"),
  ("MEMPs","Megakaryocytes"),
  ("MEMPs","Mast cells"),
  ("MEMPs","Erythroid cells"),
  ("GPs","Granulocytes"),
  ("Pro-B cells","Pre-B cells"),
  ("Pre-B cells","Mature B cells")
)

data = []
for item in diffs:
    print(item)
    h = item[1]
    c = item[0]
    h_r = item[1].replace("-",".").replace(" ",".")
    c_r = item[0].replace("-",".").replace(" ",".")
    markers = pd.read_csv(f"{input}/output/markers_top200-[{h_r}-{c_r}].csv",index_col=0)
    markers["cluster"] = f"{h}\n-{c}"
    markers["gene"] = markers.index
    data.append(markers)

data = pd.concat(data)
data = data[data["gene"].map(lambda x:x in tr_set)]
data['isTR'] = data.apply(lambda row: 1 if (row['gene'] in tr_set and abs(row["avg_log2FC"]) > 0.5 and row["p_val_adj"] < 0.05) else 0, axis=1)
data["gene"] = data["gene"].map(lambda x:x.split("(")[0])
data.to_csv(output)