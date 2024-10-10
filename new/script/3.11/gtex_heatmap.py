import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

data = []
tissues = []
tissue_top10 = []
top10_TRs = set()
top100_TRs = set()
for tissue in os.listdir("new/result/3.11/Fig. 6/GTEx/output"):
    tissues.append(tissue)
    summary_data = pd.read_csv('new/result/3.11/Fig. 6/GTEx/output/%s/TR_detail.txt' % tissue, sep='\t')
    summary_data = summary_data[['tr_base']]
    summary_data = summary_data.drop_duplicates(['tr_base'],ignore_index=True).reset_index()
    summary_data.columns = ['MMR','TR']
    summary_data['MMR'] = 1/(summary_data['MMR'] + 1)
    top100 = summary_data.head(100).TR
    top100_TRs.update(top100)
    top10 = summary_data.head(10).TR
    top10_TRs.update(top10)
    top10 = ", ".join(top10)
    summary_data = summary_data.sort_values("TR")
    data.append(summary_data['MMR'].values)
    tissue_top10.append([tissue,top10])

tissue_top10 = pd.DataFrame(tissue_top10,columns=["Tissue","Top 10 TRs"])
data = pd.DataFrame(data,columns=summary_data.TR,index=tissues)
data = ((data.T - data.min(axis=1)) / (data.max(axis=1) - data.min(axis=1))).T

data = data.iloc[:,data.columns.map(lambda x:x in top10_TRs)]
tr_index = np.where(np.in1d(list(top10_TRs),data.columns))[0] + 1

cmap1 = LinearSegmentedColormap.from_list('my_cmap', ['#fff7f4','#fde5d9'], N=10)
cmap10 = LinearSegmentedColormap.from_list('my_cmap', ['#fde5d9','#ec362b','#7c2a2f'], N=100)

cmap1 = LinearSegmentedColormap.from_list('my_cmap', ['#f6fafe','#e7f1f9'], N=10)
cmap10 = LinearSegmentedColormap.from_list('my_cmap', ['#e7f1f9','#4598c4','#004b89'], N=100)
cmap = [(cmap10(i) if i >= 10 else cmap1(i)) for i in range(100)]
sns.clustermap(
    data,
    rasterized=False,
    cmap=cmap,
    figsize=(26, 10)
)
plt.savefig("new/result/3.11/Fig. 6/figure/GTEx-heatmap.svg")
plt.clf()