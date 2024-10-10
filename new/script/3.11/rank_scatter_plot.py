import os
import re
import sys
os.environ["CUDA_VISIBLE_DEVICES"] = ""
sys.path.append(".")
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.preprocessing import minmax_scale
from tqdm import tqdm
import json
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors
from TRAPT.DLFS import FeatureSelection
from TRAPT.Tools import Args, RPMatrix
from matplotlib.path import Path

import math

from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# # # # # Sunburst-

tf_family = pd.read_csv("input/Homo_sapiens_TF.txt",sep="\t")
data = pd.DataFrame(os.listdir("output/KnockTFv1/"), columns=["name"])
data["TR"] = data["name"].map(lambda x: x.split("@")[0])
data["type"] = data["name"].map(lambda x: x.split("_")[-1])
data = data.groupby(["TR", "type"]).agg(list)

def get_rank(name):
    print(name)
    summary_data = pd.read_csv(f"output/KnockTFv1/{name}/TR_detail.txt", sep="\t")
    summary_data = summary_data[["tr_base"]]
    summary_data = summary_data.drop_duplicates(ignore_index=True).reset_index()
    summary_data.columns = ["Rank", "TR"]
    summary_data["Rank"] += 1
    summary_data = summary_data.set_index("TR")
    tr = name.split("@")[0]
    if tr not in summary_data.index:
        rank_score = 1500
    else:
        rank_score = summary_data.loc[tr, "Rank"]
    return rank_score

data["rank"] = data["name"].map(lambda x: [get_rank(name) for name in x])
data["MMR"] = data["rank"].map(lambda x: [1 / float(rank) for rank in x])
data = data.reset_index().merge(tf_family,how="left",left_on="TR",right_on="Symbol")
data.loc[data["Family"].isna(),"Family"] = "Others"

color = ["#c58cad","#c48daa","#5151e7","#d66f80",
        "#670006","#99c1ff","#dad2ea","#b55b74","#c36559",
        "#e77326","#87a3e2","#e00100","#3383fe","#f282bd",
        "#a683d5","#ffca0c","#ff6f0e","#ff720e","#678c6a",
        "#56a355","#4e9c8c","#84d066","#dab92e","#c5912c",
        "#e9d1dd","#ffd924","#ddb695","#4378ae","#924966",
        "#686089","#708d81","#ff4d6d","#f7b538","#0fa3b1",
        "#8f2d56","#468189","#e26d5c","#fec3a6","#babd8d",
        "#ffd100","#758bfd","#7678ed","#be95c4","#b56576",
        "#f9dcc4","#ff70a6"]
family = data["Family"].unique()
family_color = dict(zip(family,color))
pd.DataFrame([family_color]).T.to_csv("new/result/3.11/Fig. 3/figure/family_color.txt",header=False,sep="\t")

tree = []
for _, row in data.iterrows():
    names = row['name']
    ranks = row['MMR']
    nodes = [{'name': name, 'value': rank} for name,rank in zip(names,ranks)]
    tr_node = {'name': row['TR'], 'children': nodes}
    family_node = {'name': row['Family'], 'children': [tr_node]}
    type_node = {'name': row['type'], 'children': [family_node]}
    type_match = next((item for item in tree if item['name'] == row['type']), None)
    if type_match:
        family_match = next((item for item in type_match['children'] if item['name'] == row['Family']), None)
        if family_match:
            tr_match = next((item for item in family_match['children'] if item['name'] == row['TR']), None)
            if tr_match:
                tr_match['children'].extend(nodes)
            else:
                family_match['children'].append(tr_node)
        else:
            type_match['children'].append(family_node)
    else:
        tree.append(type_node)

tree = {'data':tree}
with open("new/result/3.11/Fig. 3/figure/Sunburst-ALL.json","w") as file:
    file.write(json.dumps(tree))

for t in ["up500","down500"]:
    tree = []
    for _, row in data.query(f"type=='{t}'").iterrows():
        names = row['name']
        ranks = row['MMR']
        nodes = [{'name': name, 'value': rank} for name,rank in zip(names,ranks)]
        tr_node = {'name': row['TR'], 'children': nodes}
        family_node = {'name': row['Family'], 'children': [tr_node]}
        family_match = next((item for item in tree if item['name'] == row['Family']), None)
        if family_match:
            tr_match = next((item for item in family_match['children'] if item['name'] == row['TR']), None)
            if tr_match:
                tr_match['children'].extend(nodes)
            else:
                family_match['children'].append(tr_node)
        else:
            tree.append(family_node)
    tree = {'data':tree}
    with open(f"new/result/3.11/Fig. 3/figure/Sunburst-{t}.json","w") as file:
        file.write(json.dumps(tree))

tf_family = pd.read_csv("input/Homo_sapiens_TF.txt",sep="\t")
data = pd.DataFrame(os.listdir("output/KnockTFv1/"), columns=["name"])
data["TR"] = data["name"].map(lambda x: x.split("@")[0])
data["type"] = data["name"].map(lambda x: x.split("_")[-1])
data = data.groupby(["TR", "type"]).agg(list)

def get_rank(name):
    summary_data = pd.read_csv(f"output/KnockTFv1/{name}/TR_detail.txt", sep="\t")
    summary_data = summary_data[["tr_base"]]
    summary_data = summary_data.drop_duplicates(ignore_index=True).reset_index()
    summary_data.columns = ["Rank", "TR"]
    summary_data["Rank"] += 1
    summary_data = summary_data.set_index("TR")
    tr = name.split("@")[0]
    if tr not in summary_data.index:
        rank_score = 1500
    else:
        rank_score = summary_data.loc[tr, "Rank"]
    return rank_score

data["rank"] = data["name"].map(lambda x: [get_rank(name) for name in x])
data = data.reset_index().merge(tf_family,how="left",left_on="TR",right_on="Symbol")
data.loc[data["Family"].isna(),"Family"] = "Others"
data = data.sort_values('Family')

# 画图

data_sub = data.loc[:, ["type","Family","TR","rank"]].explode('rank')
data_sub = data_sub.query("rank!=1500")

data_sub["MMR"] = data_sub["rank"].map(lambda x: 1 / float(x))
data_sub_g = data_sub.groupby('Family').agg({"MMR":np.mean}).reset_index()
data_sub = data_sub.merge(data_sub_g.iloc[:],on="Family").sort_values('MMR_y',ascending=False)

tr_names = data_sub["TR"].drop_duplicates().tolist()
data_sub['TR'] = data_sub['TR'].map(tr_names.index)
data_sub["size"] = 6
data_ = data_sub.drop_duplicates("TR")
scale = 50
data_sub["rank"] += scale
s_max = data_sub["rank"].max()
s_min = data_sub["rank"].min()
cmap100 = LinearSegmentedColormap.from_list('my_cmap', ['#f0402f','#ffeee6'], N=100)
cmap1400 = LinearSegmentedColormap.from_list('my_cmap', ['#ffeee6','#ffffff'], N=s_max)
cmap = lambda i:cmap100(i) if i <= 100 else cmap1400(i)
data_sub["color"] = data_sub['rank'].map(cmap).map(mcolors.rgb2hex)
data_sub.loc[data_sub["type"]=="up500","rank"] *= -1

# 
fig, ax = plt.subplots(figsize=(6,6))
sc = ax.scatter(data_sub["rank"], data_sub["TR"], 
            c=data_sub["color"],s=data_sub["size"],
            cmap='viridis',alpha=0.1,edgecolors="#f1f1f1")

sc.set_edgecolor('gray')
sc.set_linewidth(0.1)
sc.set_alpha(0.8)
ax.set_xlabel(None, fontdict={'size': 16})
ax.set_ylabel(None, fontdict={'size': 16})
ax.grid(True, axis='y', linestyle='--', linewidth=0.5, alpha=0.5)
ax.grid(True, axis='x', linestyle='--', linewidth=1, alpha=0.9)
ax.set_xticks([-s_max,-s_min,s_min,s_max],[1358,1,1,1358])
ax.set_yticks(range(len(tr_names)),[None for i in tr_names],size=2)
ax.tick_params(axis='y', color='tab:green', size=0, width=0)

divider = make_axes_locatable(ax)

ax.set_ylim(-0.5-50,len(tr_names)-0.5+50)

ax.text(0.13,-0.04,"up-regulated gene sets",color="#ff7e79",fontsize=7,transform=ax.transAxes)
ax.text(0.57,-0.04,"down-regulated gene sets",color="#76d6ff",fontsize=7,transform=ax.transAxes)

data_['family_color'] = data_["Family"].map(family_color.get)
data_['x'] = 1

ax_right = divider.append_axes("right", 0.5, pad=0.1)
ax_right.xaxis.set_tick_params(labelbottom=False,color='white')

ax_right.barh(data_["TR"],data_['x'],color=data_['family_color'],edgecolor='k',linewidth=0,height=1)

cmap100 = LinearSegmentedColormap.from_list('my_cmap', ['#f0402f','#ffeee6'], N=100)
cmap1400 = LinearSegmentedColormap.from_list('my_cmap', ['#ffeee6','#ffffff'], N=s_max-100)
cmap = lambda i:cmap100(i) if i <= 100 else cmap1400(i)
cmap_ = LinearSegmentedColormap.from_list('my_cmap', [cmap(i) for i in range(s_max)], N=s_max)

ax_bar = divider.append_axes("top", 0.1, pad=.1)
ax_bar.yaxis.set_tick_params(labelbottom=False,color='white')
ax_bar.imshow(np.array([[cmap(i) for i in range(data_["rank"].max())]]), aspect='auto', cmap=cmap_)
ax_bar.set_xticks([0,99,data_["rank"].max()],[1,100,1358])
ax_bar.xaxis.tick_top()

g = data_.groupby("Family").agg({"TR":lambda x:int(np.mean(x))}).reset_index()
g = g.sort_values("TR")
ax_right.set_ylim(-0.5-50,len(tr_names)-0.5+50)


# 标记需要注释的行或列
mark_indices = list(g["TR"])
mark_labels = list(g["Family"])

t = -50
b = 50
r = len(data_['x']) - t + b
x = 1
# 绘制注释标记和连接线
for i, index in enumerate(mark_indices):
    y = int(i * (r/len(mark_labels))) + t
    color = family_color[mark_labels[i]]
    ax_right.plot([x,x+3,x+6,x+9], [index,index,y,y], color=color,linewidth=0.5)
    ax_right.text(x+9,y,mark_labels[i],color=color,fontsize=6,va="center")

ax_right.axis('off')


ax_left = divider.append_axes("left", 0.2, pad=0.01)
ax_left.xaxis.set_tick_params(labelbottom=False,color='white')
# 标记需要注释的行或列
data_sub["rank_abs"] = data_sub["rank"].abs()
sels = data_sub.query("rank_abs < 60").drop_duplicates("TR").sort_values("TR")
mark_indices = list(sels.TR)
mark_labels = [tr_names[i] for i in mark_indices]
fam_labels = list(sels.Family)

t = -50
b = 50
r = len(data_['x']) - t + b
x = 0
# 绘制注释标记和连接线
for i, index in enumerate(mark_indices):
    y = int(i * (r/len(mark_labels))) + t
    color = family_color[fam_labels[i]]
    ax_left.plot([x,x-1,x-2,x-3][::-1], [index,index,y,y][::-1], color=color,linewidth=0.5)
    ax_left.text(x-3,y,mark_labels[i],color=color,fontsize=6,va="center",ha="right")

ax_left.set_ylim(-0.5-50,len(tr_names)-0.5+50)
ax_left.axis('off')


ax_bottom = divider.append_axes("bottom", 1, pad=0.3)

def plot_colortable(colors, *, ncols=4):
    cell_width = 212
    cell_height = 22
    swatch_width = 48
    names = list(colors)
    n = len(names)
    nrows = math.ceil(n / ncols)
    ax_bottom.set_xlim(0, cell_width * 4)
    ax_bottom.set_ylim(cell_height * (nrows-0.5), -cell_height/2.)
    ax_bottom.yaxis.set_visible(False)
    ax_bottom.xaxis.set_visible(False)
    ax_bottom.set_axis_off()
    for i, name in enumerate(names):
        row = i % nrows
        col = i // nrows
        y = row * cell_height
        swatch_start_x = cell_width * col
        text_pos_x = cell_width * col + swatch_width + 7
        ax_bottom.text(text_pos_x, y, name, fontsize=7,
                horizontalalignment='left',
                verticalalignment='center')
        ax_bottom.add_patch(
            Rectangle(xy=(swatch_start_x, y-9), width=swatch_width,
                      height=18, facecolor=colors[name], edgecolor='0.7')
        )
    return fig

plot_colortable(family_color,ncols=4)

plt.savefig("new/result/3.11/Fig. 3/figure/rank_scatter.svg")
plt.clf()
