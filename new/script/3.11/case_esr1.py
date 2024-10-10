import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import networkx as nx

# ESR1 KD in MCF7
sns.set_style(None, {"font.sans-serif": "Arial"})
color_palette = [
    "#e64b35",
    "#b09c85",
    "#4dbbd5",
    "#00a087",
    "#b3de69",
    "#d9d9d9",
    "#3c5488",
    "#f39b7f",
    "#fdb462",
    "#8491b4",
]
data = pd.read_csv("new/result/3.11/Fig. 4/ESR1-KD-MCF7.csv", sep="\t")

f, ax = plt.subplots(figsize=(12, 6))
sns.despine(f, ax, left=False, bottom=False)
sns.scatterplot(
    x="down-regulated gene sets",
    y="up-regulated gene sets",
    hue="Type",
    hue_order=["TF", "TcoF", "CR"],
    size="score",
    palette=["#e64b35", "#4dbbd5", "#00a087"],
    data=data,
    ax=ax,
)
quadrant = data[data["top"]].reset_index(drop=True)
for x in range(len(quadrant)):
    plt.annotate(
        quadrant.loc[x, "tr_base"],
        (
            quadrant.loc[x, "down-regulated gene sets"],
            quadrant.loc[x, "up-regulated gene sets"],
        ),
    )
plt.title("ESR1 KD in MCF7")
plt.savefig("new/result/3.11/Fig. 4/figure/ESR1-KD-MCF7.svg")
plt.clf()

# STRING
gene_interact_data = pd.read_csv(
    f"new/result/3.11/Fig. 4/9606.protein.links.v12.0.txt.gz", sep=" "
)
string_info = pd.read_csv("new/result/3.11/Fig. 4/9606.protein.info.v12.0.txt.gz", "\t")
string_info = pd.merge(
    string_info,
    data.loc[data.top.values, ["tr_base"]],
    left_on="preferred_name",
    right_on="tr_base",
    how="inner",
)
inner_string_id = set(string_info["#string_protein_id"])
id_genes = dict(zip(*string_info[["#string_protein_id", "preferred_name"]].T.values))
gene_interact_data = gene_interact_data[
    np.all(
        [
            gene_interact_data["protein1"].map(lambda x: x in inner_string_id),
            gene_interact_data["protein2"].map(lambda x: x in inner_string_id),
        ],
        axis=0,
    )
]
gene_interact_data["protein1"] = gene_interact_data["protein1"].map(id_genes.get)
gene_interact_data["protein2"] = gene_interact_data["protein2"].map(id_genes.get)
gene_interact_data = gene_interact_data[
    ~(gene_interact_data["protein1"] == gene_interact_data["protein2"])
]
h = set()
gene_interact_data = gene_interact_data[
    gene_interact_data.filter(items=["protein1", "protein2"]).apply(
        lambda x: (
            False if "-".join(sorted(x)) in h else (h.add("-".join(sorted(x))) or True)
        ),
        axis=1,
    )
]
gene_interact_data["combined_score"] = (
    gene_interact_data["combined_score"] / gene_interact_data["combined_score"].max()
)
G = nx.from_pandas_edgelist(
    gene_interact_data, "protein1", "protein2", edge_attr="combined_score"
)
pos = nx.spring_layout(G, k=0.15, iterations=20, scale=1)
H = nx.Graph(G)
node_sizes = np.array([len(list(G.neighbors(n))) for n in H.nodes()])
node_sizes = node_sizes / node_sizes.max() * 100
edgewidth = [G.get_edge_data(u, v).get("combined_score") * 2 for u, v in H.edges()]
nodes = nx.draw_networkx_nodes(
    G, pos, node_size=node_sizes, node_color="#210070", alpha=0.7
)
edges = nx.draw_networkx_edges(
    G, pos, node_size=node_sizes, width=edgewidth, alpha=0.3, edge_color="m"
)
label_options = {"ec": "k", "fc": "white", "alpha": 0.7}
nx.draw_networkx_labels(G, pos, font_size=7, verticalalignment="bottom")
plt.savefig(f"new/result/3.11/Fig. 4/figure/ppi.svg", dpi=1000)
plt.clf()


# GTEx
"""
wget -c https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/gtex_gene_expected_count.gz
wget -c https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/GTEX_phenotype.gz
wget -c https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/probeMap%2Fgencode.v23.annotation.gene.probemap
"""
expr = pd.read_csv(f"new/result/3.11/Fig. 4/gtex_gene_expected_count.gz",sep="\t",index_col=0).dropna()

probe_map = pd.read_csv("new/result/3.11/Fig. 4/gencode.v23.annotation.gene.probemap",sep="\t")
probe_gene_dict = dict(probe_map[["id","gene"]].values)
expr.index = expr.index.map(probe_gene_dict.get)

phenotype = pd.read_csv(f"new/result/3.11/Fig. 4/GTEX_phenotype.gz", sep="\t")
sample = pd.DataFrame(expr.columns, columns=["Sample"])
sample = sample.merge(phenotype, on="Sample")
expr = expr.loc[:, sample["Sample"].values]
tissue = list(
    sample.loc[
        np.all([sample._primary_site == "Breast", sample._gender == "female"], axis=0),
        "Sample",
    ]
)
expr_genes = set(expr.index)
sorted_tfs = data[data.top.values].sort_values("score", ascending=False)["tr_base"]
sorted_tfs = sorted_tfs[sorted_tfs.map(lambda x: x in expr_genes)]
expr_tissue = expr.loc[sorted_tfs, tissue]
corr = np.corrcoef(expr_tissue)
corr = pd.DataFrame(corr, index=expr_tissue.index, columns=expr_tissue.index)
d = corr.iloc[:, :]
R = d.values.max() - d.values.min()
N1 = int(abs(d.values.min()) * 100)
N2 = int(abs(d.values.max()) * 100)
cmap1 = LinearSegmentedColormap.from_list(
    "my_cmap", ["blue", "#f0f0f0", "#f0f0f0"], N=100
)
cmap2 = LinearSegmentedColormap.from_list(
    "my_cmap", ["#f0f0f0", "#f0f0f0", "red"], N=100
)
cmap = [(cmap1(i + 100 - N1) if i <= N1 else cmap2(i - N1)) for i in range(N1 + N2)]

plt.figure(figsize=(10, 10))
sns.clustermap(
    d,
    rasterized=False,
    cmap=cmap,
    metric="correlation", 
    method="single"
)
plt.savefig("new/result/3.11/Fig. 4/GTEx.svg")
plt.clf()

# TCGA
"""
wget -c https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BRCA.htseq_counts.tsv.gz
wget -c https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BRCA.htseq_fpkm.tsv.gz
wget -c https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v22.annotation.gene.probeMap
wget -c https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BRCA.GDC_phenotype.tsv.gz
"""
expr = pd.read_csv("new/result/3.11/Fig. 4/TCGA-BRCA.htseq_counts.tsv.gz",sep="\t",index_col=0)
probe_map = pd.read_csv("new/result/3.11/Fig. 4/gencode.v22.annotation.gene.probeMap",sep="\t")
phenotype = pd.read_csv("new/result/3.11/Fig. 4/TCGA-BRCA.GDC_phenotype.tsv.gz",sep="\t",index_col=0)

sample = pd.DataFrame(expr.columns, columns=["Sample"])
sample = sample.merge(phenotype, left_on="Sample",right_on="submitter_id.samples")
expr = expr.loc[:, sample["Sample"].values]

phenotype_female = sample.loc[(sample["gender.demographic"]=="female"),"Sample"]
expr = expr[phenotype_female]

probe_gene_dict = dict(probe_map[["id","gene"]].values)
expr.index = expr.index.map(probe_gene_dict.get)

col_names = list(expr.columns)

def process_colname(x):
    s = x.split("-")  # 分割列名
    submitter_id = "-".join(s[:3])  # 提取前3个部分作为 submitter_id
    tumor = int(s[3][:-1]) < 9  # 判断第四部分中的数字是否小于9，表示是否为肿瘤
    return [submitter_id, tumor]

phenotype = pd.DataFrame([process_colname(x) for x in col_names], columns=["submitter_id", "tumor"])
expr.columns = phenotype["submitter_id"]
expr_tissue = expr[phenotype["submitter_id"][phenotype["tumor"]]]

expr_genes = set(expr.index)
sorted_tfs = data[data.top.values].sort_values("score", ascending=False)["tr_base"]
sorted_tfs = sorted_tfs[sorted_tfs.map(lambda x: x in expr_genes)]
expr_tissue = expr_tissue.loc[sorted_tfs]

# expr_tissue = np.log2(expr_tissue + 1)

corr = np.corrcoef(expr_tissue)
corr = pd.DataFrame(corr, index=expr_tissue.index, columns=expr_tissue.index)
d = corr.iloc[:, :]
R = d.values.max() - d.values.min()
N1 = int(abs(d.values.min()) * 100)
N2 = int(abs(d.values.max()) * 100)
cmap1 = LinearSegmentedColormap.from_list(
    "my_cmap", ["blue", "#f0f0f0", "#f0f0f0"], N=100
)
cmap2 = LinearSegmentedColormap.from_list(
    "my_cmap", ["#f0f0f0", "#f0f0f0", "red"], N=100
)
cmap = [(cmap1(i + 100 - N1) if i <= N1 else cmap2(i - N1)) for i in range(N1 + N2)]

plt.figure(figsize=(10, 10))
sns.clustermap(
    d,
    rasterized=False,
    cmap=cmap,
    metric="correlation", 
    method="single"
)
plt.savefig("new/result/3.11/Fig. 4/TCGA.svg")
plt.clf()