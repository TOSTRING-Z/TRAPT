import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import seaborn as sns

dataset = pd.read_csv("new/result/3.11/Fig. 5/dataset.csv")

dataset_melt = pd.melt(
    frame=dataset[["TR", "Fine_mapping_snps", "Background_snps"]],
    id_vars=["TR"],
    value_vars=["Fine_mapping_snps", "Background_snps"],
    var_name="Type",
    value_name="Count",
    ignore_index=True,
)


import matplotlib.pyplot as plt
import seaborn as sns

fig, ax1 = plt.subplots(figsize=(5, 10))
sns.barplot(
    data=dataset_melt,
    x="Count",
    y="TR",
    hue="Type",
    order=dataset["TR"],
    ax=ax1,
    palette="pastel",
)

ax1.set_xlim(ax1.get_xlim())
ax1.xaxis.tick_top()
ax2 = ax1.twiny()
ax2.plot(
    dataset["Enrichment_score"],
    dataset["TR"],
    "--",
    linewidth=1.5,
    color="#ff4040",
)
ax2.plot(
    dataset["Enrichment_score"],
    dataset["TR"],
    "o",
    markersize=6,
    linewidth=2,
    markeredgecolor="black",
    markeredgewidth=0.5,
    color=sns.color_palette("pastel")[3],
)

plt.savefig("new/result/3.11/Fig. 5/figure/enrichment.svg")
plt.clf()