from scipy import stats
from matplotlib import pyplot as plt
import seaborn as sns
from os import makedirs

import polars as pl
import polars.selectors as cs


hmms_and_names = pl.read_csv("hmms_and_names", separator="\t")
hmms_and_names = hmms_and_names.with_columns(pl.col("type").alias("domain"))

# Note that one file is comma-separated and the other is tab-separated

otu_sm = "OTU/PRJNA983538.csv"
otu_sm = pl.read_csv(otu_sm, separator=",")
otu_sm = otu_sm.join(hmms_and_names["name", "domain"], left_on="gene", right_on="name")

otu_sh = "OTU/shanghai_soil_otu.tsv"
otu_sh = pl.read_csv(otu_sh, separator="\t")
otu_sh = otu_sh.join(hmms_and_names["name", "domain"], left_on="gene", right_on="name")

# Use only bacterial and both domains for Chao1 computation
KEEP_DOMAINS = ["both", "bacteria"]
otu_sm = otu_sm.filter(pl.col("domain").is_in(KEEP_DOMAINS))
otu_sh = otu_sh.filter(pl.col("domain").is_in(KEEP_DOMAINS))


def compute_chao1(otu, aggregation):
    chao1 = (
        otu.group_by("sample", "gene")
        .agg(
            [
                pl.len().alias("S_obs"),
                (pl.col("num_hits") == 1).sum().alias("F1"),
                (pl.col("num_hits") == 2).sum().alias("F2"),
            ]
        )
        .with_columns(
            pl.when(pl.col("F2") == 0)
            .then(pl.col("S_obs") + (pl.col("F1") * (pl.col("F1") - 1)) / 2)
            .otherwise(pl.col("S_obs") + (pl.col("F1") ** 2) / (2 * pl.col("F2")))
            .alias("Chao1")
        )
    )

    if aggregation == "mean":
        agg_f = pl.col("Chao1").mean()
    elif aggregation == "tmean10":
        agg_f = (
            pl.col("Chao1")
            .clip(pl.col("Chao1").quantile(0.1), pl.col("Chao1").quantile(0.9))
            .mean()
        )
    elif aggregation == "median":
        agg_f = pl.col("Chao1").median()
    else:
        raise ValueError(f"Unknown aggregation method: {aggregation}")
    return chao1.group_by("sample").agg(agg_f.alias(f"Chao1"))


fig, axes = plt.subplots(1, 3, figsize=(8, 4), sharey=True)
data = []
for ax, agg in zip(axes, ["mean", "tmean10", "median"]):
    ax.clear()
    print(f"Aggregation method: {agg}")
    chao1_smag = compute_chao1(otu_sm, agg)
    chao1_sh = compute_chao1(otu_sh, agg)
    stat, p_value = stats.mannwhitneyu(
        chao1_smag["Chao1"], chao1_sh["Chao1"], alternative="two-sided"
    )
    print(f"Mann-Whitney U test: statistic={stat}, p-value={p_value}")

    cur = pl.concat(
        [
            chao1_smag.with_columns(pl.lit("Smag").alias("Group")),
            chao1_sh.with_columns(pl.lit("Urban soil").alias("Group")),
        ]
    ).with_columns(pl.lit(agg).alias("Aggregation"))
    sns.boxplot(
        x="Group",
        y=f"Chao1",
        ax=ax,
        data=cur,
        fill=False,
        color="black",
        fliersize=0,
        order=["Urban soil", "Smag"],
    )
    sns.stripplot(
        x="Group",
        y=f"Chao1",
        ax=ax,
        data=cur,
        hue="Group",
        size=2,
        jitter=True,
        order=["Urban soil", "Smag"],
        alpha=0.7,
    )

    ax.set_title(f"{agg}")
    data.append(cur)
sns.despine(fig, trim=True)
fig.tight_layout()
makedirs("figures", exist_ok=True)
plt.savefig(f"figures/chao1.svg")

data = pl.concat(data)
data = data.pivot(index="sample", on="Aggregation", values="Chao1")
print(data)

makedirs("results", exist_ok=True)
data.write_csv("results/chao1.tsv", separator="\t")

print("\n# Spearman correlation between different aggregation methods\n")
print("Agg1\tAgg2\tSpearman correlation")
for col in ["mean", "median", "tmean10"]:
    for col2 in ["mean", "median", "tmean10"]:
        if col > col2:
            print(
                f"{col}\t{col2}\t{stats.spearmanr(data[col], data[col2]).statistic:.4f}"
            )


print("# Spot check computation")

SPOT_SAMPLE = "CPSNJ01_350"
shs = []
sel_1sample = otu_sh.filter(pl.col("sample") == SPOT_SAMPLE)
for g in set(sel_1sample["gene"]):
    sel = sel_1sample.filter(pl.col("gene") == g)
    shs.append((g, stats.entropy(sel["num_hits"])))
shs = pl.DataFrame(shs, schema=["gene", "Chao1"], orient="row")
recomp = shs.select(
    pl.col("Chao1").median().alias("median"),
    pl.col("Chao1").mean().alias("mean"),
    pl.col("Chao1")
    .clip(pl.col("Chao1").quantile(0.1), pl.col("Chao1").quantile(0.9))
    .mean()
    .alias("tmean10"),
)
orig = data.filter(pl.col("sample") == SPOT_SAMPLE).select(cs.float())
recomp = recomp[orig.columns]
print("Difference between original and recomputed values:")
print(orig - recomp)
