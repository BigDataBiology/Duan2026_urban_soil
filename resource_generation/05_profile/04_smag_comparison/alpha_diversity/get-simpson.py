import sys
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

# Use only bacterial and both domains for Simpson computation
KEEP_DOMAINS = ["both", "bacteria"]
otu_sm = otu_sm.filter(pl.col("domain").is_in(KEEP_DOMAINS))
otu_sh = otu_sh.filter(pl.col("domain").is_in(KEEP_DOMAINS))


def compute_simpson(otu, aggregation):
    # aggregate counts per sample x taxonomy
    df_abund = otu.group_by(["sample", "gene"]).agg(
        pl.col("num_hits").sum().alias("abundance")
    )

    # compute total counts per sample and relative abundance p_i
    df_rel = df_abund.with_columns(
        [
            pl.col("abundance").sum().over("sample").alias("total_count"),
        ]
    ).with_columns([(pl.col("abundance") / pl.col("total_count")).alias("p_i")])

    # compute sum of squared relative abundances per sample
    df_simpson_terms = df_rel.with_columns([(pl.col("p_i") ** 2).alias("p_i_squared")])

    # Simpson index (1 - sum p_i^2). Higher values -> higher diversity
    simpson = df_simpson_terms.group_by("sample").agg(
        [(1 - pl.col("p_i_squared").sum()).alias("Simpson")]
    )

    if aggregation == "mean":
        agg_f = pl.col("Simpson").mean()
    elif aggregation == "tmean10":
        agg_f = (
            pl.col("Simpson")
            .clip(pl.col("Simpson").quantile(0.1), pl.col("Simpson").quantile(0.9))
            .mean()
        )
    elif aggregation == "median":
        agg_f = pl.col("Simpson").median()
    else:
        raise ValueError(f"Unknown aggregation method: {aggregation}")
    return simpson.group_by("sample").agg(agg_f.alias(f"Simpson"))


fig, axes = plt.subplots(1, 3, figsize=(8, 4), sharey=True)
data = []
for ax, agg in zip(axes, ["mean", "tmean10", "median"]):
    ax.clear()
    print(f"Aggregation method: {agg}")
    simpson_smag = compute_simpson(otu_sm, agg)
    simpson_sh = compute_simpson(otu_sh, agg)
    stat, p_value = stats.mannwhitneyu(
        simpson_smag["Simpson"], simpson_sh["Simpson"], alternative="two-sided"
    )
    print(f"Mann-Whitney U test: statistic={stat}, p-value={p_value}")

    cur = pl.concat(
        [
            simpson_smag.with_columns(pl.lit("Smag").alias("Group")),
            simpson_sh.with_columns(pl.lit("Urban soil").alias("Group")),
        ]
    ).with_columns(pl.lit(agg).alias("Aggregation"))
    sns.boxplot(
        x="Group",
        y=f"Simpson",
        ax=ax,
        data=cur,
        fill=False,
        color="black",
        fliersize=0,
        order=["Urban soil", "Smag"],
    )
    sns.stripplot(
        x="Group",
        y=f"Simpson",
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
plt.savefig(f"figures/simpson.svg")

data = pl.concat(data)
data = data.pivot(index="sample", on="Aggregation", values="Simpson")

makedirs("results", exist_ok=True)
data.write_csv("results/simpson.tsv", separator="\t")

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
shs = pl.DataFrame(shs, schema=["gene", "Simpson"], orient="row")
recomp = shs.select(
    pl.col("Simpson").median().alias("median"),
    pl.col("Simpson").mean().alias("mean"),
    pl.col("Simpson")
    .clip(pl.col("Simpson").quantile(0.1), pl.col("Simpson").quantile(0.9))
    .mean()
    .alias("tmean10"),
)
orig = data.filter(pl.col("sample") == SPOT_SAMPLE).select(cs.float())
recomp = recomp[orig.columns]
print("Difference between original and recomputed values:")
print(orig - recomp)
