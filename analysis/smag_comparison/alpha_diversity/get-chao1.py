import sys

import polars as pl

otu = sys.argv[1]

otu = pl.read_csv(otu)

chao_per_gene = (
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
        .alias("chao1")
    )
)

trim_fraction = 0.1

trimmed = chao_per_gene.group_by("sample").agg(
    [
        # 10th and 90th percentile
        pl.col("chao1").quantile(trim_fraction).alias("lower_cut"),
        pl.col("chao1").quantile(1 - trim_fraction).alias("upper_cut"),
    ]
)

# join back to original to filter per group
trimmed_mean = (
    chao_per_gene.join(trimmed, on="sample")
    # keep only values within [lower_cut, upper_cut]
    .filter(
        (pl.col("chao1") >= pl.col("lower_cut"))
        & (pl.col("chao1") <= pl.col("upper_cut"))
    )
    # then group again to compute mean
    .group_by("sample")
    .agg([pl.col("chao1").mean().alias("trimmed_mean_chao1")])
)

trimmed_mean.write_csv("chao1.csv")
