import sys

import polars as pl

otu = sys.argv[1]

# follow the same convention as get-shannon.py (tab-separated)
otu = pl.read_csv(otu)

# aggregate counts per sample x taxonomy
df_abund = otu.group_by(["sample", "taxonomy"]).agg(
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
simpson_index = df_simpson_terms.group_by("sample").agg(
    [(1 - pl.col("p_i_squared").sum()).alias("Simpson_index")]
)

simpson_index.write_csv("simpson_index.csv")
