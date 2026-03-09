import sys

import polars as pl

otu = sys.argv[1]

otu = pl.read_csv(otu, separator="\t")

df_abund = otu.group_by(["sample", "taxonomy"]).agg(
    pl.col("num_hits").sum().alias("abundance")
)

# Step 2: For each sample, compute total counts and relative abundance p_i
df_shannon = df_abund.with_columns(
    [
        # total count per sample
        pl.col("abundance").sum().over("sample").alias("total_count")
    ]
).with_columns([(pl.col("abundance") / pl.col("total_count")).alias("p_i")])

# Step 3: Calculate Shannon term p_i * ln(p_i), treating p_i=0 as zero
df_shannon_terms = df_shannon.with_columns(
    [(-pl.col("p_i") * pl.col("p_i").log()).alias("shannon_term")]
)

# Step 4: Sum Shannon terms per sample
shannon_index = df_shannon_terms.group_by("sample").agg(
    [pl.col("shannon_term").sum().alias("Shannon_index")]
)

shannon_index.write_csv("shannon_index.csv")
