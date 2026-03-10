import sys

import pandas as pd
import requests

study = sys.argv[1]

samples = [
    sample["acc"]
    for sample in requests.get(
        f"https://sandpiper.qut.edu.au/api/project?model_bioproject={study}"
    ).json()["projects"]
]

chunk_size = 10_000_000

first_chunk = True

for chunk in pd.read_csv(
    "/work/microbiome/db/singlem/sandpiper1.0.0.otu_table.tsv.gz",
    sep="\t",
    chunksize=chunk_size,
    low_memory=True,
):
    filtered = chunk[chunk["sample"].isin(samples)]

    filtered.to_csv(
        f"{study}.csv",
        mode="w" if first_chunk else "a",
        header=first_chunk,
        index=False,
    )

    first_chunk = False
