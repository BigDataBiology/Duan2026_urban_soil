import polars as pl
import sys

shanghai = pl.read_csv("./shanghai_chao1.csv")
other = pl.read_csv(sys.argv[1])


print(
    f"MEAN CHAO1 per study:\nShanghai:\t{shanghai.mean()['chao1'][0]}\nOther:\t{other.mean()['chao1'][0]}"
)
