import glob
import pandas as pd

def read_and_set_index(n_f):
    n, f = n_f[0], n_f[1]
    return pd.read_table(f, sep='\t').set_index('taxonomy')


INPUT_FILE = glob.glob(r'~\soil\profile\result\*.tsv') 
OUTPUT_FILE = r"~\soil\profile\singlem_result.tsv"

dfs = map(read_and_set_index, enumerate(INPUT_FILE))
pd.concat(dfs, axis=1, sort=True, join='outer').to_csv(OUTPUT_FILE, sep='\t')