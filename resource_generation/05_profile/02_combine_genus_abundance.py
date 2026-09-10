import glob
import pandas as pd

def read_and_set_index(n_f):
    n, f = n_f[0], n_f[1]
    return pd.read_table(f, sep='\t').set_index('taxonomy')


INPUT_FILE = glob.glob(r'~\soil\profile\result\*.tsv') 
OUTPUT_FILE = r"genus_result.tsv"

dfs = map(read_and_set_index, enumerate(INPUT_FILE))
pd.concat(dfs, axis=1, sort=True, join='outer').to_csv(OUTPUT_FILE, sep='\t')

df = pd.read_csv(r'abundance_profile_genus.tsv', sep='\t',index_col='taxonomy')
df = df.fillna(0)
df.columns = ['CPSNJ01', 'CPSNJ02','CPSNJ03','CPSNJ04', 'CPSNJ08','CPSNJ09','CPSNJ10', 'CPSNJ11','CPSNJ13','CPSNJ14', 'CPSNJ15','CPSNJ16','CPSNJ19','CPSNJ20','s10','s11','s12','s13','s14','s15','s16','s17','s18','s19','s1','s20','s21','s22','s23','s24','s25','s26','s27','s28','s29','s2','s30','s3','s4','s5','s6','s7','s8','s9','SNJ01','SNJ02','SNJ03','SNJ04','SNJ08','SNJ09','SNJ10','SNJ11','SNJ13','SNJ14','SNJ15','SNJ16','SNJ19','SNJ20']
df['Mean'] = df.mean(axis=1)
df = df.loc[df['Mean'].between(0.00001, 1)]
df = df.drop(['Mean'], axis=1) 
df.columns = [col.split(";")[-1] for col in df.columns]
df.to_csv(r'genus_abundance_filter.tsv',sep='\t')

