def add_length(infile1,infile2,outfile):
    cdd_dict = {}
    with open(infile1,'rt') as f1:
        for line in f1:
            linelist = line.strip().split('\t')
            cdd_dict[linelist[0]] = linelist[4]

    with open(outfile,'wt') as out:
        with open(infile2,'rt') as f2:
            for line in f2:
                line = line.strip()
                linelist = line.split('\t')
                pssm = linelist[1].split('|')[2]
                if cdd_dict[pssm] != '0':
                    out.write(f'{line}\t{cdd_dict[pssm]}\n')

def filter_cov(infile,outfile):
    import pandas as pd
    result = pd.read_csv(infile,sep='\t',header=None,names=['smorf','cdd','query_length','score','align_length','identity','evalue','target_length'])
    result['tcov'] = result['align_length']/result['target_length']
    result = result[result['tcov'] >0.8]
    result.to_csv(outfile,sep='\t',index=None)

def filter_family(infile3,outfile2,outfile3):
    cluster = set()
    with open(infile3,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            cluster.add(linelist[2])
    with open(outfile3,'wt') as out:
        with open(outfile2,'rt') as f:
            for line in f:
                linelist = line.strip().split('\t')
                if linelist[0] in cluster:
                    out.write(line)

def add_info(infile,outfile):
    import pandas as pd
    result = pd.read_csv(infile,sep='\t',header=None,names=['smorf','cdd','query_length','score','align_length','identity','evalue','target_length','tcov'])
    df = result.smorf.value_counts()
    bg_count = df.shape[0]
    bg_multi_count = df[df>1].shape[0]
    bg_perc = bg_count/2206902
    bg_multi_perc = bg_multi_count/2206902

    print(f'{bg_perc}({bg_count}/2206902) of 90AA smORFs are annotated with CDD.')
    print(f'{bg_multi_perc}({bg_multi_count}/2206902) of 90AA smORFs are annotated with multiple CDD.')

infile1 = 'cddid_all.tbl'
infile2 = 'all.macrel.mapped.smorfs_dedup.tsv'
infile3 = 'cluster_count.tsv'
outfile1 = 'all.macrel.mapped.smorfs_dedup_cdd.tsv'
outfile2 = 'all.macrel.mapped.smorfs_dedup_cdd_tcov.tsv'
outfile3 = 'all.macrel.mapped.smorfs_dedup_cdd_tcov_family.tsv'
outfile4 = 'all.macrel.mapped.smorfs_dedup_cdd_tcov_family_bg.tsv'

add_length(infile1,infile2,outfile1)
filter_cov(outfile1,outfile2)
filter_family(infile3,outfile2,outfile3)
add_info(outfile3,outfile4)