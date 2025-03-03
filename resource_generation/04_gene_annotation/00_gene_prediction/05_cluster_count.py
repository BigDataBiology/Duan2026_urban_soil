def cluster(infile,outfile):
    cluster_dict = {}
    with open(infile,'rt') as f:
        for line in f:
            if line.startswith('>'):
                cluster_name = line.strip().replace('>','')
                cluster_dict[cluster_name] = ['',[]]
            else:
                smorf = line.strip().split('\t')[1].split('>')[1].split('...')[0]
                if line.startswith('0'):
                    cluster_dict[cluster_name][0] = smorf
                    cluster_dict[cluster_name][1].append(smorf)
                else:
                    cluster_dict[cluster_name][1].append(smorf)
    with open(outfile,'wt') as out:
        for key,value in cluster_dict.items():
            seqs = ','.join(value[1])
            out.write(f'{key}\t{len(value[1])}\t{value[0]}\t{seqs}\n')

infile = 'all.macrel.mapped.smorfs_dedup_cdhit.faa.clstr'
outfile = 'cluster_count.tsv'
cluster(infile,outfile)