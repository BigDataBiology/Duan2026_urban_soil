sample = 'sample1'
with open(f'{sample}_tax.tsv','wt') as out:
    with open(f'{sample}_gtdb.tsv','rt') as f:
        for line in f:
            if line.startswith('user_genome'):
                continue
            else:
                linelist = line.strip().split('\t')
                out.write(f'{sample}_{linelist[0]}\t{linelist[1]}\n')