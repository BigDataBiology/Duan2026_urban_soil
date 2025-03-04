def filter_eggnog(infile,outfile):
    with open(outfile,'wt') as out:
        with open(infile,'rt') as f:
            for line in f:
                linelist = line.strip().split('\t')
                if float(linelist[2]) <0.00001:
                    out.write(line)

infile = 'all.emapper.annotations.tsv'
outfile = 'all.emapper.annotations.filter.tsv'
filter_eggnog(infile,outfile)

def map_eggnog(infile1,infile2,outfile1):
    eggnog_dict = {}
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            og = linelist[4].split(',')
            cogs = ','.join([item.split('@')[0] for item in og])
            info = linelist[7]
            pfams = linelist[20]
            eggnog = f'{cogs}@{info}@{pfams}'
            eggnog_dict[linelist[0]] = eggnog
        eggnog_dict['NA'] = 'NA'

    with open(outfile1,'wt') as out1:
        with open(infile2,'rt') as f:
            for line in f:
                linelist = line.strip().split('\t')
                upstream_list = linelist[1].split(',')
                downstream_list = linelist[2].split(',')
                up = []
                for item in upstream_list:
                    gene = item.split(' # ')[0]
                    if gene in eggnog_dict.keys():
                        eggnog = eggnog_dict[gene]
                    else:
                        eggnog = 'NA'
                    upstream = f'{item} # {eggnog}'
                    up.append(upstream)

                final_up = ' *** '.join(up)

                down = []
                for item in downstream_list:
                    gene = item.split(' # ')[0]
                    if gene in eggnog_dict.keys():
                        eggnog = eggnog_dict[gene]
                    else:
                        eggnog = 'NA'
                    downstream = f'{item} # {eggnog}'
                    down.append(downstream)

                final_down = ' *** '.join(down)

                out1.write(f'{linelist[0]}\t{final_up}\t{final_down}\n')

infile1 = 'all.emapper.annotations.filter.tsv'
infile2 = 'smorf_genes_10.tsv'
outfile1 = 'smorf_genes_10_eggnog.tsv'
map_eggnog(infile1,infile2,outfile1) 