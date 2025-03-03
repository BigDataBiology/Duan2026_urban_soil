def get_defense(infile1,infile2,outfile):
    with open(infile1,'rt') as f:
        defense = set()
        for line in f:
            defense.add(line.strip())

    with open(outfile,'wt') as out:
        with open(infile2,'rt') as f:
            for line in f:
                smorf,upstream,downstream = line.strip().split('\t')
                if upstream != 'NA # NA':
                    up_gene_defense_list = []
                    up_genes = upstream.split(' *** ')
                    for gene in up_genes:
                        up_defense = set()
                        parts = gene.split(' # ')
                        gene_name = f'{parts[0]} # {parts[1]} # {parts[2]} # {parts[3]}'
                        eggnog = gene.split(' # ')[-1]
                        if eggnog.startswith('NA'):
                            continue
                        else:
                            cogs = eggnog.split('@')[0]
                            pfams = eggnog.split('@')[2]
                            cogs_list = cogs.split(',')
                            pfams_list = pfams.split(',')
                            for cog in cogs_list:
                                if cog in defense:
                                    up_defense.add(cog)
                            for pfam in pfams_list:
                                if pfam in defense:
                                    up_defense.add(pfam)
                            up_defenses = ','.join(list(up_defense))
                            if up_defenses != '':
                                up_gene_defense = f'{gene_name} # {up_defenses}'
                                up_gene_defense_list.append(up_gene_defense)
                            else:
                                continue
                    if up_gene_defense_list != []:
                        up_gene_defenses = ' *** '.join(up_gene_defense_list)
                    else:
                        up_gene_defenses = 'NA'
                else:
                    up_gene_defenses = 'NA'
                if downstream != 'NA # NA':
                    down_genes = downstream.split(' *** ')
                    down_gene_defense_list = []
                    for gene in down_genes:
                        down_defense = set()
                        parts = gene.split(' # ')
                        gene_name = f'{parts[0]} # {parts[1]} # {parts[2]} # {parts[3]}'
                        eggnog = gene.split(' # ')[-1]
                        if eggnog.startswith('NA'):
                            continue
                        else:
                            cogs = eggnog.split('@')[0]
                            pfams = eggnog.split('@')[2]
                            cogs_list = cogs.split(',')
                            pfams_list = pfams.split(',')
                            for cog in cogs_list:
                                if cog in defense:
                                    down_defense.add(cog)
                            for pfam in pfams_list:
                                if pfam in defense:
                                    down_defense.add(pfam)
                            down_defenses = ','.join(list(down_defense))
                            if down_defenses != '':
                                down_gene_defense = f'{gene_name} # {down_defenses}'
                                down_gene_defense_list.append(down_gene_defense)
                            else:
                                continue
                    if down_gene_defense_list != []:
                        down_gene_defenses = ' *** '.join(down_gene_defense_list)
                    else:
                        down_gene_defenses = 'NA'
                else:
                    down_gene_defenses = 'NA'
                out.write(f'{smorf}\t{up_gene_defenses}\t{down_gene_defenses}\n')

def filter_defense(infile,outfile):
    with open(outfile,'wt') as out:
        with open(infile,'rt') as f:
            for line in f:
                smorf,upstream,downstream = line.strip().split('\t')
                if upstream != 'NA' and downstream != 'NA':
                    out.write(line)

infile1 = 'defense_id.tsv'
infile2 = 'smorf_genes_10_eggnog_new.tsv'
outfile1 = 'smorf_genes_10_eggnog_defense.tsv'
get_defense(infile1,infile2,outfile1)

outfile2 = 'smorf_genes_10_eggnog_defense_filter.tsv'
filter_defense(outfile1,outfile2)