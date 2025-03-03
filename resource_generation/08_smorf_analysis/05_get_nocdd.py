def newname(infile1,infile2,infile3,outfile):
    from fasta import fasta_iter

    seq_name = {}
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            seq_name[linelist[1]] = linelist[0]

    name_seq = {}
    for h,seq in fasta_iter(infile2):
        name_seq[h] = seq

    with open(outfile,'wt') as out:
        with open(infile3,'rt') as f:
            linelist = line.strip().split('\t')
            rep = seq_name[name_seq[linelist[2]]]
            member_list = linelist[3].split(',')
            members = ','.join([seq_name[name_seq[item]]for item in member_list])
            out.write(f'{linelist[0]}\t{linelist[1]}\t{rep}\t{members}\n')

infile1 = 'smorf_mapped_contig.tsv'
infile2 = 'all.mapped.smorfs_dedup.faa'
infile3 = 'cluster_count.tsv'
outfile = 'cluster_count_newname.tsv'
newname(infile1,infile2,infile3,outfile)

def add_family(infile1,infile2,outfile):
    cluster = {}
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            member_list = linelist[3].split(',')
            for item in member_list:
                cluster[item] = linelist[0]
    
    with open(outfile,'wt') as out:
        with open(infile2,'rt') as f:
            for line in f:
                linelist = line.strip().split('\t')
                out.write(f'{cluster[linelist[0]]}\t{line}')

infile1 = 'cluster_count_newname.tsv'
infile2 = 'smorf_genes_10_eggnog_defense_filter.tsv'
outfile = 'smorf_genes_10_eggnog_defense_filter_family.tsv'
add_family(infile1,infile2,outfile)

def cdd_newname(infile1,infile2,outfile):
    cluster_name = {}
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            cluster_name[linelist[2]] = linelist[0]
    with open(outfile,'wt') as out:
        with open(infile2,'rt') as f:
            for line in f:
                linelist = line.strip().split('\t')
                out.write(f'{cluster_name[linelist[0]]}\t{line}')

infile1 = 'cluster_count.tsv'
infile2 = 'all.mapped.smorfs_dedup_cdd_tcov_family.tsv'
outfile = 'all.mapped.smorfs_dedup_cdd_tcov_family_newname.tsv'
cdd_newname(infile1,infile2,outfile)


def nocdd(infile1,infile2,outfile):
    cdd = set()
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            cdd.add(linelist[0])
    with open(outfile,'wt') as out:
        with open(infile2,'rt') as f:
            for line in f:
                linelist = line.strip().split('\t')
                if linelist[0] not in cdd:
                    out.write(f'{line}')

infile1 = 'all.mapped.smorfs_dedup_cdd_tcov_family_newname.tsv'
infile2 = 'smorf_genes_10_eggnog_defense_filter_family.tsv'
outfile = 'smorf_genes_10_eggnog_defense_filter_family_nocdd.tsv'
nocdd(infile1,infile2,outfile)

def cal_cluster(infile1,infile2,outfile):
    length = {}
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            length[linelist[0]] = linelist[1]
    defense = {}
    with open(infile2,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            if linelist[0] not in defense.keys():
                defense[linelist[0]] = set()
            defense[linelist[0]].add(linelist[1])
    with open(outfile,'wt') as out:
        for key,value in defense.items():
            defense_len = len(value)
            fraction = defense_len/int(length[key])
            out.write(f'{key}\t{length[key]}\t{defense_len}\t{fraction}\n')
    
infile1 = 'cluster_count_newname.tsv'
infile2 = 'smorf_genes_10_eggnog_defense_filter_family_nocdd.tsv'
outfile = 'cluster_defense.tsv'
cal_cluster(infile1,infile2,outfile)

def map_tax(infile1,infile2,outfile):
    tax_dict = {}
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            tax = f'{linelist[1]}\t{linelist[2]}'
            tax_dict[linelist[0]] = tax
    with open(outfile,'wt') as out:
        with open(infile2,'rt') as f:
            for line in f:
                linelist = line.strip().split('\t')
                parts = linelist[2].split('_')
                contig = f'{parts[0]}_{parts[1]}_{parts[2]}'
                out.write(f'{line.strip()}\t{tax_dict[contig]}\n')

infile1 = 'contigs_mags_taxonomy.tsv'
infile2 = 'smorf_genes_10_eggnog_defense_filter_family_nocdd.tsv'
outfile = 'smorf_genes_10_eggnog_defense_filter_family_nocdd_tax.tsv'
map_tax(infile1,infile2,outfile)