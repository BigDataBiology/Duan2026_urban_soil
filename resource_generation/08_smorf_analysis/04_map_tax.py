def filter_mq(infile1,infile2,outfile):
    mq_set = set()
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            mq = linelist[1]
            mq_set.add(mq)

    with open(outfile,'wt') as out:
        with open(infile2,'rt') as f:
            for line in f:
                contig,bins = line.strip().split('\t')
                if bins in mq_set:
                    out.write(f'{contig}\t{bins}\n')

infile1 = 'MQ_mags.tsv'
infile2 = 'all_contig_bins.tsv'
outfile = 'all_contig_bins_mq.tsv'
filter_mq(infile1,infile2,outfile)

def map_tax(infile1,infile2,infile3,outfile):
    contig_mq = {}
    with open(infile1,'rt') as f:
        for line in f:
            contig,bins = line.strip().split('\t')
            contig_mq[contig] = bins

    mags_tax = {}
    with open(infile2,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            mags_tax[linelist[1]] = linelist[16]

    with open(outfile,'wt') as out:
        with open(infile3,'rt') as f:
            for line in f:
                linelist = line.strip().split('\t')
                contig = linelist[0].replace('_polypolish','')
                if contig in contig_mq.keys():
                    mag = contig_mq[contig]
                    if mags_tax[mag] != '':
                        tax = mags_tax[mag]
                    else:
                        tax = 'NA'
                    
                    if tax=='NA':
                        rank = 'NA'
                    elif tax.endswith('p__;c__;o__;f__;g__;s__') or tax=='Unclassified Bacteria':
                        rank = 'superkingdom'
                    elif tax.endswith('c__;o__;f__;g__;s__'):
                        rank = 'phylum'
                    elif tax.endswith('o__;f__;g__;s__'):
                        rank = 'class'
                    elif tax.endswith('f__;g__;s__'):
                        rank = 'order'
                    elif tax.endswith('g__;s__'):
                        rank = 'family'
                    elif tax.endswith('s__'):
                        rank = 'genus'
                    else:
                        rank = 'species'

                    out.write(f'{contig}\t{mag}\t{tax}\t{rank}\n')
                else:
                    tax = f'{linelist[1]} # {linelist[2]} # {linelist[3]}'
                    if linelist[2] == 'no rank':
                        rank = 'NA'
                    else:
                        rank = linelist[2]
                    out.write(f'{contig}\tNA\t{tax}\t{rank}\n')

infile1 = 'all_contig_bins_mq.tsv'
infile2 = 'MQ_mags.tsv'
infile3 = 'contigs_taxonomy.tsv'
outfile = 'contigs_mags_taxonomy.tsv'
map_tax(infile1,infile2,infile3,outfile)

def smorf_tax(infile1,infile2,outfile):
    rank_dict = {'NA':0,'superkingdom':1,'phylum':2,'class':3,'order':4,'family':5,'genus':6,'species':7}
    rank_T_dict = {0:'NA',1:'superkingdom',2:'phylum',3:'class',4:'order',5:'family',6:'genus',7:'species'}
    contig_tax = {}
    mag_contig = set()
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            contig_tax[linelist[0]] = linelist[3]
            if linelist[1] != 'NA':
                mag_contig.add(linelist[0])

    with open(outfile,'wt') as out:
        with open(infile2,'rt') as f:
            for line in f:
                linelist = line.strip().split('\t')
                contig_list = linelist[3].split(',')
                originate_set = set()
                for item in contig_list:
                    if item in mag_contig:
                        originate_set.add('mag')
                    else:
                        originate_set.add('contig')
                originate = ','.join(sorted(list(originate_set)))

                rank_list = [str(rank_dict[contig_tax[item]]) for item in contig_list]
                final_rank = max(rank_list)
                rank = ','.join(rank_list)
                out.write(f'{linelist[0]}\t{originate}\t{rank}\t{rank_T_dict[int(final_rank)]}\n')

infile1 = 'contigs_mags_taxonomy.tsv'
infile2 = 'smorf_mapped_contig.tsv'
outfile = 'smorf_mapped_tax.tsv'
smorf_tax(infile1,infile2,outfile)

def cal(infile,outfile1,outfile2,outfile3):
    mag = 0
    contig = 0
    both = 0
    mag_superkingdom = 0
    mag_phylum = 0
    mag_class = 0
    mag_order = 0
    mag_family = 0
    mag_genus = 0
    mag_species = 0
    mag_norank = 0
    contig_superkingdom = 0
    contig_phylum = 0
    contig_class = 0
    contig_order = 0
    contig_family = 0
    contig_genus = 0
    contig_species = 0
    contig_norank = 0
    all_superkingdom = 0
    all_phylum = 0
    all_class = 0
    all_order = 0
    all_family = 0
    all_genus = 0
    all_species = 0
    all_norank = 0

    out1 = open(outfile1,'wt')
    out2 = open(outfile2,'wt')
    out3 = open(outfile3,'wt')
    with open(infile,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            if linelist[1] == 'mag':
                mag += 1
                if linelist[3] == 'superkingdom':
                    mag_superkingdom += 1
                elif linelist[3] == 'phylum':
                    mag_phylum += 1
                elif linelist[3] == 'class':
                    mag_class += 1
                elif linelist[3] == 'order':
                    mag_order += 1
                elif linelist[3] == 'family':
                    mag_family += 1
                elif linelist[3] == 'genus':
                    mag_genus += 1
                elif linelist[3] == 'species':
                    mag_species += 1
                else:
                    mag_norank += 1

            elif linelist[1] == 'contig':
                contig += 1
                if linelist[3] == 'superkingdom':
                    contig_superkingdom += 1
                elif linelist[3] == 'phylum':
                    contig_phylum += 1
                elif linelist[3] == 'class':
                    contig_class += 1
                elif linelist[3] == 'order':
                    contig_order += 1
                elif linelist[3] == 'family':
                    contig_family += 1
                elif linelist[3] == 'genus':
                    contig_genus += 1
                elif linelist[3] == 'species':
                    contig_species += 1
                else:
                    contig_norank += 1
            else:
                both += 1
            if linelist[3] == 'superkingdom':
                all_superkingdom += 1
            elif linelist[3] == 'phylum':
                all_phylum += 1
            elif linelist[3] == 'class':
                all_class += 1
            elif linelist[3] == 'order':
                all_order += 1
            elif linelist[3] == 'family':
                all_family += 1
            elif linelist[3] == 'genus':
                all_genus += 1
            elif linelist[3] == 'species':
                all_species += 1
            else:
                all_norank += 1    
    out1.write(f'Taxonomy\tSuper_kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tNo_rank\n')
    out1.write(f'mags\t{mag_superkingdom}\t{mag_phylum}\t{mag_class}\t{mag_order}\t{mag_family}\t{mag_genus}\t{mag_species}\t{mag_norank}\n')
    out2.write(f'Taxonomy\tSuper_kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tNo_rank\n')
    out2.write(f'contigs\t{contig_superkingdom}\t{contig_phylum}\t{contig_class}\t{contig_order}\t{contig_family}\t{contig_genus}\t{contig_species}\t{contig_norank}\n')
    out3.write(f'Taxonomy\tSuper_kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tNo_rank\n')
    out3.write(f'all\t{all_superkingdom}\t{all_phylum}\t{all_class}\t{all_order}\t{all_family}\t{all_genus}\t{all_species}\t{all_norank}\n')
    
    out1.close()
    out2.close()
    out3.close()

    print(f'mag:{mag}\tcontig:{contig}\tboth:{both}\n')
            

infile = 'smorf_mapped_tax.tsv'
outfile1 = 'mags_tax.tsv'
outfile2 = 'contigs_tax.tsv'
outfile3 = 'all_tax.tsv'
cal(infile,outfile1,outfile2,outfile3)

def source(infile1,infile2,outfile):
    smorf_source = {}
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            smorf_source[linelist[0]] = linelist[1].split(',')
    cluster_source = {}
    with open(infile2,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            member_list = linelist[3].split(',')
            for item in member_list:
                source_list = smorf_source[item]
                if linelist[0] not in cluster_source.keys():
                    cluster_source[linelist[0]] = set()
                for i in source_list:
                    cluster_source[linelist[0]].add(i)
    with open(outfile,'wt') as out:
        for key,value in cluster_source.items():
            source = ','.join(sorted(list(value)))
            out.write(f'{key}\t{source}\n')

infile1 = 'smorf_mapped_tax.tsv'
infile2 = 'cluster_count_newname.tsv'
outfile = 'cluster_source.tsv'
source(infile1,infile2,outfile)