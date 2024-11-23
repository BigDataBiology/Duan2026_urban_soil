'''
extract ARGs annotated as strict.
add MAGs where they are from
'''
def strict(infile1,infile2,outfile):
    mags = {}
    with open(infile1,'rt') as f:
        for line in f:
            contig,bin = line.strip().split('\t')
            mags[contig] = bin

    with open(infile2,'rt') as f:
        with open(outfile,'wt') as out:
            for line in f:
                linelist = line.strip().split('\t')
                if line.startswith('ORF'):
                    out.write(f'MAG\t{line}')
                else:
                    if linelist[5] == 'Strict' or linelist[5] == 'Perfect':
                        parts = linelist[0].split(' # ')[0].split('_')
                        if len(parts) == 5:
                            contig_name = f'{parts[1]}_{parts[2]}_{parts[3]}'
                        else:
                            contig_name = f'{parts[1]}_{parts[2]}'
                        out.write(f'sample1_SemiBin_{mags[contig_name]}\t{line}')

def stats(infile,outfile):
    mags = {}
    with open(infile,'rt') as f:
        for line in f:
            if line.startswith('MAG'):
                continue
            else:
                linelist = line.strip().split('\t')
                if linelist[0] not in mags.keys():
                    mags[linelist[0]] = [set(),set(),set(),set()]
                mags[linelist[0]][0].add(linelist[11])
                mags[linelist[0]][1].add(linelist[15])
                mags[linelist[0]][2].add(linelist[16])
                mags[linelist[0]][3].add(linelist[17])
    with open(outfile,'wt') as out:
        out.write(f'MAG\tARO\tDrug Class\tResistance Mechanism\tAMR Gene Family\n')
        for key,value in mags.items():
            aro = ';'.join(list(value[0]))
            drug = ','.join(list(value[1]))
            mechanism = ','.join(list(value[2]))
            family = ','.join(list(value[3]))
            out.write(f'{key}\t{aro}\t{drug}\t{mechanism}\t{family}\n')

def map_gene(infile1,infile2,outfile):
    from fasta import fasta_iter
    strict = {}
    with open(infile1,'rt') as f:
        for line in f:
            if line.startswith('MAG'):
                continue
            else:
                linelist = line.strip().split('\t')
                strict[linelist[1]] = line

    with open(outfile,'wt') as out:
        for h,seq in fasta_iter(infile2):
            #name = h.replace('mag_')
            if h in strict.keys():
                out.write(f'mag\t{strict[h]}')

def cal(infile,outfile1,outfile2,outfile3,outfile4):
    out1 = open(outfile1,'wt')
    out2 = open(outfile2,'wt')
    out3 = open(outfile3,'wt')
    out4 = open(outfile4,'wt')

    aro = {}
    drug = {}
    mechanism = {}
    family = {}

    with open(infile,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            if linelist[12] not in aro:
                aro[linelist[12]] = [1,linelist[10]]
            else:
                aro[linelist[12]][0] += 1
            if linelist[16] not in drug:
                drug[linelist[16]] = 1
            else:
                drug[linelist[16]] += 1
            if linelist[17] not in mechanism:
                mechanism[linelist[17]] = 1
            else:
                mechanism[linelist[17]] += 1
            if linelist[18] not in family:
                family[linelist[18]] = 1
            else:
                family[linelist[18]] += 1
    for key,value in aro.items():
        out1.write(f'{key}\t{value[0]}\t{value[1]}\n')
    for key,value in drug.items():
        out2.write(f'{key}\t{value}\n')
    for key,value in mechanism.items():
        out3.write(f'{key}\t{value}\n')
    for key,value in family.items():
        out4.write(f'{key}\t{value}\n')

infile1 = 'contig_bins.tsv'
infile2 = 'sample1_mags_rgi.tsv'
outfile1 = 'sample1_mags_rgi_strict.tsv'
strict(infile1,infile2,outfile1)

infile3 = 'all_rgi_mags_strict.tsv'
outfile2 = 'rgi_mags_strict_stats.tsv'
stats(infile3,outfile2)

infile4 = 'mags_dedup.faa'
outfile3 = 'mags_dedup_rgi.tsv'
map_gene(infile3,infile4,outfile3)

outfile4 = 'aro_number.tsv'
outfile5 = 'drug_number.tsv'
outfile6 = 'mechanism_number.tsv'
outfile7 = 'family_number.tsv'
cal(outfile3,outfile4,outfile5,outfile6,outfile7)