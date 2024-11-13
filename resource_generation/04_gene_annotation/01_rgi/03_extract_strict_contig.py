def filter_strict(infile,outfile,name):
    with open(outfile,'wt') as out:
        #out.write(f'ORF_ID	Contig	Start	Stop	Orientation	Cut_Off	Pass_Bitscore	Best_Hit_Bitscore	Best_Hit_ARO	Best_Identities	ARO	Model_type	SNPs_in_Best_Hit_ARO	Other_SNPs	Drug Class	Resistance Mechanism	AMR Gene Family	Predicted_DNA	Predicted_Protein	CARD_Protein_Sequence	Percentage Length of Reference Sequence	ID	Model_ID	Nudged	Note	Hit_Start	Hit_End	Antibiotic\n')
        with open(infile,'rt') as f:
            for line in f:
                linelist = line.strip().split('\t')
                if linelist[5] == 'Strict':
                    out.write(f'{name}_{line}')

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

def cal_aro_number(infile,outfile):
    sample_aro = {}
    all_aro = set()
    sh_aro = set()
    nj_aro = set()

    with open(infile,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            sample = linelist[0].split('_')[0]
            aro = linelist[10]
            if sample not in sample_aro.keys():
                sample_aro[sample] = set()
            else:
                sample_aro[sample].add(aro)
            all_aro.add(aro)
            if sample.startswith('sample'):
                sh_aro.add(aro)
            else:
                nj_aro.add(aro)
    
    print(len(all_aro))
    intersection = sh_aro&nj_aro
    print(f'SH:{len(sh_aro)}\tNJ:{len(nj_aro)}\tintersection:{len(intersection)}')

    with open(outfile,'wt') as out:
        for key,value in sample_aro.items():
            number = len(value)
            out.write(f'{key}\t{number}\n')

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
            if linelist[10] not in aro:
                aro[linelist[10]] = [1,linelist[8]]
            else:
                aro[linelist[10]][0] += 1
            if linelist[14] not in drug:
                drug[linelist[14]] = 1
            else:
                drug[linelist[14]] += 1
            if linelist[15] not in mechanism:
                mechanism[linelist[15]] = 1
            else:
                mechanism[linelist[15]] += 1
            if linelist[16] not in family:
                family[linelist[16]] = 1
            else:
                family[linelist[16]] += 1
    for key,value in aro.items():
        out1.write(f'ARO:{key} {value[1]}\t{value[0]}\n')
    for key,value in drug.items():
        out2.write(f'{key}\t{value}\n')
    for key,value in mechanism.items():
        out3.write(f'{key}\t{value}\n')
    for key,value in family.items():
        out4.write(f'{key}\t{value}\n')

name = 'snj01' 
infile1 = f'{name}_contig.rgi.txt'
outfile1 = f'{name}_contig_rgi_strict.tsv'
filter_strict(infile1,outfile1,name)

infile2 = 'contig_rgi_strict.tsv'
infile3 = 'contigs_dedup.faa'
outfile2 = 'contig_rgi_strict_dedup.tsv'
map_gene(infile2,infile3,outfile2)

outfile3 = 'contig_sample_aro_number.tsv'
cal_aro_number(outfile2,outfile3)

outfile4 = 'contig_aro_number.tsv'
outfile5 = 'contig_drug_number.tsv'
outfile6 = 'contig_mechanism_number.tsv'
outfile7 = 'contig_family_number.tsv'
cal(outfile2,outfile4,outfile5,outfile6,outfile7)