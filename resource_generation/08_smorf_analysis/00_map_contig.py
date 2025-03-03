def store_length(infile1):
    length_dict = {}
    with open(infile1,'rt') as f:
        for line in f:
            contig,length = line.strip().split('\t')
            parts = contig.split('_')
            contig_name = f'{parts[0]}_{parts[1]}_{parts[2]}'
            length_dict[contig_name] = length
    return length_dict

def map_contig(length_dict,infile2,infile3,outfile):
    from fasta import fasta_iter
    dedup_smorf = set()
    for h,seq in fasta_iter(infile2):
        dedup_smorf.add(seq)

    smorf_contig = {}
    smorf_header = {}
    for h,seq in fasta_iter(infile3):
        parts = h.split('_')
        contig = f'{parts[0]}_{parts[1]}_{parts[2]}'
        #if seq not in smorf_contig.keys():
        if seq in dedup_smorf:
            if seq not in smorf_contig.keys():
                smorf_contig[seq] = set()
                smorf_header[seq] = []
            smorf_contig[seq].add(contig)
            smorf_header[seq].append(h)
    smorf_longest = {}
    for key,value in smorf_contig.items():
        pairs = [(contig_name,int(length_dict[contig_name])) for contig_name in value]
        smorf_longest[key] = max(pairs, key=lambda pair: pair[1])[0]

    with open(outfile,'wt') as out:
        n = 0
        for key,value in smorf_longest.items():
            nf = f'{n:07}'
            name = f'{nf[:1]}_{nf[1:4]}_{nf[4:7]}'
            contigs_list = list(smorf_contig[key])
            contigs = ','.join(contigs_list)
            headers = ','.join(smorf_header[key])
            out.write(f'{name}\t{key}\t{value}\t{contigs}\t{headers}\n')
            n += 1

infile1 = '/data/yiqian/soil/pipeline/02_polish/length/contig_length.tsv'
infile2 = '/data/yiqian/soil/pipeline/10_smorf/combine_smorf/all_macrel.mapped.smorfs.faa'
infile3 = '/data/yiqian/soil/pipeline/10_smorf/combine_all_smorf/all_macrel.out.all_smorfs_sample.faa'
outfile = 'smorf_mapped_contig.tsv'
length_dict= store_length(infile1)
map_contig(length_dict,infile2,infile3,outfile)