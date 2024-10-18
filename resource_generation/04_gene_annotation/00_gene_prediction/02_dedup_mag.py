def dedup(infile,outfile1,outfile2):
    from fasta import fasta_iter
    dedup = {}
    with open(outfile1,'wt') as out1:
        for h,seq in fasta_iter(infile):
            if seq not in dedup.keys():
                out1.write(f'>{h}\n{seq}\n')
                dedup[seq] = 1
            else:
                dedup[seq] += 1

    with open(outfile2,'wt') as out2:
        for key,value in dedup.items():
            out2.write(f'{key}\t{value}\n')

infile = 'mags.faa.gz'
outfile1 = 'mags_dedup.faa'
outfile2 = 'mags_cpnumber.tsv'
dedup(infile,outfile1,outfile2)