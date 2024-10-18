def contig(infile0,infile1,infile2,outfile):
    from fasta import fasta_iter
    mags = set()
    with open(infile0,'rt') as f:
        for line in f:
            mag = line.strip().split('_')[2].split('.')[0]
            mags.add(mag)

    contigs = set()
    with open(infile1,'rt') as f:
        for line in f:
            contig,bin = line.strip().split('\t')
            if bin in mags:
                contigs.add(contig)

    with open(outfile,'wt') as out:
        for h,seq in fasta_iter(infile2):
            header = h.split('_')
            name = f'{header[0]}_{header[1]}_{header[2]}'
            if name not in contigs:
                out.write(f'>cpsnj01_{h}\n{seq}\n')

infile0 = 'mags.tsv'
infile1 = 'contig_bins.tsv'
infile2 = 'cpsnj01_contig.faa.gz'
outfile = 'cpsnj01_only_contig.faa'
contig(infile0,infile1,infile2,outfile)