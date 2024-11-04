def map_gene(infile1,infile2,outfile):
    from fasta import fasta_iter
    eggnog = {}
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            eggnog[linelist[1]] = line

    with open(outfile,'wt') as out:
        for h,seq in fasta_iter(infile2):
            name = h.replace('mag_','').split(' # ')[0]
            if name in eggnog.keys():
                out.write(f'mag\t{eggnog[name]}')

infile1 = 'mags.emapper.annotations.tsv'
infile2 = 'mags_dedup.faa'
outfile = 'mags_dedup_eggnog.tsv'
#infile1 = 'contig.emapper.annotations.tsv'
#infile2 = 'contig_dedup.faa'
#outfile = 'contig_dedup_eggnog.tsv'
map_gene(infile1,infile2,outfile)