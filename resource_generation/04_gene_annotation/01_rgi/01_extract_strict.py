'''
extract ARGs annotated as strict.
add MAGs where they are from
'''
def mags(infile1,infile2,outfile):
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
                    if linelist[5] == 'Strict':
                        parts = linelist[0].split(' # ')[0].split('_')
                        if len(parts) == 5:
                            contig_name = f'{parts[1]}_{parts[2]}_{parts[3]}'
                        else:
                            contig_name = f'{parts[1]}_{parts[2]}'
                        out.write(f'sample1_SemiBin_{mags[contig_name]}\t{line}')

infile1 = 'contig_bins.tsv'
infile2 = 'sample1_mags_rgi.tsv'
outfile = 'sample1_mags_rgi_strict.tsv'
mags(infile1,infile2,outfile)