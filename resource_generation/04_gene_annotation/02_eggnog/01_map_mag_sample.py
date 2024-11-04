def mags(infile1,infile2,outfile,name):
    mags = {}
    with open(infile1,'rt') as f:
        for line in f:
            contig,bin = line.strip().split('\t')
            mags[contig] = bin

    with open(infile2,'rt') as f:
        with open(outfile,'wt') as out:
            for line in f:
                linelist = line.strip().split('\t')
                if line.startswith('#'):
                    continue
                else:
                    parts = linelist[0].split('_')
                    if len(parts) == 4:
                        contig_name = f'{parts[0]}_{parts[1]}_{parts[2]}'
                    else:
                        contig_name = f'{parts[0]}_{parts[1]}'
                    out.write(f'{name}_SemiBin_{mags[contig_name]}\t{name}_{line}')

infile1 = '/data/Projects/urban_soil/data/UrbanSoilTables/all_mags_info/sample1_contig_bins.tsv'
infile2 = 'sample1.emapper.annotations'
outfile = 'sample1.mag.emapper.annotations.tsv'
name = 'sample1'
mags(infile1,infile2,outfile,name)