cogs = {}
def cal_cog(infile,outfile):
    with open(infile,'rt') as f:
        for line in f:
            coglist = line.split('\t')[8]
            for i in coglist:
                if i not in cogs.keys():
                    cogs[i] = 1
                else:
                    cogs[i] += 1
    with open(outfile,'wt') as out:
        for key,value in cogs.items():
            out.write(f'{key}\t{value}\n')

infile = 'mags_dedup_eggnog.tsv'
#infile = 'contig_dedup_eggnog.tsv'
outfile = 'cog_number.tsv'
cal_cog(infile,outfile)