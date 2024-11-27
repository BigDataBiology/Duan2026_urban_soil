def add_sample(infile,outfile,name):
    with open(outfile,'wt') as out:
        with open(infile,'rt') as f:
            for line in f:
                if line.startswith('seq'):
                    continue
                else:
                    out.write(f'{name}_{line}')
name = 's30'
fix = 'sample30'
infile = f'{name}_medaka_polypolish.fasta.PolcaCorrected_plasmid_summary.tsv'
outfile = f'{name}_plasmid_summary.tsv'
add_sample(infile,outfile,fix)

def map_mge(infile1,infile2,outfile):
    mge = {}
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            info = f'{linelist[1]}\t{linelist[5]}\t{linelist[7]}\t{linelist[9]}\t{linelist[10]}'
            mge[linelist[0]] = info
    
    with open(outfile,'wt') as out:
        with open(infile2,'rt') as f:
            for line in f:
                contigpart = line.split('\t')[0].split(' # ')[0].split('_')
                if len(contigpart) == 5:
                    contig = f'{contigpart[0]}_{contigpart[1]}_{contigpart[2]}_{contigpart[3]}'
                else:
                    contig = f'{contigpart[0]}_{contigpart[1]}_{contigpart[2]}'
                line = line.replace('\n','\t')
                if contig in mge.keys():
                    out.write(f'{line}Y\t{mge[contig]}\n')
                else:
                    out.write(f'{line}N\n')
                

infile1 = 'contigs_plasmid_summary.tsv'
infile2 = 'contig_rgi_pf_st_dedup.tsv'
outfile = 'contig_rgi_pf_st_dedup_mge_info.tsv'
map_mge(infile1,infile2,outfile)