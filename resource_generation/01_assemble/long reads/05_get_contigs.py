'''
Get same contigs from flye assemblies,medaka polished assemblies, and polypolish polished assemblies of MQ Mags
'''
def get_contigs(indir,assembly,outdir,ppl):
    import os
    from fasta import fasta_iter

    fna = {}
    for h,seq in fasta_iter(assembly):
        fna[h] = seq

    for dirname, subdirs, files in os.walk(indir):
        for filename in files:
            if filename.endswith('.fa'):
                mag = filename.split('.')[0]
                outfile = f'{outdir}{mag}.fa'
                with open(outfile,'wt') as out:
                    infile = os.path.join(dirname, filename)
                    for h,seq in fasta_iter(infile):
                        if ppl:
                            contig = h
                        else:
                            contig = h.replace('_polypolish','')
                        out.write(f'>{contig}\n{fna[contig]}\n')

indir = './final_bins'
assembly = './flye_result/assembly.fasta'
outdir1 = './flye_s1/'
get_contigs(indir,assembly,outdir1,0)

medaka = './medaka/final_result/1_medaka.fasta'
outdir2 = './medaka_s1/'
get_contigs(indir,medaka,outdir2,0)

polypolish = './polypolish/medaka_polypolish.fasta'
outdir3 = './polypolish_s1/'
get_contigs(indir,polypolish,outdir3,1)