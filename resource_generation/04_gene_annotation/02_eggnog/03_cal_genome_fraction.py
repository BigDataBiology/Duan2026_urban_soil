def cal(infile1,infile2,infile3,outfile,name):
    from fasta import fasta_iter

    sample_annot = {}
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            if linelist[0] not in sample_annot:
                sample_annot[linelist[0]] = 1
            else:
                sample_annot[linelist[0]] += 1

    mags = {}
    sample_gene = {}
    with open(infile2,'rt') as f:
        for line in f:
            contig,bin = line.strip().split('\t')
            mags[contig] = bin

    with open(outfile,'wt') as out:
        for h,seq in fasta_iter(infile3):
            parts = h.split('_')
            if len(parts) == 4:
                contig_name = f'{parts[0]}_{parts[1]}_{parts[2]}'
            else:
                contig_name = f'{parts[0]}_{parts[1]}'
            mag_name = f'{name}_SemiBin_{mags[contig_name]}'
            if mag_name not in sample_gene:
                sample_gene[mag_name] = 1
            else:
                sample_gene[mag_name] += 1
        for key,value in sample_gene.items():
            out.write(f'{key}\t{value}\t{sample_annot[key]}\n')

name = 'sample1'
infile1 = f'{name}.mag.emapper.annotations.tsv'
infile2 = f'~/08_gene_prediction/{name}/contig/contig_bins.tsv'
infile3 = f'~/08_gene_prediction/mags/{name}_mag.faa'
outfile = f'{name}_eggnog_genome_fraction.tsv'
cal(infile1,infile2,infile3,outfile,name)