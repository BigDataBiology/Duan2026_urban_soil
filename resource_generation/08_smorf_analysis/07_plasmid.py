def plasmid(infile0,infile1,infile2,outfile):
    tax = {}
    with open(infile0,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            tax[linelist[0]] = linelist[2]

    circ = set()
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            if linelist[3] == 'Y':
                circ.add(linelist[0])
    
    out = open(outfile,'wt')
    with open(infile2,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            parts = linelist[0].split('_')
            contig = f'{parts[0]}_{parts[1]}_{parts[2]}'
            if float(linelist[5]) >0.95:
                if contig in circ:
                    out.write(f'{contig}\tcircular\t{linelist[1]}\t{linelist[5]}\t{tax[contig]}\n')
                else:
                    out.write(f'{contig}\tlinear\t{linelist[1]}\t{linelist[5]}\t{tax[contig]}\n')
    out.close()
                
infile0 = 'contigs_mags_taxonomy.tsv'
infile1 = 'assembly_info.tsv'
infile2 = 'contigs_plasmid_summary.tsv'
outfile = 'contigs_plasmid_all.tsv'
plasmid(infile0,infile1,infile2,outfile)

def contig_cpnumber(infile,outfile):
    cp = {}
    with open(infile,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            smorfs = linelist[4].split(',')
            for item in smorfs:
                parts = item.split(' # ')[0].split('_')
                contig = f'{parts[0]}_{parts[1]}_{parts[2]}'
                if contig not in cp.keys():
                    cp[contig] = 0
                cp[contig] += 1
    with open(outfile,'wt') as out:
        for key,item in cp.items():
            out.write(f'{key}\t{item}\n')

infile = 'smorf_mapped_contig.tsv'
outfile = 'contig_cpnumber.tsv'
contig_cpnumber(infile,outfile)

def map_cp(infile1,infile2,outfile):
    cp = {}
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            cp[linelist[0]] = linelist[1]

    with open(outfile,'wt') as out:
        with open(infile2,'rt') as f:
            for line in f:
                linelist = line.strip().split('\t')
                sample = linelist[0].split('_')[0]
                if linelist[0] in cp.keys():
                    out.write(f'{line.strip()}\t{sample}\t{cp[linelist[0]]}\n')
                else:
                    out.write(f'{line.strip()}\t{sample}\t0\n')

infile1 = 'contig_cpnumber.tsv'
infile2 = 'contigs_plasmid_all.tsv'
outfile = 'contigs_plasmid_all_cpnumber.tsv'
map_cp(infile1,infile2,outfile)

def select(infile0,infile1,infile2,outfile):
    import random
    plasmid = set()
    with open(infile0,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            plasmid.add(linelist[0])

    length = {}
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            if linelist[0] not in plasmid:
                length[linelist[0]] = linelist[1]

    cp = {}
    with open(infile2,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            cp[linelist[0]] = linelist[1]

    with open(outfile,'wt') as out:
        selected_items = random.sample(list(length.items()), 13695)
        for item in selected_items:
            sample = item[0].split('_')[0]
            if item[0] in cp.keys():
                out.write(f'{item[0]}\t{item[1]}\t{cp[item[0]]}\t{sample}\n')
            else:
                out.write(f'{item[0]}\t{item[1]}\t0\t{sample}\n')

infile0 = 'contigs_plasmid_all_cpnumber.tsv'
infile1 = 'assembly_info.tsv'
infile2 = 'contig_cpnumber.tsv'
outfile = 'contig_selected.tsv'
select(infile0,infile1,infile2,outfile)

def assign_smorf(infile1,infile2,outfile):
    plasmid = set()
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            plasmid.add(linelist[0])
    plasmid_smorf = set()
    with open(infile2,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            contigs = linelist[3].split(',')
            for item in contigs:
                if item in plasmid:
                    plasmid_smorf.add(linelist[0])
    with open(outfile,'wt') as out:
        for item in plasmid_smorf:
            out.write(f'{item}\n')

infile1 = 'D:\soil\mge\contigs_plasmid_all_cpnumber.tsv'
infile2 = 'smorf_mapped_contig.tsv'
outfile = 'smorf_plasmid.tsv'
assign_smorf(infile1,infile2,outfile)

def assign_cluster(infile1,infile2,outfile):
    smorf = set()
    with open(infile1,'rt') as f:
        for line in f:
            smorf.add(line.strip())
    cluster = set()
    with open(infile2,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            smorfs = linelist[3].split(',')
            for item in smorfs:
                if item in smorf:
                    cluster.add(linelist[0])
    with open(outfile,'wt') as out:
        for item in cluster:
            out.write(f'{item}\n')

infile1 = 'smorf_plasmid.tsv'
infile2 = 'cluster_count_newname.tsv'
outfile = 'cluster_plasmid.tsv'
assign_cluster(infile1,infile2,outfile)

def fraction(infile1,infile2,outfile):
    plasmid = set()
    with open(infile1,'rt') as f:
        for line in f:
            plasmid.add(line.strip())
    cluster = {}
    with open(infile2,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            if linelist[4] not in cluster.keys():
                cluster[linelist[4]] = [linelist[0]]
            else:
                cluster[linelist[4]].append(linelist[0])
    with open(outfile,'wt') as out:
        for key,value in cluster.items():
            p = 0
            for item in value:
                if item in plasmid:
                    p += 1
            fraction = p/len(value)
            out.write(f'{key}\t{fraction}\n')

infile1 = 'cluster_plasmid.tsv'
infile2 = 'cluster_habitat_general.tsv'
outfile = 'habitat_plasmid_fraction_general.tsv'
fraction(infile1,infile2,outfile)

def split(infile1,infile2,outfile):
    single = set()
    with open(infile1,'rt') as f:
        for line in f:
            cluster,n = line.strip().split('\t')
            if int(n) < 2:
                single.add(cluster)
    out = open(outfile,'wt')

    with open(infile2,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            if linelist[0] in single:
                out.write(f'{line.strip()}\tsingle sample\n')
            else:
                out.write(f'{line.strip()}\tmultiple samples\n')

    out.close()

infile1 = 'cluster_sample.tsv'
infile2 = 'cluster_habitat_general.tsv'
outfile = 'cluster_habitat_general_sample.tsv'
split(infile1,infile2,outfile)

def cal_cdd(infile1,infile2,outfile):
    plasmid = set()
    with open(infile1,'rt') as f:
        for line in f:
            plasmid.add(line.strip())
    with open(outfile,'wt') as out:
        with open(infile2,'rt') as f:
            for line in f:
                linelist = line.strip().split('\t')
                if linelist[0] in plasmid:
                    out.write(line)

infile1 = 'cluster_plasmid.tsv'
infile2 = 'all.macrel.mapped.smorfs_dedup_cdd_tcov_family_newname.tsv'
outfile = 'cluster_plasmid_cdd.tsv'
cal_cdd(infile1,infile2,outfile)

def add_info(infile):
    import pandas as pd
    result = pd.read_csv(infile,sep='\t',header=None,names=['cluster','smorf','cdd','query_length','score','align_length','identity','evalue','target_length','tcov'])
    df = result.smorf.value_counts()
    bg_count = df.shape[0]
    bg_multi_count = df[df>1].shape[0]
    bg_perc = bg_count/14448
    bg_multi_perc = bg_multi_count/14448

    print(f'{bg_perc}({bg_count}/14448) of plasmid 90AA smORFs are annotated with CDD.')
    print(f'{bg_multi_perc}({bg_multi_count}/1448) of plasmid 90AA smORFs are annotated with multiple CDD.')

infile = 'cluster_plasmid_cdd.tsv'
add_info(infile)