def get_sequences(sample,outfile):
    from glob import glob
    from fasta import fasta_iter
    fs = glob('~/barrnap/*.faa')
    with open(outfile,'wt') as out:
        for f in fs:
            for h,seq in fasta_iter(f):
                mag = f.split('/')[-1].split('.')[0]
                if h.startswith('16S_rRNA::'):
                    h = f'{sample}_{mag}_{h}'
                    out.write(f'>{h}\n{seq}\n')

#sample = 'cpsnj04'
#outfile = f'~/barrnap/{sample}_16S.fna'
#get_sequences(sample,outfile)

def filter_sgb(infile1,infile2,outfile):
    from fasta import fasta_iter
    hq = set()
    with open(infile1,'rt') as f:
        for line in f:
            hq.add(line.strip())
    with open(outfile,'wt') as out:
        for h,seq in fasta_iter(infile2):
            mag = h.split('#')[0]
            if mag in hq:
                out.write(f'>{h}\n{seq}\n')
       
#infile1 = 'sgb_hq.txt'
#infile2 = 'mag_16S.fna'
#outfile = 'mag_hq_16S.fna'

def map_tax(infile1,infile2,outfile):
    gtdbdict = {}
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split(',')
            gtdbdict[linelist[7]] = line.strip()
    with open(outfile,'wt') as out:
        with open(infile2,'rt') as f:
            for line in f:
                if line.startswith('MAG'):
                    continue
                else:
                    linelist = line.strip().split(',')
                    out.write(f'{gtdbdict[linelist[0]]},{linelist[1]},{linelist[2]},{linelist[3]}\n')


#infile1 = 'sgb_gtdb.csv'
#infile2 = 'microbe-atlas.csv'
#outfile = 'sgb_16s_gtdb.csv'
#map_tax(infile1,infile2,outfile)

def check_tax(infile):
    count = 0
    mag_dict = set()
    with open(infile,'rt') as f:
        for line in f:
            linelist = line.strip().split(',')
            otulist = linelist[11].split(';')
            if linelist[7] not in mag_dict:
                #if otulist[-1].startswith('90') or otulist[-1].startswith('96'):
                    #if linelist[6] == 's__':
                        #count += 1
                if linelist[6] == 's__' or linelist[6] == '':
                    count += 1
            mag_dict.add(linelist[7])
    print(count)

infile = 'sgb_16s_gtdb.csv'
check_tax(infile)

def map_ref(infile1,infile2,outfile):
    '''
    otu_mag = {}
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split(',')
            otu = linelist[-1]
            otu_mag[otu] = line.strip()
    with open(infile2,'rt') as f:
        with open(outfile,'wt') as out:
            for line in f:
                linelist = line.strip().split('\t')
                if linelist[0] in otu_mag.keys():
                    out.write(f'{otu_mag[linelist[0]]},{linelist[5]}\n')
'''
    otu_ref = {}
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            otu_ref[linelist[0]] = linelist[5]

    with open(infile2,'rt') as f:
        with open(outfile,'wt') as out:
            for line in f:
                linelist = line.strip().split(',')
                if linelist[-1] in otu_ref.keys():
                    out.write(f'{line.strip()},{otu_ref[linelist[-1]]}\n')
                else:
                    out.write(f'{line.strip()},NA\n')
infile1 = r'~\mapref-3.0\otus.info'
infile2 = r'sgb_16s_gtdb_nos.csv'
outfile = r'sgb_16s_gtdb_nos_ref.csv'
map_ref(infile1,infile2,outfile)