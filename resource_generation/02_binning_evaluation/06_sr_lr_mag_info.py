def get_info(infile1,infile2,outfile,sample):
    mq = set()
    with open(infile1,'rt') as f:
        for line in f:
            mq.add(line.strip())
    with open(outfile,'wt') as out:
        with open(infile2,'rt') as f:
            for line in f:
                if line.startswith('filename'):
                    continue
                linelist = line.strip().split('\t')
                mag = linelist[0].split('/')[-1].split('.')[0]
                if mag in mq:
                    line = line.replace('semibin_output/output_bins/','').replace('.fa.gz','')
                    out.write(f'{sample}_{line}')

##for i in range(1,5):
sample = f'sample26'
infile1 = f'D:/soil/binning/sr/sample26_checkm/pass_checkm.txt'
infile2 = f'D:/soil/binning/sr/sample26_checkm/recluster_bins_info.tsv'
outfile = fr'D:\soil\binning\sr/{sample}_bins_info.tsv'
get_info(infile1,infile2,outfile,sample)

def cal_count(infile):
    s16 = 0
    s23 = 0
    s5 = 0
    with open(infile,'rt') as f:
        for line in f:
            if line.startswith('Sample'):
                continue
            else:
                line = line.rstrip()
                linelist= line.split(',')
                if linelist[11] != '0':
                    s16 +=1 
                if linelist[12] != '0':  
                    s23 +=1     
                if linelist[13] != '0':  
                    s5 +=1     
    s16 = s16/7949
    s23 = s23/7949
    s5 = s5/7949
    print(f'16S: {s16}\n23S: {s23}\n5S: {s5}')

infile1 = r'D:\soil\barrnap\sr\all_barrnap.tsv'
infile2 =  r'D:\soil\gtdb\r226\mag_gtdb_r226_quality_S16.csv'
cal_count(infile1)
cal_count(infile2)