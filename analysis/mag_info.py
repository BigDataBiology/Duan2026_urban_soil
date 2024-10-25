'''
map MAG info
'''
def extract_info(infile):
    with open(infile,'rt') as f:
        for line in f:
            if line.startswith('filename'):
                continue
            else:
                name,info = line.strip().split('\t',1)
                bin = name.split('/')[2].split('.')[0]
    return info,bin

def store_checkm(infile1):
    checkm = {}
    with open(infile1,'rt') as f:
        for line in f:
            name,info = line.strip().split(',',1)
            name = name.split('.')[0]
            checkm[name] = info
    return checkm

def add_checkm(checkm,name):
    com = checkm[name].split(',')[0]
    con = checkm[name].split(',')[1]
    if float(com) > 90 and float(con) <5:
        quality = 'High quality'
    else:
        quality = 'Medium quality'
    return com,con,quality

def store_drep(infile2,infile3):
    ref_95 = set()
    ref_99 = set()
    with open(infile2,'rt') as f:
        for line in f:
            linelist = line.strip().split(',',1)
            bin = linelist[0].split('.')[0]
            ref_95.add(bin)
    with open(infile3,'rt') as f:
        for line in f:
            linelist = line.strip().split(',',1)
            bin = linelist[0].split('.')[0]
            ref_99.add(bin)
    return ref_95,ref_99

def add_drep(ref_95,ref_99,name):
    if name in ref_95:
        drep_95='T'
    else:
        drep_95='F'
    if name in ref_99:
        drep_99='T'
    else:
        drep_99='F'
    return drep_95,drep_99

def store_barrnap(infile4):
    barrnap = {}
    with open(infile4,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t',1)
            barrnap[linelist[0]] = linelist[1]
    return barrnap

def store_trnascan(infile5):
    trnascan = {}
    with open(infile5,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            trnascan[linelist[0]] = linelist[1]
    return trnascan

def store_gtdb(infile6):
    gtdb = {}
    with open(infile6,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            gtdb[linelist[0]] = linelist[1]
    return gtdb

def add_gtdb(gtdb,name):
    if name in gtdb.keys():
        tax = gtdb[name]
    else:
        tax = 'NA'
    return tax

def store_sgb(infile7):
    sgb = {}
    with open(infile7,'rt') as f:
        for line in f:
            linelist = line.strip().split(',')
            mag = linelist[0].replace('.fa','')
            sgb[mag] = linelist[1]
    return sgb

def store_rgi(infile8):
    rgi = {}
    with open(infile8,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t',1)
            rgi[linelist[0]] = linelist[1]
    return rgi

def add_rgi(rgi,name):
    if name in rgi:
        aro = rgi[name]
    else:
        aro = f'NA\tNA\tNA\tNA'
    return aro

def combine_info(inputdir1,samplename,checkm,ref_95,ref_99,gtdb,rgi):
    infile = f'{inputdir1}/{samplename}_bins_info.tsv'
    info,bin = extract_info(infile)
    newname = f'{samplename}_{bin}'
    com,con,quality = add_checkm(checkm,newname)
    drep_95,drep_99 = add_drep(ref_95,ref_99,newname)
    tax = add_gtdb(gtdb,newname)
    aro = add_rgi(rgi,newname)
    return newname,info,com,con,quality,drep_95,drep_99,tax,aro
            
def add_bins_info(inputdir1,infile1,infile2,infile3,infile4,infile5,infile6,infile7,infile8,outfile1):
    snj = ['snj01','snj02','snj03','snj04','snj08','snj09','snj10','snj11','snj13','snj14','snj15','snj16','snj19','snj20']
    cpsnj = ['cpsnj01','cpsnj02','cpsnj03','cpsnj04','cpsnj08','cpsnj09','cpsnj10','cpsnj11','cpsnj13','cpsnj14','cpsnj15','cpsnj16','cpsnj19','cpsnj20']

    checkm = store_checkm(infile1)
    ref_95,ref_99 = store_drep(infile2,infile3)
    barrnap = store_barrnap(infile4)
    trnascan = store_trnascan(infile5)
    gtdb = store_gtdb(infile6)
    sgb = store_sgb(infile7)
    rgi = store_rgi(infile8)

    with open(outfile1,'wt') as out:
        out.write(f'Sample\tMag ID\tnbps\tNr contigs\tN50\tL50\tCompleteness\tContamination\tQuality\tDrep_95\tDrep_99\tS16\tS23\tS5\tUnique tRNA\tGTDB\tSGB\tARO\tDrug Class\tResistance Mechanism\tAMR Gene Family\n')
        for samplename in snj:
            newname,info,com,con,quality,drep_95,drep_99,tax,aro = combine_info(inputdir1,samplename,checkm,ref_95,ref_99,gtdb,rgi)
            out.write(f'{samplename}\t{newname}\t{info}\t{com}\t{con}\t{quality}\t{drep_95}\t{drep_99}\t{barrnap[newname]}\t{trnascan[newname]}\t{tax}\t{sgb[newname]}\t{aro}\n')
        for samplename in cpsnj:
            newname,info,com,con,quality,drep_95,drep_99,tax,aro = combine_info(inputdir1,samplename,checkm,ref_95,ref_99,gtdb,rgi)
            out.write(f'{samplename}\t{newname}\t{info}\t{com}\t{con}\t{quality}\t{drep_95}\t{drep_99}\t{barrnap[newname]}\t{trnascan[newname]}\t{tax}\t{sgb[newname]}\t{aro}\n')
        for i in range(1,31):
            samplename = f'sample{i}'
            newname,info,com,con,quality,drep_95,drep_99,tax,aro = combine_info(inputdir1,samplename,checkm,ref_95,ref_99,gtdb,rgi)
            out.write(f'{samplename}\t{newname}\t{info}\t{com}\t{con}\t{quality}\t{drep_95}\t{drep_99}\t{barrnap[newname]}\t{trnascan[newname]}\t{tax}\t{sgb[newname]}\t{aro}\n')

def mimag(infile,outfile):
    with open(outfile,'wt') as out:
        with open(infile,'rt') as f:
            for line in f:
                if line.startswith('Sample'):
                    out.write(f'{line.strip()}\tMIMAG\n')
                else:
                    linelist = line.strip().split('\t')
                    if linelist[8] == 'High quality' and int(linelist[11])>0 and int(linelist[12])>0 and int(linelist[13])>0 and int(linelist[15])>=18:
                        out.write(f'{line.strip()}\tT\n')
                    else:
                        out.write(f'{line.strip()}\tF\n')

def cal_aro(infile,outfile):
    sample_aro = {}
    with open(infile,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            sample = linelist[0]
            aro = linelist[17]
            arolist = aro.split(';')
            if sample not in sample_aro.keys():
                sample_aro[sample] = set()
            else:
                for item in arolist:
                    sample_aro[sample].add(item)

    with open(outfile,'wt') as out:
        for key,value in sample_aro.items():
            number = len(value)
            out.write(f'{key}\t{number}\n')
            
def sgb(infile,outfile):
    sgb_dict = {}
    with open(infile,'rt') as f:
        for line in f:
            linelist = line.strip().split(',')
            if line.startswith('genome'):
                continue
            else:
                if linelist[1] not in sgb_dict.keys():
                    sgb_dict[linelist[1]] = 1
                else:
                    sgb_dict[linelist[1]] += 1
    with open(outfile,'wt') as out:
        for key,value in sgb_dict.items():
            out.write(f'{key}\t{value}\n')

inputdir1 = '/data/Projects/urban_soil/data/UrbanSoilTables/all_mags_info/'
infile1 = './pre-calculated_data/all_mag_info/checkm.tsv'
infile2 = './pre-calculated_data/all_mag_info/Widb_95.csv'
infile3 = './pre-calculated_data/all_mag_info/Widb_99.csv'
infile4 = './pre-calculated_data/all_mag_info/barrnap.tsv'
infile5 = './pre-calculated_data/all_mag_info/trnascan.tsv'
infile6 = './pre-calculated_data/all_mag_info/gtdb.tsv'
infile7 = './pre-calculated_data/all_mag_info/Cdb_95.csv'
infile8 = './pre-calculated_data/all_mag_info/rgi_mags_strict_stats.tsv'
outfile1 = './pre-calculated_data/sample_bins_info.tsv'
outfile2 = './pre-calculated_data/sample_bins_info_mimag.tsv'

add_bins_info(inputdir1,infile1,infile2,infile3,infile4,infile5,infile6,infile7,infile8,outfile1)
mimag(outfile1,outfile2)

outfile3 = './pre-calculated_data/sample_aro_number.tsv'
cal_aro(outfile2,outfile3)

outfile4 = 'sgb.tsv'
sgb(infile7,outfile4)