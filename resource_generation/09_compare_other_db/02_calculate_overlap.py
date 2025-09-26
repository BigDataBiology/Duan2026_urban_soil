def cal_spire(infile1,infile2,outfile):
    sgb = set()
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            if linelist[9] == 'T':
                sgb.add(linelist[1])
    print(len(sgb))
    
    with open(outfile,'wt') as out:
        with open(infile2,'rt') as f:
            for line in f:
                linelist = line.strip().split(',')
                if linelist[0].startswith(('sample', 'snj', 'cpsnj')):
                    name = linelist[0].replace('.fa.gz','')
                    #print(name)
                    if name in sgb:
                        out.write(line)
                else:
                    out.write(line)


infile1 = r'D:\urban_soil\analysis\pre-calculated_data\sample_bins_info_mimag_gc.tsv'
infile2 = 'Cdb.csv'
outfile = 'Cdb_filter.csv'
cal_spire(infile1,infile2,outfile)

def spire_sgb(infile,outfile):
    with open(outfile,'wt') as out:
        with open(infile,'rt') as f:
            group = {}
            for line in f:
                linelist = line.strip().split(',')
                if linelist[1] not in group.keys():
                    group[linelist[1]] = [linelist[0]]
                else:
                    group[linelist[1]].append(linelist[0])
            n = 0
            y = 'no'
            for key,value in group.items():
                if len(value) > 1:
                    for item in value:
                        if item.startswith('spire') :
                            y = 'yes'  
                    if y == 'yes':
                        for item in value:
                            if item.startswith('sample') or item.startswith('snj') or item.startswith('cpsnj'):
                                n += 1 
                                overlap = ','.join(value)   
                                out.write(f'{item}\n')
                    if y == 'no':
                        print(f'{key}: {value}\n')
                    y = 'no'
    print(n)
infile = 'Cdb_filter.csv'
outfile = 'spire_overlap_sgb.tsv'
spire_sgb(infile,outfile)

def cal_soil(infile,outfile):
    with open(outfile,'wt') as out:
        with open(infile,'rt') as f:
            group = {}
            for line in f:
                linelist = line.strip().split(',')
                if linelist[1] not in group.keys():
                    group[linelist[1]] = [linelist[0]]
                else:
                    group[linelist[1]].append(linelist[0])
            n = 0
            y = 'no'
            for key,value in group.items():
                if len(value) > 1:
                    y = 'yes' if any(not s.startswith(('sample', 'snj', 'cpsnj')) for s in value) else 'no' 
                    if y == 'yes':
                        for item in value:
                            if item.startswith('sample') or item.startswith('snj') or item.startswith('cpsnj'):
                                n += 1    
                                overlap = ','.join(value)   
                                out.write(f'{key}\t{overlap}\n')
                    if y == 'no':
                        print(value)     
                    y = 'no'
        print(n)
infile = 'Cdb_soil.csv'
outfile = 'soil_overlap.tsv'
cal_soil(infile,outfile)

def cal_overlap(infile1,infile2):
    with open(infile1,'rt') as f:
        soil = set()
        for line in f:
            soil.add(line.strip())
    print(soil)
    with open(infile2,'rt') as f:
        spire = set()
        for line in f:
            spire.add(line.strip().replace('.gz',''))
    insert = soil.intersection(spire)
    print(f'Overlap: {len(insert)}')
    
infile1 = 'soil_overlap_sgb.tsv'
infile2 = 'spire_overlap_sgb.tsv'
cal_overlap(infile1,infile2)