def cal_spire_soil(infile):
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
                if y == 'no':
                    continue   
                y = 'no'
    print(n)
infile = '~/soil/pipeline/12_compare/result_spire/data_table/Cdb.csv'
cal_spire_soil(infile)

def cal_soil(infile):
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
                if y == 'no':
                    print(value)     
                y = 'no'
    print(n)
infile = '~/soil/pipeline/12_compare/result/data_table/Cdb.csv'
cal_soil(infile)