sample = 'sample9'

with open(f'./barrnap/{sample}.tsv','wt') as out:
    with open(f'./barrnap/{sample}.err','rt',encoding='gbk',errors='ignore') as f:
        s16 = 0
        s23 = 0
        s5 = 0
        for line in f:
            if line.startswith('[barrnap] Found: 16S'):
                s16 +=1
            if line.startswith('[barrnap] Found: 23S'):
                s23 +=1
            if line.startswith('[barrnap] Found: 5S'):
                s5 +=1
            if line.startswith('[barrnap] Writing'):
                #bin = line.split(': ')[1].split('.')[0]
                bin = line.split('//')[1].split('.')[0]
                name = f'{sample}_{bin}'
                out.write(f'{name}\t{s16}\t{s23}\t{s5}\n')
                s16 = 0
                s23 = 0
                s5 = 0