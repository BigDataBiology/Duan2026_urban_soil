def cal():
    sample = 'sample9'
    aalist = ['Ala','Gly','Pro','Thr','Val','Ser','Leu','Arg','Phe','Asn','Lys','Asp','Glu','His','Gln','Ile','Met/fMet','Tyr','Cys','Trp']
    with open(f'./trnascan/{sample}.tsv','wt') as out:
        with open(f'./trnascan/{sample}.txt','rt') as f:
            aa = 0
            for line in f:
                if line.startswith('Search statistics saved in:'):
                    bin = line.strip().split('/')[-1].split('.')[0]
                    name = f'{sample}_{bin}'
                #if line.startswith('Total tRNAs:'):
                    #trna = line.strip().replace('Total tRNAs:                                ','')
                for item in aalist:
                    if line.startswith(item):
                        trna = line.split(': ')[1][0]
                        if trna != '0':
                            aa += 1
                if line.startswith('SelCys'):
                    trna = line.split(': ')[1][0]
                    if trna != '0':
                        aa += 1
                    out.write(f'{name}\t{aa}\n')
                    aa = 0
cal()