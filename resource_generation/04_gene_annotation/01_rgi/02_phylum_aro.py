def cal_aro(infile):
    sh = set()
    nj = set()
    with open(infile,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            if linelist[1].startswith('cpsnj') or linelist[1].startswith('snj'):
                arolist = linelist[19].split(';')
                for item in arolist:
                    if item != 'NA' and item not in nj:
                        nj.add(item)
            if linelist[1].startswith('sample'):
                arolist = linelist[19].split(';')
                for item in arolist:
                    if item != 'NA' and item not in sh:
                        sh.add(item)
    print(len(sh))
    print(len(nj))
    print(len(sh.intersection(nj)))

def cal_tax(infile):
    phylum = {}
    with open(infile,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            if linelist[9] == 'T' and linelist[19] != 'NA':
                tax = linelist[16]
                taxlist = tax.split(';')
                if len(taxlist) > 1:
                    if taxlist[1] not in phylum.keys():
                        phylum[taxlist[1]] = 1
                    else:
                        phylum[taxlist[1]] += 1
    print(phylum)

def cal_tax_aro(infile,outfile):
    phylum = {}
    with open(infile,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            if linelist[9] == 'T' and linelist[19] != 'NA':
                tax = linelist[16]
                taxlist = tax.split(';')
                if len(taxlist) > 1:
                    if taxlist[1] not in phylum.keys():
                        phylum[taxlist[1]] = set()
                    arolist = linelist[19].split(';')
                    for item in arolist:
                        if item not in phylum[taxlist[1]]:
                            phylum[taxlist[1]].add(item)

    with open(outfile,'wt') as out:
        for key,value in phylum.items():
            out.write(f'{key}\t{len(value)}\n')

infile = r'sample_bins_info_mimag.tsv'
outfile = 'phylum_aro.tsv'
#cal_aro(infile)
#cal_tax(infile)
cal_tax_aro(infile,outfile)