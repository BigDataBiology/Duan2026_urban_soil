def bgc_length(infile1):
    length_dict = {}
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            if line.startswith('\t'):
                continue
            else:
                contig = linelist[5]
                number = linelist[7].replace('\'','').split('.')[1]
                start = int(linelist[9])
                end = int(linelist[10])
                genome = linelist[12]
                name = f'{genome}.fa_{contig}.region{number.zfill(3)}'
                length = end - start + 1
                length_dict[name] = f'{start}\t{end}\t{length}'
    return length_dict  

def map_length(length_dict,infile2,outfile):
    with open(outfile,'wt') as out:
        with open(infile2,'rt') as f:
            for line in f:
                if line.startswith('Record'):
                    continue
                else:
                    linelist = line.strip().split('\t')
                    out.write(f'{linelist[1]}\t{linelist[4]}\t{linelist[5]}\t{length_dict[linelist[1]]}\n')

infile1 = r'D:\soil\bgc\bigscape\bgc_region.txt'
infile2 = r'D:\soil\bgc\bigscape\record_annotations.tsv'
outfile = r'D:\soil\bgc\bigscape\bgc_region_length.tsv'
length_dict = bgc_length(infile1)
map_length(length_dict,infile2,outfile)

def map_gtdb(infile1,infile2,outfile):
    gtdb_dict = {}
    with open(infile1, 'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            gtdb_dict[linelist[1]] = linelist[16]
    with open(outfile, 'wt') as out:
        with open(infile2, 'rt') as f:
            for line in f:
                linelist = line.strip().split('\t')
                genome = linelist[0].split('.')[0]
                if len(gtdb_dict[genome].split(';')) == 7:
                    d = gtdb_dict[genome].split(';')[0]
                    p = gtdb_dict[genome].split(';')[1]
                    c = gtdb_dict[genome].split(';')[2]
                    o = gtdb_dict[genome].split(';')[3]
                    fa = gtdb_dict[genome].split(';')[4]
                    g = gtdb_dict[genome].split(';')[5]
                    s = gtdb_dict[genome].split(';')[6]
                    out.write(f'{line.strip()}\t{genome}\t{gtdb_dict[genome]}\t{d}\t{p}\t{c}\t{o}\t{fa}\t{g}\t{s}\n')
                else:
                    out.write(f'{line.strip()}\t{genome}\t{gtdb_dict[genome]}\n')

infile1 = r'D:\soil\bgc\bigscape\mag_gtdb_r226.txt'
infile2 = r'D:\soil\bgc\bigscape\bgc_region_length.tsv'
outfile = r'D:\soil\bgc\bigscape\bgc_region_len_gtdb.tsv'
map_gtdb(infile1, infile2, outfile)

def map_arg(infile1,infile2,outfile):
    arg = set()
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split('.gbk:')
            arg.add(linelist[0])
    with open(outfile,'wt') as out:
        with open(infile2,'rt') as f:
            for line in f:
                linelist = line.strip().split('\t')
                if linelist[0] in arg:
                    out.write(f'{line.strip()}\tARG\n')

infile1 = r'D:\soil\bgc\bigscape\arg\bgc_arg.tsv'
infile2 = r'D:\soil\bgc\bigscape\bgc_region_len_gtdb.tsv'
outfile = r'D:\soil\bgc\bigscape\bgc_region_len_gtdb_arg.tsv'
map_arg(infile1,infile2,outfile)