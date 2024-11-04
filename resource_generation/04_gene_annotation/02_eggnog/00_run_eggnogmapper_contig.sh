prefix="sample1"
path="./08_gene_prediction/${prefix}/contig/${prefix}_contig.faa"
result="./08_gene_prediction/${prefix}/eggnog/contig/"
mkdir ${result}

emapper.py -i ${path} -o ${result}/${prefix} --cpu 10

