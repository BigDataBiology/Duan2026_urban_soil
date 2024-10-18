prefix="sample1"
new_path="./08_gene_prediction/${prefix}/medaka_polypolish.fasta.PolcaCorrected.fa"
result="./08_gene_prediction/${prefix}/contig/"

mkdir ${result}

prodigal -i ${new_path} -o ${result}/${prefix}_contig.gff -a ${result}/${prefix}_contig.faa -d ${result}/${prefix}_contig.fna -p meta