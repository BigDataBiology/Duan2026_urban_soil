prefix="sample1"
old_path="./08_gene_prediction/${prefix}/contig/${prefix}_contig.faa"
result="./08_gene_prediction/${prefix}/rgi/contig/"
mkdir ${result}
new_path="${result}/${prefix}_contig.faa"

sed 's/*//g' ${old_path}  > ${new_path}

rgi main -i ${new_path} -o ${result}/${prefix}_contig.rgi -n 8 -t protein --include_loose