prefix="sample1"
folder_path="./${prefix}/final_bins/"
new_path="./08_gene_prediction/${prefix}/final_bins/"
result="./08_gene_prediction/${prefix}/genes/"

mkdir ${new_path}
mkdir ${result}

for file in $(ls ${folder_path})
  do
    cp ${folder_path}/${file} ${new_path}/${prefix}_${file}
  done

for file in $(ls ${new_path})
  do
    prodigal -i ${new_path}/${file} -o ${result}/${file}.gff -a ${result}/${file}.faa -d ${result}/${file}.fna
  done