prefix="sample1"
folder_path="./08_gene_prediction/${prefix}/genes/"
result="./08_gene_prediction/${prefix}/eggnog/result"

cd ${folder_path}
for file in *.faa
  do
    emapper.py -i ${folder_path}/${file} -o ${result}/${file} --cpu 10
  done