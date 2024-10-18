prefix="sample1"
folder_path="./08_gene_prediction/${prefix}/genes/"
new_path="./08_gene_prediction/${prefix}/rgi/result"

mkdir rgi
mkdir ${new_path}

cd ${folder_path}
for file in *.faa
  do
    sed 's/*//g' ${folder_path}/${file}  > ${new_path}/${file}
  done

for file in $(ls ${new_path})
  do
    rgi main -i ${new_path}/${file} -o ${new_path}/${file}.rgi -t protein -n 8 --include_loose
  done