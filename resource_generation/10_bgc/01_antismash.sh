prefix="snj01"
folder_path="~/soil/08_gene_prediction/${prefix}/final_bins"
result="~/soil/15_bgc/${prefix}/"

mkdir ${result}

for file in $(ls ${folder_path})
  do
    cd ${result}
    every="${result}/${file}"
    antismash ${folder_path}/${file} --genefinding-tool prodigal --output-dir ${every} -c 8
  done
