# Use tRNAscan to predict tRNA for bacteria
prefix="sample1"
folder_path="./${prefix}/final_bins/"
result="./${prefix}/trnascan/result"

mkdir ${result}

for file in $(ls ${folder_path})
  do
    if [[ ${file} == *.fa ]]; then
      tRNAscan-SE ${folder_path}/${file} -B -o ${result}/${file}.trnascan.txt -f ${result}/${file}.trnascan.structure.txt -m ${result}/${file}.trnascan.stat.txt
    fi
  done