# Use tRNAscan to predict tRNA for archaea
prefix="cpsnj08"
folder_path="./${prefix}/trnascan/archaea"
result="./${prefix}/trnascan/result_arc"

mkdir ${result}

for file in $(ls ${folder_path})
  do
    if [[ ${file} == *.fa ]]; then
      tRNAscan-SE ${folder_path}/${file} -A -o ${result}/${file}.trnascan.txt -f ${result}/${file}.trnascan.structure.txt -m ${result}/${file}.trnascan.stat.txt
    fi
  done