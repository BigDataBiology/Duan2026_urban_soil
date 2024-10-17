# Use barrnap to predict rRNA
prefix="sample1"
folder_path="./${prefix}/final_bins/"
result="./${prefix}/barrnap/"

mkdir ${result}

for file in $(ls ${folder_path})
  do
    barrnap --threads 4 --outseq ${result}/${file}.rRNA.faa ${folder_path}/${file} >${result}/${file}.barrnap.gff
  done