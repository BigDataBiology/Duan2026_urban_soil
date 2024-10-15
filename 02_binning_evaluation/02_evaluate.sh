# Use checkm2 and GUNC to evaluate MAG quality

checkm2 predict -i ./semibin_output/output_bins -o checkm_results --database_path uniref100.KO.1.dmnd -t 24 -x .fa.gz

conda activate gunc
mkdir gunc_result
gunc run -d ./semibin_output/output_bins -o gunc_result -t 20 -r gunc_db_progenomes2.1.dmnd -e .fa.gz

awk '$2 > 50 && $3 < 10{print $1}' ./checkm_results/quality_report.tsv|sed 's/.fa//' > pass_checkm.txt
awk '$13 == "True"{print $1}' ./gunc_result/GUNC.progenomes_2.1.maxCSS_level.tsv >pass_gunc.txt
sort pass_checkm.txt pass_gunc.txt|uniq -d >pass_checkm_gunc.txt
mkdir final_bins
for line in `cat pass_checkm_gunc.txt`
  do
    cp ./semibin_output/output_bins/${line}.fa.gz ./final_bins
  done