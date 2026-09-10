singlem summarise --input-otu-tables ~/06_profile/*/otutable.tsv --output-otu-table all_otu_table.csv
singlem summarise --input-otu-table all_otu_table.csv --unifrac-by-otu all_otu_table
for f in *.unifrac
  do
    echo $f >> total_counts_by_otu_marker.txt
    awk '{gsub(/[^0-9]/, "", $3); sum+=$3} END {print sum}' $f >> total_counts_by_otu_marker.txt
  done
convertToEBD.py all_otu_table.S3.5.ribosomal_protein_S2_rpsB.unifrac all_otu_table.S3.5.rib_prot_S2_rpsB.ebd
convertToEBD.py all_otu_table.S3.18.EIF_2_alpha.unifrac all_otu_table.S3.18.EIF_2_alpha.ebd
convertToEBD.py all_otu_table.S3.42.RNA_pol_A_bac.unifrac all_otu_table.S3.42.RNA_pol_A_bac.ebd
convertToEBD.py all_otu_table.S3.37.uS4_arch.unifrac all_otu_table.S3.37.uS4_arch.ebd