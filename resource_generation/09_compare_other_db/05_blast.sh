makeblastdb -in ~/soil/13_compare/marker/marker_result/soil_cog0018.fa -dbtype nucl -out soil_cog0018_db
blastn -db soil_cog0018_db -query ~/soil/13_compare/marker/marker_result/sgb_cog0018.fa -out cog_0018.tsv -perc_identity 95 -outfmt "6 qseqid qlen sseqid slen qstart qend sstart send length qcovs pident evalue" -num_threads 20
awk '$10 > 90 {print $1}' cog_0018.tsv | sort | uniq | wc -l