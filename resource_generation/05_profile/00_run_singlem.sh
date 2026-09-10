# Use SingleM to profile
export SINGLEM_METAPACKAGE_PATH='S3.2.1.GTDB_r214.metapackage_20231006.smpkg.zb'

singlem pipe -1 ~/00_quality_control/00_raw_clean_data/ngs/CPSNJ01/SNJ01_350.fq1.gz -2 ~/00_quality_control/00_raw_clean_data/ngs/CPSNJ01/SNJ01_350.fq2.gz -p output.profile.tsv --threads 4 --otu-table otutable.tsv --taxonomic-profile-krona krona.html

#singlem microbial_fraction -p output.profile.tsv  --output-tsv read_fraction.tsv --input-metagenome-sizes size.tsv

singlem summarise --input-taxonomic-profiles output.profile.tsv --output-species-by-site-relative-abundance genus-by-site-relative-abundance.tsv --output-species-by-site-level genus