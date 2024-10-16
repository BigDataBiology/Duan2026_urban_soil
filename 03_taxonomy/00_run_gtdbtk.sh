#Use gtdb-tk to assign taxonomy to MAGs
gtdbtk classify_wf --cpus 20 --genome_dir ./final_bins --out_dir gtdb_result --extension .fa.gz --mash_db .