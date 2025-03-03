export PATH=/data/yiqian/software/ncbi-blast-2.14.0+/bin:$PATH

rpsblast -query all.mapped.smorfs_dedup.faa -out all.mapped.smorfs_dedup.tsv -db ~/istbi/GMSC/predict/negative/new/data/Cdd -num_threads 64 -evalue 0.01 -outfmt "6 qseqid sseqid qlen tlen score length pident evalue"