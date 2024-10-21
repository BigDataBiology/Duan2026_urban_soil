# Use dRep to deduplicate MQ MAGs
dRep dereplicate ./09_drep/result --genomeInfo ./09_drep/checkm/checkm_quality_report.tsv -g ./09_drep/data/*.fa -p 40 -comp 50 -con 10 -sa 0.95
dRep dereplicate ./09_drep/result_99 --genomeInfo ./09_drep/checkm/checkm_quality_report.tsv -g ./09_drep/data/*.fa -p 40 -comp 50 -con 10 -sa 0.99