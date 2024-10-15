#Use fastqc for quality checking raw short reads
#!/bin/bash

ngsdata="/data/Public/soil_data/ngs"
resultpath="/data/yiqian/soil/pipeline/00_quality_control/00_raw_clean_data/ngs"

for n in {1..30}
  do
    gunzip -c ${ngsdata}/sample${n}/${n}_350.fq1.gz > ${resultpath}/${n}_350_1.fq
    gunzip -c ${ngsdata}/sample${n}/${n}_350.fq2.gz > ${resultpath}/${n}_350_2.fq
    fastqc -f fastq -o ${resultpath} ${resultpath}/${n}_350_1.fq ${resultpath}/${n}_350_2.fq -t 64
  done