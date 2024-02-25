export PATH=/data/yiqian/software/bwa-0.7.17:$PATH

bwa index medaka_2.fasta
bwa mem -t 24 -a medaka_2.fasta /data/yiqian/soil/pipeline/00_quality_control/01_trim_filter_data/cpsnj16_350.fq_trim_filter.pair.1.fq >medaka_2_mapped_1.sam
bwa mem -t 24 -a medaka_2.fasta /data/yiqian/soil/pipeline/00_quality_control/01_trim_filter_data/cpsnj16_350.fq_trim_filter.pair.2.fq >medaka_2_mapped_2.sam

export PATH=/data/yiqian/software:$PATH

polypolish_insert_filter.py --in1 medaka_2_mapped_1.sam --in2 medaka_2_mapped_2.sam --out1 medaka_2_mapped_filtered_1.sam --out2 medaka_2_mapped_filtered_2.sam
polypolish medaka_2.fasta medaka_2_mapped_filtered_1.sam medaka_2_mapped_filtered_2.sam >medaka_2_polypolish.fasta