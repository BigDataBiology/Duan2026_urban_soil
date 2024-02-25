#!/bin/bash

export PATH=/data/yiqian/software/target/release:$PATH

gunzip -c /data/Public/soil_Release_X101SC23043779-Z01-J003_20230828/ONT/SNJ13/fastq_pass.gz >snj13_fastq_pass.fastq

less snj13_fastq_pass.fastq | chopper -q 10 -l 500 --threads 20 --contam snj13_ont_trim_filter_500.fastq

porechop_abi -abi -i snj13_ont_trim_filter_500.fastq --format fastq --verbosity 1 --discard_middle -t 40 -o snj13_ont_trim_filter_500_porechop.fastq > snj13_ont_trim_filter_500_porechop.txt