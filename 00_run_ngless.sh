#!/bin/bash

ngsdata="/data/Public/soil_Release_X101SC23043779-Z01-J003_20230828/NGS/"
resultpath="/data/yiqian/soil/pipeline/00_quality_control/01_trim_filter_data"

for n in {15..15}
  do
    ngless -j 24 trim_filter_cpsnj${n}.ngl  ${ngsdata}/CPSNJ${n}/SNJ${n}_350.fq1.gz ${ngsdata}/CPSNJ${n}/SNJ${n}_350.fq2.gz ${resultpath}/cpsnj${n}_350.fq -t ./tmp_cpsnj${n}
  done