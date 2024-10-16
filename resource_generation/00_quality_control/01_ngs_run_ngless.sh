#trim low quality and filter human reads
#!/bin/bash

ngsdata="/data/Public/soil_data/ngs"
resultpath="/data/yiqian/soil/pipeline/00_quality_control/01_trim_filter_data"

for n in {1..30}
  do
    ngless -j 24 trim_filter_sample${n}.ngl  ${ngsdata}/sample${n}/sample${n}_350.fq1.gz ${ngsdata}/sample${n}/sample${n}_350.fq2.gz ${resultpath}/sample${n}_350.fq -t ./tmp_sample${n}
  done