# Use NGLess to map short reads to contigs and MAGs

#!/bin/bash

datapath="/data/yiqian/soil/pipeline/00_quality_control/01_trim_filter_data/"
name="sample1"

ngless -t . -j 5 run_${name}.ngl ${datapath}/${name}_350.fq_trim_filter.pair.1.fq ${datapath}/${name}_350.fq_trim_filter.pair.2.fq ${datapath}/${name}_350.fq_trim_filter.singles.fq ./08_mapping/mag/final_bins_fa/${name}_final_bins.fa.gz ${name}_mag.sam