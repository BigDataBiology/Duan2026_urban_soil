# Use CoverM to map reads to contigs and MAGs

#!/bin/bash

coverm contig --single /data/Projects/urban_soil/intermediate_data/00_quality_control/01_trim_filter_data/ont/porechop/snj19_ont_trim_filter_500_porechop.fastq -r /data/Projects/urban_soil/data/UrbanSoilAssemblies/snj19_medaka_polypolish.fasta.PolcaCorrected.fa.gz --output-file ./snj19.tsv --output-format dense --mapper minimap2-ont --methods mean -t 40
coverm genome --single /data/Projects/urban_soil/intermediate_data/00_quality_control/01_trim_filter_data/ont/porechop/23_ont_trim_filter_500_porechop.fastq --genome-fasta-directory /data/Projects/urban_soil/data/UrbanSoilMAGs/sample23 --genome-fasta-extension gz --output-file ./sample23.tsv --output-format dense --mapper minimap2-ont --methods mean relative_abundance trimmed_mean -t 64