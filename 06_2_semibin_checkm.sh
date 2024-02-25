export PATH=~/miniconda3/bin:$PATH

source activate /data/yiqian/miniconda3/envs/soilnew
ngless map.ngl -j 24 -t tmp

source activate /data/yiqian/mambaforge/envs/semibin

samtools view -h -b -S medaka_polypolish.fasta.PolcaCorrected.sam -o medaka_polypolish.fasta.PolcaCorrected.bam
samtools view -b -F 4 medaka_polypolish.fasta.PolcaCorrected.bam -o medaka_polypolish.fasta.PolcaCorrected_mapped.bam
samtools sort medaka_polypolish.fasta.PolcaCorrected_mapped.bam -o medaka_polypolish.fasta.PolcaCorrected_mapped_sorted.bam

samtools index medaka_polypolish.fasta.PolcaCorrected_mapped_sorted.bam

SemiBin2 single_easy_bin --environment soil -i medaka_polypolish.fasta.PolcaCorrected.fa -b medaka_polypolish.fasta.PolcaCorrected_mapped_sorted.bam -o semibin_output --sequencing-type long_read -p 24

source activate /data/yiqian/mambaforge/envs/checkm2

checkm2 predict -i ./semibin_output/output_bins -o checkm_results --database_path /data/yiqian/software/checkm_database/uniref100.KO.1.dmnd -t 24 -x .fa.gz