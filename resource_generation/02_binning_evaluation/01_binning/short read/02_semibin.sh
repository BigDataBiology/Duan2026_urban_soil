export PATH=~/miniconda3/bin:$PATH

source activate ~/miniconda3/envs/ngless
ngless map.ngl -j 40

export PATH=~/software/samtools-1.18:$PATH

samtools view -h -b -S final_contigs.sam -o final_contigs.bam
samtools view -b -F 4 final_contigs.bam -o final_contigs_mapped.bam
samtools sort final_contigs_mapped.bam -o final_contigs_mapped_sorted.bam

samtools index final_contigs_mapped_sorted.bam

source activate semibin
SemiBin2 single_easy_bin --environment soil -i ~/soil/02_assembly/megahit/cpsnj01/final.contigs.fa -b final_contigs_mapped_sorted.bam -o semibin_output