# Use NGless to map reads to polished contigs
# Use SemiBin2 to generate bins

ngless map.ngl -j 24 -t tmp

samtools view -h -b -S medaka_polypolish.fasta.PolcaCorrected.sam -o medaka_polypolish.fasta.PolcaCorrected.bam
samtools view -b -F 4 medaka_polypolish.fasta.PolcaCorrected.bam -o medaka_polypolish.fasta.PolcaCorrected_mapped.bam
samtools sort medaka_polypolish.fasta.PolcaCorrected_mapped.bam -o medaka_polypolish.fasta.PolcaCorrected_mapped_sorted.bam

samtools index medaka_polypolish.fasta.PolcaCorrected_mapped_sorted.bam

SemiBin2 single_easy_bin --environment soil -i medaka_polypolish.fasta.PolcaCorrected.fa -b medaka_polypolish.fasta.PolcaCorrected_mapped_sorted.bam -o semibin_output --sequencing-type long_read -p 24

rm *.sam
rm *.bam