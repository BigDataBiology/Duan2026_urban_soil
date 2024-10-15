# Use Polypolish to polish with short reads
# All the contigs are too large to run, so split contigs into to files

for n in {1..2}
  do
    bwa index medaka_${n}.fasta
    bwa mem -t 24 -a medaka_${n}.fasta sample1_350.fq_trim_filter.pair.1.fq >medaka_${n}_mapped_1.sam
    bwa mem -t 24 -a medaka_${n}.fasta sample1_350.fq_trim_filter.pair.2.fq >medaka_${n}_mapped_2.sam

    polypolish_insert_filter.py --in1 medaka_${n}_mapped_1.sam --in2 medaka_${n}_mapped_2.sam --out1 medaka_${n}_mapped_filtered_1.sam --out2 medaka_${n}_mapped_filtered_2.sam
    polypolish medaka_${n}.fasta medaka_${n}_mapped_filtered_1.sam medaka_${n}_mapped_filtered_2.sam >medaka_${n}_polypolish.fasta
    
    rm *.sam
  done