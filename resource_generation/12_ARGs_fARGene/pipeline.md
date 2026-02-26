### Split mag sequences in several files
```bash
cd /work/microbiome/users/juan/Urban_soil/Genes/mags/data
python ../split_data.py /work/microbiome/urban_soil/data/UrbanSoilGenes/gene_prediction/mags/mags.faa.gz 10000000 &
```

### Run fARGene in mags data 
```bash
cd /work/microbiome/users/juan/Urban_soil/Genes/mags/
conda activate snakemake-env
snakemake -s Snakefile_mags -p -j 20 --use-conda
```

### Run fARGene in mags data 
```bash
cd /work/microbiome/users/juan/Urban_soil/Genes/mags/
conda activate snakemake-env
snakemake -s Snakefile_mags -p -j 20 --use-conda
```

### Run fARGene in mags data 
```bash
cd /work/microbiome/users/juan/Urban_soil/contigs/
snakemake -s Snakefile_mags -p -j 20 --use-conda
```

### Rename the sequences to avoid duplicate names
```bash
conda activate seqkit

cd /work/microbiome/users/juan/Urban_soil/contigs/fargene_predicted
for file in *.fasta; do
    seqkit rename -n "$file" -o file.tmp
    mv file.tmp "../fargene_predicted_suffix/$file"
done


cd /work/microbiome/users/juan/Urban_soil/Genes/mags/fargene_predicted
for file in *.fasta; do
    seqkit rename -n "$file" -o file.tmp
    mv file.tmp "../fargene_predicted_suffix/$file"
done
```

### Extract the name of the sequences, add the gene class and source file
```bash
cd /work/microbiome/users/juan/Urban_soil/contigs
for filename in fargene_predicted_suffix/*.fasta; do
    tag="${filename##*/}"
    #tag="${tag#*.gz-}"
    #tag="${tag%.fasta}"
    awk -v tag="$tag" '
        /^>/ { print ">" substr($0, 2) "@@@" tag; next }
        { print }
    ' "$filename" >> predicted_genes_with_class_clustering.fasta
done

cd /work/microbiome/users/juan/Urban_soil/Genes/mags
for filename in fargene_predicted_suffix/*.fasta; do
    tag="${filename##*/}"
    #tag="${tag#*.gz-}"
    #tag="${tag%.fasta}"
    awk -v tag="$tag" '
        /^>/ { print ">" substr($0, 2) "@@@" tag; next }
        { print }
    ' "$filename" >> predicted_genes_with_class_clustering.fasta
done
```


### Run CD-HIT on the predicted genes from MAGs
```bash
conda activate cdhit-env
cd /work/microbiome/users/juan/Urban_soil/Genes/mags
cd-hit -i ../predicted_genes_with_class_clustering.fasta -o cluster95 -c 0.95 -n 5 -M 0 -T 8
```

### BLAST the predicted genes to the centroids 
```bash
cd /work/microbiome/users/juan/Urban_soil/Genes/mags/diamond
diamond makedb --in ../cdhit/cluster95 -d ref_db
diamond blastp -q ../predicted_genes_with_class_clustering.fasta -d ref_db.dmnd -o closest_hits.tsv --max-target-seqs 1 --outfmt 6 qseqid sseqid pident evalue bitscore
```

### Run vsearch on th predicted genes from assembly
```bash
conda activate vsearch-env
cd /work/microbiome/users/juan/Urban_soil/contigs/
vsearch --cluster_fast ../predicted_genes_with_class_clustering.fasta --id 0.95 --centroids centroids95.fasta --uc clusters95.uc --threads
```

### Run RGI on the genes from MAGs
```bash
cd /work/microbiome/users/juan/Urban_soil/Genes/mags/rgi
zcat /work/microbiome/urban_soil/data/UrbanSoilGenes/gene_prediction/mags/mags.faa.gz | sed 's/\*//g' | gzip > mags_clean.faa.gz
conda activate RGI
rgi main -a DIAMOND -i mags_clean.faa.gz -o rgi  --local --clean -t protein -n 16
```


### Running AMRFinderPlus on fARGene results
```bash
cd /work/microbiome/users/juan/Urban_soil/Genes/contigs/clusters
conda activate seqkit
seqkit translate --frame 6 centroids.fasta |seqkit rename -n |seqkit sort -l -r | seqkit replace -p "_[\d]+$" -r "" | seqkit rmdup -n > centroids_longest_orf.faa

cd ../amrfinder
conda activate amrfinder
amrfinder -p ../clusters/centroids_longest_orf.faa > amrfinder_on_fargene_centroids95.tsv

conda activate fargene
prodigal -i ../clusters/centroids.fasta -a prodigal_centroids.faa  -p meta
amrfinder -p prodigal_centroids.faa > amrfinder_on_fargene_centroids95_prodigal.tsv


cd /work/microbiome/users/juan/Urban_soil/Genes/mags/amrfinder
amrfinder -p ../cdhit/cluster95 > amrfinder_on_fargene_centroids95.tsv

```


