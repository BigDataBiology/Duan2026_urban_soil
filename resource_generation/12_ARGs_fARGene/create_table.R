library(dplyr)
library(Biostrings)



setwd("~/Documents/GitHub/urban_soil/analysis/")

gene_classes <- data.frame(rbind(
  c("aac2p",	"aminoglycoside acetyltransferase",	"aac(2')", "aminoglycoside"),
  c("aac3_1",	"aminoglycoside acetyltransferase",	"aac(3)", "aminoglycoside"), c("aac3_2",	"aminoglycoside acetyltransferase",	"aac(3)", "aminoglycoside"),
  c("aac6p_1",	"aminoglycoside acetyltransferase",	"aac(6')", "aminoglycoside"), c("aac6p_2",	"aminoglycoside acetyltransferase",	"aac(6')", "aminoglycoside"),
  c("aac6p_3",	"aminoglycoside acetyltransferase",	"aac(6')", "aminoglycoside"), c("aph2b",	"aminoglycoside phosphotransferase", "aph(2'')", "aminoglycoside"),
  c("aph3p",	"aminoglycoside phosphotransferase", "aph(3')", "aminoglycoside"), c("aph6p",	"aminoglycoside phosphotransferase", "aph(6)", "aminoglycoside"),
  c("class_a",	"beta-lactamase A",	"beta-lactamase A", "beta-lactamase"), c("class_b1_b2",	"beta-lactamase B",	"beta-lactamase B1-B2", "beta-lactamase"),
  c("class_b3",	"beta-lactamase B",	"beta-lactamase B3", "beta-lactamase"), c("class_c",	"beta-lactamase C",	"beta-lactamase C", "beta-lactamase"),
  c("class_d1",	"beta-lactamase D",	"beta-lactamase D", "beta-lactamase"), c("class_d2",	"beta-lactamase D",	"beta-lactamase D", "beta-lactamase"),
  c("erm_1",	"macrolide erm", "erm", "macrolide"), c("erm_2",	"macrolide erm", "erm",	'macrolide'),
  c("mph",	"macrolide mph", "mph",	"macrolide"), c("qnr",	"qnr", "qnr", "quinolone"), c("tet_efflux",	"tetracycline efflux pump", "tet efflux", "tetracycline"),
  c("tet_enzyme",	"tetracycline inactivating enzymes", "tet enzyme", "tetracycline"), c("tet_rpg",	"tetracycline ribosomal protection genes", "tet RPG", "tetracycline")))

metadata <- readxl::read_excel("../resource_generation/12_ARGs_fARGene/supp1_sample_metadata.xlsx")
metadata <- metadata %>% mutate(Sample_id = gsub("ample", "", Sample_id)) 

centroids <- read.table("../resource_generation/12_ARGs_fARGene/contigs/orfs-amino-centroids.txt", quote="\"", comment.char="")
genes <- read.table("../resource_generation/12_ARGs_fARGene/contigs/closest_hits_orfs_amino.tsv", quote="\"", comment.char="")
genes <- genes %>% rename(orf = V1, centroid = V2) %>% select(orf, centroid)
genes <- genes %>% mutate(centroid = ifelse(orf %in% centroids$V1, orf, centroid))

genes$sample <- sapply(strsplit(genes$orf, split = "@@@"), function(x) x[1])
genes$hmm <- sapply(strsplit(genes$orf, split = "@@@"), function(x) x[2])
genes$c_hmm <- sapply(strsplit(genes$centroid, split = "@@@"), function(x) x[2])
genes$sequence <- sapply(strsplit(genes$orf, split = "@@@"), function(x) x[3])
genes$contig <- sapply(strsplit(genes$sequence, split = "_"), function(x) paste(x[1], x[2], x[3], sep = "_"))
genes <- genes %>% mutate(hmm = c_hmm)
genes$c_hmm <- NULL


######################################
######################################

## add metadata 


genes <- genes %>% mutate(description = gene_classes$X2[match(hmm, gene_classes$X1)],
                          class = gene_classes$X3[match(hmm, gene_classes$X1)],
                          hclass = gene_classes$X4[match(hmm, gene_classes$X1)],
                          City = metadata$City[match(sample, metadata$Sample_id)],
                          Location = metadata$Location[match(sample, metadata$Sample_id)],
                          Date = metadata$Date[match(sample, metadata$Sample_id)],
                          Longitude = metadata$Longitude[match(sample, metadata$Sample_id)],
                          Latitude = metadata$Latitude[match(sample, metadata$Sample_id)])


genes <- genes %>% 
  mutate(
    hclass = 
      factor(hclass, levels = c("aminoglycoside", "beta-lactamase", 
                                "macrolide", "tetracycline", "quinolone")))

## remove sample s4 and s11

genes <- genes %>% 
  filter(!sample %in% c("s4", "s11"))

colnames(genes)

genes <- genes %>% select(orf, centroid, sample, contig, class, hclass)


fasta <- read.table("../resource_generation/12_ARGs_fARGene/contigs/predicted-orfs-amino-nonredundant.fasta", quote="\"", comment.char="")
idx <- which(grepl("^>", fasta$V1))
 
fasta_df <- data.frame(
  orf = sub("^>", "", fasta$V1[idx]),
  sequence = sapply(seq_along(idx), function(i) {
    start <- idx[i] + 1
    end <- if (i < length(idx)) idx[i + 1] - 1 else nrow(fasta)
    paste(fasta$V1[start:end], collapse = "")
  }),
  stringsAsFactors = FALSE
)


genes <- genes %>% left_join(fasta_df, by = "orf")

resfinder <- read.delim("~/Documents/GitHub/urban_soil/resource_generation/12_ARGs_fARGene/contigs/ResFinder_results_tab.txt")
genes <- genes %>% mutate(resfinder = resfinder$Resistance.gene[match(orf , resfinder$Contig)])

write.table(
  genes,
  file = "../resource_generation/12_ARGs_fARGene/ARGs.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)





