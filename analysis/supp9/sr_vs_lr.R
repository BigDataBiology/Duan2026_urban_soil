library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)
library(ggpattern)
library(grid)
library(cowplot)
library(scales)
library(ComplexHeatmap)
library(geosphere)
library(ggbreak)
library(ggforce)
library(ggpmisc)
library(Cairo)
library(e1071)
library(car)
library(emmeans)
library(MASS)
library(glmmTMB)
library(ggrepel)
library(ggpp)
library(broom.mixed)

pal_8 <- brewer.pal(8, "Dark2")
pal_8 <- pal_8[c(1,6,2:5,7,8)]


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

# when i re-run in short-reads i corrected aph(6') to aph(6)
gene_classes_2 <- data.frame(rbind(
  c("aac2p",	"aminoglycoside acetyltransferase",	"aac(2')", "aminoglycoside"),
  c("aac3_1",	"aminoglycoside acetyltransferase",	"aac(3)", "aminoglycoside"), c("aac3_2",	"aminoglycoside acetyltransferase",	"aac(3)", "aminoglycoside"),
  c("aac6p_1",	"aminoglycoside acetyltransferase",	"aac(6')", "aminoglycoside"), c("aac6p_2",	"aminoglycoside acetyltransferase",	"aac(6')", "aminoglycoside"),
  c("aac6p_3",	"aminoglycoside acetyltransferase",	"aac(6')", "aminoglycoside"), c("aph2b",	"aminoglycoside phosphotransferase", "aph(2'')", "aminoglycoside"),
  c("aph3p",	"aminoglycoside phosphotransferase", "aph(3')", "aminoglycoside"), c("aph6",	"aminoglycoside phosphotransferase", "aph(6)", "aminoglycoside"),
  c("class_a",	"beta-lactamase A",	"beta-lactamase A", "beta-lactamase"), c("class_b1_b2",	"beta-lactamase B",	"beta-lactamase B1-B2", "beta-lactamase"),
  c("class_b3",	"beta-lactamase B",	"beta-lactamase B3", "beta-lactamase"), c("class_c",	"beta-lactamase C",	"beta-lactamase C", "beta-lactamase"),
  c("class_d1",	"beta-lactamase D",	"beta-lactamase D", "beta-lactamase"), c("class_d2",	"beta-lactamase D",	"beta-lactamase D", "beta-lactamase"),
  c("erm_1",	"macrolide erm", "erm", "macrolide"), c("erm_2",	"macrolide erm", "erm",	'macrolide'),
  c("mph",	"macrolide mph", "mph",	"macrolide"), c("qnr",	"qnr", "qnr", "quinolone"), c("tet_efflux",	"tetracycline efflux pump", "tet efflux", "tetracycline"),
  c("tet_enzyme",	"tetracycline inactivating enzymes", "tet enzyme", "tetracycline"), c("tet_rpg",	"tetracycline ribosomal protection genes", "tet RPG", "tetracycline")))


metadata <- readxl::read_excel("../resource_generation/12_ARGs_fARGene/supp1_sample_metadata.xlsx")
metadata <- metadata %>% mutate(Sample_id = gsub("ample", "", Sample_id)) 


## centroids_sr each of the genes found in sr clustered
centroids_sr <- read.table("../resource_generation/12_ARGs_fARGene/sr_contigs/cluster_membership_sr.tsv", quote="\"", comment.char="")


## blast of the long-read ARGs to short-read ARGs
blasted_lr_vs_sr <- read.table("../resource_generation/12_ARGs_fARGene/long_orfs_nonredundant_vs_short_merged_orfs.tsv", quote="\"", comment.char="")
blasted_lr_vs_sr$sample_lr <- sapply(strsplit(blasted_lr_vs_sr$V1, split = "@@@"), function(x) x[1])
blasted_lr_vs_sr$sample_sr <- sapply(strsplit(blasted_lr_vs_sr$V2, split = "@@@"), function(x) x[2])
blasted_lr_vs_sr$sample_sr <- gsub("sample", "s", blasted_lr_vs_sr$sample_sr)
blasted_lr_vs_sr0 <- blasted_lr_vs_sr

# since we are looking for genes detected in short-read and long-read, the sample lr and sr in blast should match
blasted_lr_vs_sr <- blasted_lr_vs_sr[blasted_lr_vs_sr$sample_lr == blasted_lr_vs_sr$sample_sr,]

## all short-read ARGs

sr_genes <- centroids_sr %>% rename(orf = V1, centroid = V2) %>% dplyr::select(orf, centroid)
sr_genes$sequence <- sapply(strsplit(sr_genes$orf, split = "@@@"), function(x) x[1])
sr_genes$sample <- sapply(strsplit(sr_genes$orf, split = "@@@"), function(x) x[2])
sr_genes$sample <- gsub("sample", "s", sr_genes$sample)
sr_genes$hmm <- sapply(strsplit(sr_genes$orf, split = "@@@"), function(x) x[3])
sr_genes$c_hmm <- sapply(strsplit(sr_genes$centroid, split = "@@@"), function(x) x[3])
sr_genes$contig <- sapply(strsplit(sr_genes$sequence, split = "_"), function(x) paste(x[1], x[2], x[3], sep = "_"))
sr_genes <- sr_genes %>% mutate(hmm = c_hmm)
sr_genes$c_hmm <- NULL


sr_genes <- sr_genes %>% mutate(description = gene_classes_2$X2[match(hmm, gene_classes_2$X1)],
                          class = gene_classes_2$X3[match(hmm, gene_classes_2$X1)],
                          hclass = gene_classes_2$X4[match(hmm, gene_classes_2$X1)],
                          City = metadata$City[match(sample, metadata$Sample_id)],
                          Location = metadata$Location[match(sample, metadata$Sample_id)],
                          Date = metadata$Date[match(sample, metadata$Sample_id)],
                          Longitude = metadata$Longitude[match(sample, metadata$Sample_id)],
                          Latitude = metadata$Latitude[match(sample, metadata$Sample_id)])


sr_genes <- sr_genes %>% 
  mutate( hclass =  factor(hclass, levels = c("aminoglycoside", "beta-lactamase", 
              "macrolide", "tetracycline", "quinolone")))


# contig statistics for long-read assemblies
contig_stat_lr <- read_table("../resource_generation/12_ARGs_fARGene/stat_lr.txt", 
                             locale = locale(decimal_mark = ".", grouping_mark = ","))
contig_stat_lr$sample <- sapply(strsplit(sapply(strsplit(contig_stat_lr$file, split = "/"), function(x) x[7]), split = "_"), function(x) x[1])
contig_stat_lr$sample <- gsub("sample", "s", contig_stat_lr$sample)

# contig statistics for short-read assemblies
contig_stat <- read.table("../resource_generation/12_ARGs_fARGene/sr_contigs/stat_sr.txt", quote="\"", comment.char="", header = T)
contig_stat$sample <- sapply(strsplit(sapply(strsplit(contig_stat$file, split = "/"), function(x) x[7]), split = "_"), function(x) x[1])
contig_stat$sample <- gsub("sample", "s", contig_stat$sample)


# length of the ARGs from long-read ARGs
lines <- readLines("../resource_generation/12_ARGs_fARGene/contigs/predicted-orfs-amino-nonredundant.fasta")

header_idx <- grep("^>", lines)
headers <- sub("^>", "", lines[header_idx])

seq_lengths <- sapply(seq_along(header_idx), function(i) {
  start <- header_idx[i] + 1
  end <- if (i < length(header_idx)) header_idx[i + 1] - 1 else length(lines)
  seq <- paste(lines[start:end], collapse = "")
  seq <- sub("\\*$", "", seq)  # drop trailing stop codon marker, not a real residue
  nchar(seq)
})

# all genes with lengths for lr ARGs
gene_length_lr <- data.frame(header = headers, length = seq_lengths, stringsAsFactors = FALSE)

gene_length_lr$class <- sapply(strsplit(gene_length_lr$header, split = "@@@"), function(x) x[2])
gene_length_lr$class <- gene_classes$X3[match(gene_length_lr$class, gene_classes$X1)]
gene_length_lr$sample <- sapply(strsplit(gene_length_lr$header, split = "@@@"), function(x) x[1])
gene_length_lr$sample <- gsub("sample", "s", gene_length_lr$sample)
gene_length_lr <- gene_length_lr %>% mutate(Location = metadata$Location[match(sample, metadata$Sample_id)])
gene_length_lr <- gene_length_lr %>% mutate(City = metadata$City[match(sample, metadata$Sample_id)])

# length of the ARGs from short-read ARGs
lines <- readLines("../resource_generation/12_ARGs_fARGene/sr_contigs/merged-orf-amino.fasta")
header_idx <- grep("^>", lines)
headers <- sub("^>", "", lines[header_idx])

seq_lengths <- sapply(seq_along(header_idx), function(i) {
  start <- header_idx[i] + 1
  end <- if (i < length(header_idx)) header_idx[i + 1] - 1 else length(lines)
  seq <- paste(lines[start:end], collapse = "")
  seq <- sub("\\*$", "", seq)  # drop trailing stop codon marker, not a real residue
  nchar(seq)
})

gene_length_sr <- data.frame(header = headers, length = seq_lengths, stringsAsFactors = FALSE)

gene_length_sr$class <- sapply(strsplit(gene_length_sr$header, split = "@@@"), function(x) x[3])
gene_length_sr$class <- gene_classes_2$X3[match(gene_length_sr$class, gene_classes_2$X1)]
gene_length_sr$sample <- sapply(strsplit(gene_length_sr$header, split = "@@@"), function(x) x[2])
gene_length_sr$sample <- gsub("sample", "s", gene_length_sr$sample)
gene_length_sr$City <- metadata$City[match(gene_length_sr$sample, metadata$Sample_id)]
gene_length_sr$Location <- metadata$Location[match(gene_length_sr$sample, metadata$Sample_id)]


# gene_length_lr_2 add the dummy, was it detected in short reads?
gene_length_lr_2 <- gene_length_lr %>% 
  mutate(sr80 = ifelse(header %in% blasted_lr_vs_sr$V1[blasted_lr_vs_sr$V3 > 80 & blasted_lr_vs_sr$V13 > 80], 1, 0),
         sr50 = ifelse(header %in% blasted_lr_vs_sr$V1[blasted_lr_vs_sr$V3 > 50 & blasted_lr_vs_sr$V13 > 50], 1, 0))



# summary per gene class across all samples
gn_summary1 <- gene_length_lr_2 %>% 
  group_by(class) %>%
  summarise(n_lr = n(), mean_len = mean(length), detect_rate80 = mean(sr80), detect_rate50 = mean(sr50)) %>%
  arrange(desc(mean_len))

gn_summary1 <- gn_summary1 %>% left_join(gene_length_sr %>% group_by(class) %>% summarise(n_sr = n()) %>% ungroup(), by = "class") %>%
  mutate(n_sr =ifelse(is.na(n_sr), 0, n_sr))

# summary per gene class and city across all within-city samples
gn_summary2 <- gene_length_lr_2 %>% 
  group_by(class, City) %>%
  summarise(n_lr = n(), mean_len = mean(length), detect_rate80 = mean(sr80), detect_rate50 = mean(sr50)) %>%
  arrange(desc(mean_len))

gn_summary2 <- gn_summary2 %>% left_join(gene_length_sr %>% group_by(class, City) %>% summarise(n_sr = n()) %>% ungroup(), by = c("class","City")) %>%
  mutate(n_sr =ifelse(is.na(n_sr), 0, n_sr))

# summary per gene class, city, and location across all within-city samples
gn_summary3 <- gene_length_lr_2 %>% 
  group_by(class, Location, City, sample) %>%
  summarise(n_lr = n(), mean_len = mean(length), detect_rate80 = mean(sr80), detect_rate50 = mean(sr50)) %>%
  arrange(desc(mean_len))

gn_summary3 <- gn_summary3 %>% left_join(gene_length_sr %>% group_by(class, Location, City) %>% summarise(n_sr = n()) %>% ungroup(), by = c("class","Location", "City")) %>%
  mutate(n_sr =ifelse(is.na(n_sr), 0, n_sr))

# number of args per class
tbl_data <- gn_summary1 %>% arrange(desc(n_lr)) %>% dplyr::select(class, n_lr, n_sr)

# number of args per class and city
tbl_data2 <- gn_summary2 %>% arrange(class, desc(n_lr)) %>% dplyr::select(class, City, n_lr, n_sr)


p1 <- ggplot(gn_summary1, aes(x = mean_len, y = detect_rate80 , color = !class %in% c("qnr"))) +
  geom_smooth(data = gn_summary1[!gn_summary1$class %in% c("qnr"),],
              method = "lm", se = TRUE, color = "gray40", linetype = "dashed") +
  geom_point(aes(size = n_lr), alpha = 0.7)+# , color = "steelblue") +
  geom_text_repel(aes(label = class), size = 3, color = "black") +
  geom_table(data = data.frame(x = Inf, y = Inf, tbl = I(list(tbl_data))),
             aes(x = x, y = y, label = tbl),
             hjust = 1, vjust = 1, size = 2, table.theme = ttheme_gtlight) +
  scale_size_continuous(range = c(1, 15),
                        breaks = c(10, 50, 100, 500, 1000, 5000),
                        name = "n (long-read genes)") +
  labs(x = "Mean gene length (aa)",
       y = "Detection rate (80 ID / 80 Coverage threshold)",
       title = "Class-level detection rate vs. gene length",
       subtitle = "Line and shaded band: linear model",
       color = "In the regresion model") +
  scale_y_continuous(limits = c(0, .15)) +
  scale_color_manual(values = c(pal_8[7:8])) +
  theme_minimal()

p1



p2 <-  gene_length_lr_2 %>% mutate(sr80 = ifelse(sr80==0, "Not detected","Detected")) %>% 
  ggplot(aes(y = class, x = length, fill = sr80)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.95)) +
  geom_jitter(alpha = 0.1, size = 0.1, color = "black",
              position = position_jitterdodge(jitter.height = 0, jitter.width = 0.15, dodge.width = 0.75)) +
  ylab("ARG class") +
  scale_fill_manual(values = c(pal_8)) +
  xlab("amino acids") +
  ggtitle("Distribution of LR-ARGs length by class and if the gene was detected in SR-seq") + 
  theme_minimal() 

p2 


p4 <-  bind_rows(gene_length_sr %>% dplyr::select(length, class) %>% mutate(read = "short-read ARGs"), gene_length_lr %>% dplyr::select(length, class) %>% mutate(read = "long-read ARGs")) %>%
  ggplot(aes(y = class, x = length, fill = read)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.95)) +
  geom_jitter(alpha = 0.1, size = 0.1, color = "black",
              position = position_jitterdodge(jitter.height = 0, jitter.width = 0.15, dodge.width = 0.75)) +
  ylab("ARG class") +
  scale_fill_manual(values = c(pal_8)) +
  xlab("amino acids") +
  ggtitle("Distribution of gene length by class") + 
  theme_minimal() 

p4 

lr_genes_and_assembly <- gene_length_lr_2 %>%
  mutate(len =  length) %>%  
  left_join(contig_stat[,c("sample","N50","sum_len","Q2","Q1","Q3", "num_seqs", "avg_len")], by = "sample") %>% 
  mutate(lenrat =  len/Q3) 


# test if the long-read ARGs that are also detected in short-read depend on several factors:
# total number of bp in the assembly (sr), average N50 (sr), length of the ARG (lr), ARG class, city, location

lr_genes_and_assembly$class <- relevel(factor(lr_genes_and_assembly$class), ref = "aph(6)")
lr_genes_and_assembly <- lr_genes_and_assembly %>% left_join(contig_stat_lr[,c("sample","N50","sum_len","Q2","Q1","Q3", "num_seqs", "avg_len")], by = "sample")

# GLM
mod_gene_full  <- glmmTMB(sr80 ~ log(len) + log(sum_len.x) + log(Q2.x) + class , family = binomial, data = lr_genes_and_assembly[!lr_genes_and_assembly$class %in% c("qnr"),])
mod_gene_full2  <- glmmTMB(sr80 ~ log(len) + log(sum_len.x) + log(Q2.x) + class + Location, family = binomial, data = lr_genes_and_assembly[!lr_genes_and_assembly$class %in% c("qnr"),])
mod_gene_full3  <- glmmTMB(sr80 ~ log(len) + log(sum_len.x) + log(Q2.x) + class + City, family = binomial, data = lr_genes_and_assembly[!lr_genes_and_assembly$class %in% c("qnr"),])
mod_gene_len   <- glmmTMB(sr80 ~ log(len) + log(sum_len.x) + log(Q2.x) , family = binomial, data = lr_genes_and_assembly[!lr_genes_and_assembly$class %in% c("qnr"),])
mod_gene_sample   <- glmmTMB(sr80 ~ sample , family = binomial, data = lr_genes_and_assembly[!lr_genes_and_assembly$class %in% c("qnr"),])
mod_gene_class <- glmmTMB(sr80 ~ class + log(sum_len.x) + log(Q2.x) , family = binomial, data = lr_genes_and_assembly[!lr_genes_and_assembly$class %in% c("qnr"),])
mod_gene_location   <- glmmTMB(sr80 ~ Location , family = binomial, data = lr_genes_and_assembly[!lr_genes_and_assembly$class %in% c("qnr"),])
mod_gene_location2   <- glmmTMB(sr80 ~ Location + log(sum_len.x) + log(Q2.x), family = binomial, data = lr_genes_and_assembly[!lr_genes_and_assembly$class %in% c("qnr"),])
mod_gene_city   <- glmmTMB(sr80 ~ City , family = binomial, data = lr_genes_and_assembly[!lr_genes_and_assembly$class %in% c("qnr"),])
mod_gene_city2   <- glmmTMB(sr80 ~ City  + log(len) + log(sum_len.x) + log(Q2.x), family = binomial, data = lr_genes_and_assembly[!lr_genes_and_assembly$class %in% c("qnr"),])

summary(mod_gene_full)
summary(mod_gene_full2)

tidy_model <- function(model) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    mutate(
      sig = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        p.value < 0.1   ~ ".",
        TRUE ~ ""
      ),
      across(c(estimate, std.error, statistic), ~ signif(.x, 4)),
      p.value = signif(p.value, 3)
    ) %>%
    dplyr::select(term, estimate, std.error, statistic, p.value, sig)
}


results_mod_gene_full <- tidy_model(mod_gene_full)
results_mod_gene_full2 <- tidy_model(mod_gene_full2)
results_mod_gene_city2 <- tidy_model(mod_gene_city2)
results_mod_gene_location2 <- tidy_model(mod_gene_location2)

#Anova(mod_gene_full) 
#AIC(mod_gene_full, mod_gene_full2, mod_gene_full3) 
#BIC(mod_gene_full, mod_gene_full2, mod_gene_full3) 
#Anova(mod_gene_len)

# mod_gene_rand_sample  <- glmmTMB(sr80 ~ log(len) + (1 | sample) + class, family = binomial, data = lr_genes_and_assembly[!lr_genes_and_assembly$class %in% c("qnr"),])
# mod_gene_rand_len_sample   <- glmmTMB(sr80 ~ log(len) + (1 | sample) , family = binomial, data = lr_genes_and_assembly[!lr_genes_and_assembly$class %in% c("qnr"),])
# mod_gene_rand_class_sample <- glmmTMB(sr80 ~ class + (1 | sample) , family = binomial, data = lr_genes_and_assembly[!lr_genes_and_assembly$class %in% c("qnr"),])
# 
# summary(mod_gene_rand_sample)
# summary(mod_gene_rand_len_sample)
# summary(mod_gene_rand_class_sample)


# Test if the short-read ARGs NOT matching any long-read ARG are biased towards any factor





blasted_sr_vs_lr <- read.table("../resource_generation/12_ARGs_fARGene/sr_contigs/merged-orf-amino_vs_predicted-orfs-amino-nonredundant.tsv", quote="\"", comment.char="")
blasted_sr_vs_lr$sample_sr <- sapply(strsplit(blasted_sr_vs_lr$V1, split = "@@@"), function(x) x[2])
blasted_sr_vs_lr$sample_sr <- gsub("sample", "s", blasted_sr_vs_lr$sample_sr)
blasted_sr_vs_lr$sample_lr <- sapply(strsplit(blasted_sr_vs_lr$V2, split = "@@@"), function(x) x[1])
blasted_sr_vs_lr$class <- sapply(strsplit(blasted_sr_vs_lr$V2, split = "@@@"), function(x) x[2])
blasted_sr_vs_lr$class <- gene_classes$X3[match(blasted_sr_vs_lr$class, gene_classes$X1)]
blasted_sr_vs_lr0 <- blasted_sr_vs_lr

# since we are looking for genes detected in short-read and long-read, the sample lr and sr in blast should match
blasted_sr_vs_lr <- blasted_sr_vs_lr[blasted_sr_vs_lr$sample_lr == blasted_sr_vs_lr$sample_sr,]
blasted_sr_vs_lr <- blasted_sr_vs_lr %>% group_by(V1) %>% slice_head(n = 1) %>% ungroup()
dim(blasted_sr_vs_lr)


blasted_sr_vs_lr0 %>% filter(!V1 %in% blasted_sr_vs_lr$V1) %>% group_by(V1) %>% slice_head(n = 1) %>% ungroup()  %>% group_by(class) %>% summarise(n = n())
unmatched_sample <- blasted_sr_vs_lr0 %>% filter(!V1 %in% blasted_sr_vs_lr$V1) %>% group_by(V1) %>% slice_head(n = 1)


sr_genes$length <- gene_length_sr$length[match(sr_genes$orf, gene_length_sr$header)]

sr_genes$undetected <- ifelse(!sr_genes$orf %in% blasted_sr_vs_lr$V1, 1, 0)

sr_genes$low_quality <- ifelse(sr_genes$orf %in% 
                        blasted_sr_vs_lr$V1[blasted_sr_vs_lr$V3 < 80 | blasted_sr_vs_lr$V13 < 80 |
                        blasted_sr_vs_lr$V5/blasted_sr_vs_lr$V4>.05],1, 0)

sr_genes$high_quality <- ifelse(sr_genes$orf %in% 
                         blasted_sr_vs_lr$V1[
                          blasted_sr_vs_lr$V3 >= 80 & blasted_sr_vs_lr$V13 >= 80 &
                          blasted_sr_vs_lr$V5/blasted_sr_vs_lr$V4 <= 0.05], 1, 0)

sr_genes$highest_quality <- ifelse(sr_genes$orf %in% 
                            blasted_sr_vs_lr$V1[
                              blasted_sr_vs_lr$V3 >= 99 &  blasted_sr_vs_lr$V13>=99 & 
                              blasted_sr_vs_lr$V5/blasted_sr_vs_lr$V4 <= 0.05], 1, 0)

sr_genes_summary <- sr_genes %>% mutate(h = high_quality | highest_quality, 
                    l = low_quality ) %>% ungroup() %>% 
  group_by(class) %>% summarise(n = n(), h = sum(h), l = sum(l), u = sum(undetected), mean_len = mean(length)) %>% 
  mutate(sh = h /n, sl = l/n, su = u/n)



sr_genes <- sr_genes %>% mutate(label = ifelse(highest_quality | high_quality, "H",
                                           ifelse(low_quality, "L", "U")))

sr_genes_summary2 <- sr_genes %>% 
  group_by(class) %>% mutate(N = n()) %>% ungroup() %>% 
  group_by(class, label) %>% summarise(N = N[1], n = n(), mean_len = mean(length)) %>% ungroup() %>%
  mutate(r = n/N)


p5.1 <- ggplot(sr_genes_summary2 %>% mutate(label = ifelse(label == "H", "High", ifelse(label == "L", "low", "unaligned"))) , 
             aes(x = r, y = class , fill = forcats::fct_rev(label))) +
  geom_col() +
  xlab("Rate") +
  ylab("") +
  labs(fill = "") +
  scale_fill_manual(values = c(pal_8[5:8])) +
  theme_minimal()

class_totals <- sr_genes_summary2 %>%
  group_by(class) %>%
  summarise(total_n = sum(n))

p5.2 <- ggplot(sr_genes_summary2 %>% mutate(label = ifelse(label == "H", "High", ifelse(label == "L", "low", "unaligned"))) , 
               aes(x = n, y = class , fill = forcats::fct_rev(label))) +
  geom_col() +
  geom_text(data = class_totals, aes(x = total_n, y = class, label = total_n),
            inherit.aes = FALSE, hjust = -0.1, size = 3) +
  xlab("ARGs") +
  ylab("") +
  labs(fill = "") +
  scale_fill_manual(values = c(pal_8[5:8])) +
  scale_x_continuous(limits = c(0, 4000)) +
  theme_minimal()


p5 <- grid.arrange(p5.1 , p5.2 + theme(legend.position = "none", axis.text.y = element_blank()), nrow = 1) 

sr_genes_summary3 <- sr_genes_summary2 %>%
  pivot_wider(
    id_cols = class,
    names_from = label,
    values_from = c(n, mean_len, r) )

p6 <-  bind_rows(sr_genes %>% mutate(label = ifelse(label == "H", "High", ifelse(label == "L", "low", "unaligned")))) %>%
  ggplot(aes(y = class, x = length, fill = label)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.95)) +
  geom_jitter(alpha = 0.1, size = 0.1, color = "black",
              position = position_jitterdodge(jitter.height = 0, jitter.width = 0.15, dodge.width = 0.75)) +
  ylab("ARG class") +
  scale_fill_manual(values = c(pal_8[5:8])) +
  xlab("amino acids") +
  labs(fill="") + 
  ggtitle("Distribution of gene length by class and alignment to LR quality / unaligned") + 
  theme_minimal() 


dev.off()

svg("supp9/Detection_LR_genes_in_SR_sequencing.svg", width = 10, height = 8)
p1
dev.off()

svg("supp9/gene_length_LR_by_detection_in_SR.svg", width = 10, height = 8)
p2
dev.off()

svg("supp9/gene_length_LR_SR.svg", width = 10, height = 8)
p4
dev.off()

svg("supp9/gene_length_SR_by_alignment_to_LR_quality.svg", width = 10, height = 8)
p6
dev.off()

svg("supp9/SR_gene_detection_undetection_rates_and_number_of_ARGs.svg", width = 10, height = 8)
grid.arrange(p5.1 , p5.2 + theme(legend.position = "none", axis.text.y = element_blank()), nrow = 1) 
dev.off()


results_mod_gene_full <- tidy_model(mod_gene_full)
results_mod_gene_full2 <- tidy_model(mod_gene_full2)
results_mod_gene_city2 <- tidy_model(mod_gene_city2)
results_mod_gene_location2 <- tidy_model(mod_gene_location2)

write.csv(results_mod_gene_full, "supp9/DR-genelen-depth-q2-class.csv", row.names = FALSE)
write.csv(results_mod_gene_full2, "supp9/DR-genelen-depth-q2-class-location.csv", row.names = FALSE)
write.csv(results_mod_gene_city2, "supp9/DR-genelen-depth-q2-city.csv", row.names = FALSE)
write.csv(results_mod_gene_location2, "supp9/DR-genelen-depth-q2-location.csv", row.names = FALSE)



"number_ARGs_per_sample.csv"

write.csv( sr_genes %>% group_by(City, Location, sample) %>% summarise(n=n_distinct(centroid)), "supp9/number_ARGs_per_sample_SR.csv")

m <- sr_genes %>% mutate(sum_len = contig_stat$sum_len[match(sample, contig_stat$sample)],
                    Q2 = contig_stat$Q2[match(sample, contig_stat$sample)]) %>% 
  group_by(City, Location, sample) %>% summarise(n=n_distinct(centroid), Q2 = Q2[1], sum_len = sum_len[1])

plot(m$Q2, m$n)
plot(log(m$sum_len), m$n)
ggplot(m, aes(x = log(sum_len), y = n, color = City)) +
  geom_point()

ggplot(m, aes(x = log(sum_len), y = n, color = City)) +
  geom_point()

ggplot(m, aes(x = Q2, y = n, color = City)) +
  geom_point()

ggplot(m, aes(x = log(sum_len), y = Q2, color = City)) +
  geom_point()


mod_sr  <- lm(n ~ log(sum_len) + log(Q2), data = m)
summary(mod_sr)

library(plotly)
city_levels <- levels(factor(m$City))

m$n_ARGs <- m$n
plot_ly(m, x = ~log(sum_len), y = ~log(Q2), z = ~n_ARGs, color = ~City,
        colors = pal_8[1:length(city_levels)],
        type = "scatter3d", mode = "markers", marker = list(size = 4)) %>%
  plotly::layout(scene = list(
    xaxis = list(title = "log(sum_len)"),
    yaxis = list(title = "log(Q2)"),
    zaxis = list(title = "n")
  ))



library(scatterplot3d)

city_levels <- levels(factor(m$City))
city_colors <- pal_8[1:length(city_levels)]
colors <- city_colors[as.numeric(factor(m$City))]
s3d <- scatterplot3d(log(m$sum_len), log(m$Q2), m$n,
                     color = colors, pch = 19, angle = 20,
                     xlab = "log(sum_len)", ylab = "log(Q2)", zlab = "n")
legend("topright", legend = levels(factor(m$City)),
       col = 1:length(unique(m$City)), pch = 19)

results_sr_model <- tidy_model(mod_sr)
write.csv(results_sr_model, "supp9/nARGS_SR-depth-q2.csv", row.names = FALSE)
