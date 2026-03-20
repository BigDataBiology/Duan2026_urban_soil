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

setwd("~/Documents/GitHub/urban_soil/analysis/")

gene_classes <- data.frame(rbind(
  c("aac2p",	"aminoglycoside acetyltransferase",	"aac(2')", "aminoglycoside"),
  c("aac3_1",	"aminoglycoside acetyltransferase",	"aac(3)", "aminoglycoside"), c("aac3_2",	"aminoglycoside acetyltransferase",	"aac(3)", "aminoglycoside"),
  c("aac6p_1",	"aminoglycoside acetyltransferase",	"aac(6')", "aminoglycoside"), c("aac6p_2",	"aminoglycoside acetyltransferase",	"aac(6')", "aminoglycoside"),
  c("aac6p_3",	"aminoglycoside acetyltransferase",	"aac(6')", "aminoglycoside"), c("aph2b",	"aminoglycoside phosphotransferase", "aph(2'')", "aminoglycoside"),
  c("aph3p",	"aminoglycoside phosphotransferase", "aph(3')", "aminoglycoside"), c("aph6p",	"aminoglycoside phosphotransferase", "aph(6)", "aminoglycoside"),
  c("class_a",	"beta-lactamase A",	"beta-lactamase A", "beta-lactamase"), c("class_b1_b2",	"beta-lactamase B",	"beta-lactamase B1B2", "beta-lactamase"),
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


# font for plots 
general_size <- 10

# palette for plots  
pal_8 <- brewer.pal(8, "Dark2")
pal_8 <- pal_8[c(1,6,2:5,7,8)]
pal_div <- brewer.pal(8, "BrBG")
pal_seq <- brewer.pal(8, "YlOrBr")

##########

genes <- genes %>% 
  mutate(
    hclass = 
      factor(hclass, levels = c("aminoglycoside", "beta-lactamase", 
                                "macrolide", "tetracycline", "quinolone")))
##########################################################################################
##########################################################################################
##########################################################################################


plot_genes_per_class_contig <- genes %>%  
  group_by(City, centroid) %>% slice_head(n = 1) %>% 
  group_by(centroid) %>% slice_head(n = 1) %>% 
  ggplot(aes(x = class, fill = hclass, pattern = City)) +
  geom_bar_pattern( position = position_dodge2(preserve = "single", width = 0.3, padding = 0.1), 
                    width = 0.85,
                    pattern_density = 0.15,
                    pattern_color = "black",
                    pattern_fill = "black",
                    pattern_size =  0.3) +
  scale_fill_manual(values = pal_8) +
  scale_pattern_manual(values = c('Shanghai' = 'circle', 'Nanjing' = 'stripe')) + 
  theme_minimal() +
  labs( fill = "", pattern = "") +
  facet_wrap( hclass ~ ., scales = "free", nrow = 1) +
  xlab("") +
  ylab("Number of ARGs") +
  scale_x_discrete(labels = function(x) {
    x <- gsub("beta-lactamase", "", x)
    x <- gsub("tet", "", x)
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    
    sapply(x, function(lbl) {
      if (grepl("aph", lbl) | grepl("aac", lbl) | 
          grepl("erm", lbl) | grepl("mph", lbl) |
          grepl("qnr", lbl)) {
        bquote(italic(.(lbl)))
      } else {
        lbl
      }
    })
    
    }) + 
  guides(
    fill = "none",
    pattern = guide_legend(override.aes = list(fill = pal_8[8]))
  ) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = general_size),
    title = element_text(size = general_size + 2, face = "bold"),
    strip.text = element_text(size = general_size , face = "bold"),
    axis.title = element_text(size = general_size , face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.1, size = general_size),
    axis.text.y = element_text(size = general_size),
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.background = element_blank())

plot_genes_per_class_contig

ggsave("supp8/s8.svg", plot_genes_per_class_contig, width = 180, height = 70, unit = "mm")


