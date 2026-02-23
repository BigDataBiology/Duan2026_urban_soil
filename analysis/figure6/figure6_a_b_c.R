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

setwd("~/Documents/GitHub/urban_soil/analysis/")

gene_classes <- data.frame(rbind(
  c("aac2p",	"aminoglycoside acetyltransferase",	"aac", "aminoglycoside"),
  c("aac3_1",	"aminoglycoside acetyltransferase",	"aac3", "aminoglycoside"), c("aac3_2",	"aminoglycoside acetyltransferase",	"aac3", "aminoglycoside"),
  c("aac6p_1",	"aminoglycoside acetyltransferase",	"aac6'", "aminoglycoside"), c("aac6p_2",	"aminoglycoside acetyltransferase",	"aac6'", "aminoglycoside"),
  c("aac6p_3",	"aminoglycoside acetyltransferase",	"aac6'", "aminoglycoside"), c("aph2b",	"aminoglycoside phosphotransferase", "aph2''", "aminoglycoside"),
  c("aph3p",	"aminoglycoside phosphotransferase", "aph3'", "aminoglycoside"), c("aph6p",	"aminoglycoside phosphotransferase", "aph6'", "aminoglycoside"),
  c("class_a",	"beta-lactamase A",	"beta-lactamase A", "beta-lactamase"), c("class_b1_b2",	"beta-lactamase B",	"beta-lactamase B1B2", "beta-lactamase"),
  c("class_b3",	"beta-lactamase B",	"beta-lactamase B3", "beta-lactamase"), c("class_c",	"beta-lactamase C",	"beta-lactamase C", "beta-lactamase"),
  c("class_d1",	"beta-lactamase D",	"beta-lactamase D", "beta-lactamase"), c("class_d2",	"beta-lactamase D",	"beta-lactamase D", "beta-lactamase"),
  c("erm_1",	"macrolide erm", "erm", "macrolide"), c("erm_2",	"macrolide erm", "erm",	'macrolide'),
  c("mph",	"macrolide mph", "mph",	"macrolide"), c("qnr",	"qrn", "qnr", "qnr"), c("tet_efflux",	"tetracycline efflux pump", "tet efflux", "tetracycline"),
  c("tet_enzyme",	"tetracycline inactivating enzymes", "tet enzyme", "tetracycline"), c("tet_rpg",	"tetracycline ribosomal protection genes", "tet rpg", "tetracycline")))

metadata <- readxl::read_excel("../resource_generation/12_ARGs_fARGene/supp1_sample_metadata.xlsx")
metadata <- metadata %>% mutate(Sample_id = gsub("ample", "", Sample_id)) 

mag_contig <- read.table("../resource_generation/12_ARGs_fARGene/assembly_info_mag.tsv", quote="\"", comment.char="")

genes <- read.table("../resource_generation/12_ARGs_fARGene/Assembly/predicted_gene_names.txt", quote="\"", comment.char="")
mags <- read.table("../resource_generation/12_ARGs_fARGene/Genes/mags/predicted_gene_names.txt", quote="\"", comment.char="")
clusters <- read.table("../resource_generation/12_ARGs_fARGene/Assembly/clusters.uc", quote="\"", comment.char="")
clusters2 <- read.table("../resource_generation/12_ARGs_fARGene/Genes/mags/closest_hits.tsv", quote="\"", comment.char="")

genes <- genes %>% mutate(centroid = clusters$V10[match(V1, clusters$V9)])
genes <- genes %>% mutate(centroid = ifelse(centroid == "*", V1, centroid))
mags <- mags %>% mutate(centroid = clusters2$V2[match(V1, clusters2$V1)])


## contigs
split1 <- strsplit(genes$V1, split = "@@@")
split2 <- strsplit(sapply(split1, function(x) x[2]), split = ".fa.gz-")
split3 <- strsplit(sapply(split2, function(x) x[2]), split = ".fasta")
split4 <- strsplit(sapply(split1, function(x) x[2]), split = "_")

genes <- genes %>% mutate(contig = sapply(split1, function(x) x[1]))
genes <- genes %>% mutate(sample = sapply(split4, function(x) x[1]))
genes <- genes %>% mutate(hmm = sapply(split3, function(x) x[1]))

split1 <- strsplit(genes$centroid, split = "@@@")
split2 <- strsplit(sapply(split1, function(x) x[2]), split = ".fa.gz-")
split3 <- strsplit(sapply(split2, function(x) x[2]), split = ".fasta")
split4 <- strsplit(sapply(split1, function(x) x[2]), split = "_")

genes <- genes %>% mutate(c_contig = sapply(split1, function(x) x[1]))
genes <- genes %>% mutate(c_sample = sapply(split4, function(x) x[1]))
genes <- genes %>% mutate(c_hmm = sapply(split3, function(x) x[1]))
genes <- genes %>% mutate(data_type = "assembly")
genes <- genes %>% mutate(id = NA)
genes <- genes %>% mutate(c_id = NA)



## MAGs
split1 <- strsplit(mags$V1, split = "@@@")
split2 <- strsplit(sapply(split1, function(x) x[2]), split = ".faa.gz-")
split3 <- strsplit(sapply(split2, function(x) x[2]), split = ".fasta")
split4 <- strsplit(sapply(split1, function(x) x[1]), split = "_")

mags <- mags %>% mutate(contig = sapply(split4, function(x) paste(x[3],x[4],x[5], sep = "_")))
mags <- mags %>% mutate(sample = sapply(split4, function(x) x[3]))
mags <- mags %>% mutate(hmm = sapply(split3, function(x) x[1]))
mags <- mags %>% mutate(id = sapply(split4, function(x) x[length(x)]))

split1 <- strsplit(mags$centroid, split = "@@@")
split2 <- strsplit(sapply(split1, function(x) x[2]), split = ".faa.gz-")
split3 <- strsplit(sapply(split2, function(x) x[2]), split = ".fasta")
split4 <- strsplit(sapply(split1, function(x) x[1]), split = "_")

mags <- mags %>% mutate(c_contig = sapply(split4, function(x) paste(x[3],x[4],x[5], sep = "_")))
mags <- mags %>% mutate(c_sample = sapply(split4, function(x) x[3]))
mags <- mags %>% mutate(c_hmm = sapply(split3, function(x) x[1]))
mags <- mags %>% mutate(data_type = "mags")
mags <- mags %>% mutate(c_id = sapply(split4, function(x) x[length(x)]))


mags <- mags %>% mutate(id = mag_contig$V9[match(contig, mag_contig$V1)],
                        c_id = mag_contig$V9[match(c_contig, mag_contig$V1)])

rm(mag_contig, split1, split2, split3, split4)

mags <- mags %>% mutate(sample = gsub("ample", "", sample)) 
mags <- mags %>% mutate(c_sample = gsub("ample", "", c_sample)) 
mags <- mags %>% mutate(read = sapply(strsplit(V1,split = "@@@"), function(x) x[1])) 


######################################
######################################
# 17 contigs with repeated centroid:
# 21 contigs with different hmm than the centroid c_hmm, all erm1 or erm2
# remove contigs with double hmm, in this case they are the same as the ones with duplicated centroid and different centroid class to the read

genes <- genes %>%filter(hmm == c_hmm)
mags <- mags %>% filter(hmm == c_hmm)

######################################
######################################

## add metadata 

mags <- mags %>% mutate(description = gene_classes$X2[match(c_hmm, gene_classes$X1)],
                        class = gene_classes$X3[match(c_hmm, gene_classes$X1)],
                        hclass = gene_classes$X4[match(c_hmm, gene_classes$X1)],
                        City = metadata$City[match(sample, metadata$Sample_id)],
                        Location = metadata$Location[match(sample, metadata$Sample_id)],
                        Date = metadata$Date[match(sample, metadata$Sample_id)],
                        Longitude = metadata$Date[match(sample, metadata$Longitude)],
                        Latitude = metadata$Date[match(sample, metadata$Latitude)])


genes <- genes %>% mutate(description = gene_classes$X2[match(c_hmm, gene_classes$X1)],
                          class = gene_classes$X3[match(c_hmm, gene_classes$X1)],
                          hclass = gene_classes$X4[match(c_hmm, gene_classes$X1)],
                          City = metadata$City[match(sample, metadata$Sample_id)],
                          Location = metadata$Location[match(sample, metadata$Sample_id)],
                          Date = metadata$Date[match(sample, metadata$Sample_id)],
                          Longitude = metadata$Date[match(sample, metadata$Longitude)],
                          Latitude = metadata$Date[match(sample, metadata$Latitude)])

genes <- genes %>% mutate(read = paste(contig, sample, sep = "_")) 


## remove sample s4 and s11

genes <- genes %>% 
  filter(!sample %in% c("s4", "s11"))

mags <- mags %>% 
  filter(!sample %in% c("s4", "s11"))

# font for plots 
general_size <- 10
# palette for plots  
pal_8 <- brewer.pal(8, "Dark2")
pal_8 <- pal_8[c(1,6,2:5,7,8)]
pal_div <- brewer.pal(8, "BrBG")
pal_seq <- brewer.pal(8, "YlOrBr")

## Calculate the pairwise sets
sets <- genes %>% 
  ungroup()  %>% 
  arrange(centroid) %>% 
  group_by(centroid) %>% 
  mutate(n_sample = n_distinct(sample)) %>% 
  mutate(single = (n_sample ==1))  %>%  
  group_by(sample) %>%
  summarise(centroids = list(centroid), .groups = "drop") %>% # put every query in a list
  mutate(Longitude = metadata$Longitude[match(sample, metadata$Sample_id)],
         Latitude = metadata$Latitude[match(sample, metadata$Sample_id)]) 

## calculate jaccard's index
JI <- expand_grid(sample_ref = sets$sample, sample_comp = sets$sample)  %>%
  left_join(sets, by = c("sample_ref" = "sample")) %>%
  rename(values1 = centroids) %>%
  left_join(sets, by = c("sample_comp" = "sample")) %>%
  rename(values2 = centroids) %>%
  mutate(jaccard = map2_dbl(values1, values2, ~ length(intersect(.x, .y)) / length(union(.x, .y)))) %>% 
  mutate(n_overlap = map2_dbl(values1, values2, ~ length(intersect(.x, .y)))) %>% 
  mutate(n_ref = map2_dbl(values1, values2, ~ length(.x))) %>% 
  mutate(n_comp = map2_dbl(values1, values2, ~ length(.y))) %>% 
  group_by(sample_ref, sample_comp) %>% 
  mutate(distance = distHaversine(c(Longitude.x[1], Latitude.x), c(Longitude.y, Latitude.y))) %>%
  select(sample_ref, sample_comp, jaccard, distance, n_ref,  n_comp, n_overlap)  %>% 
  filter(sample_ref != sample_comp)


## transform the jaccard index to matrix to use complex::heatmap
mat <- xtabs(jaccard ~ sample_ref + sample_comp, data = JI)
mat <- mat[,gtools::mixedsort(colnames(mat))]
mat <- mat[gtools::mixedsort(rownames(mat)),]
mat <- mat[metadata$Sample_id[metadata$Sample_id %in% genes$sample], 
           metadata$Sample_id[metadata$Sample_id %in% genes$sample]]

cities <- sort(unique(metadata$City))

# annotation for heatmap
col_ann = HeatmapAnnotation(
  City = metadata$City[match(colnames(mat), metadata$Sample_id) ],
  col = list(City = setNames(pal_8[1:2], cities)))

row_ann = rowAnnotation(
  City     = metadata$City[match(rownames(mat), metadata$Sample_id)],
  col = list(City = setNames(pal_8[1:2], cities)),
  show_legend = FALSE)


# heatmap Jaccard's index
pheat11 <- ComplexHeatmap::Heatmap(
  mat, 
  rect_gp = gpar(type = "none"), column_dend_side = "top",
  cell_fun = function(j, i, x, y, w, h, fill) {
    if(as.numeric(x) >= 1 - as.numeric(y)) {
      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
    }},
  show_row_dend = FALSE,
  name = "J.I.", right_annotation = row_ann,
  top_annotation = col_ann, col = pal_seq, na_col = "white",
  show_row_names = FALSE, show_column_names = FALSE)




##########################################################################################
##########################################################################################
##########################################################################################v


genes <- genes %>% 
  mutate(
    hclass = 
      factor(hclass, levels = c("aminoglycoside", "beta-lactamase", 
                                "macrolide", "tetracycline", "qnr")))

mags <- mags %>% 
  mutate(
    hclass = 
      factor(hclass, levels = c("aminoglycoside", "beta-lactamase", 
                                "macrolide", "tetracycline", "qnr")))



##########################################################################################
##########################################################################################
##########################################################################################


two_samples <- genes %>% select(data_type, centroid, City, Location, sample) %>% 
  group_by(centroid) %>% 
  mutate(N = n_distinct(sample)) %>% 
  filter(N > 1) %>% ungroup() %>% select(centroid) %>% distinct() %>% 
  pull() 

two_locations <- genes %>% select(data_type, centroid, City, Location, sample) %>% 
  group_by(centroid) %>% 
  mutate(N = n_distinct(Location)) %>% 
  filter(N > 1) %>% ungroup() %>% select(centroid) %>% distinct() %>% 
  pull() 

count_locations <- rbind(
  genes %>% 
    filter(centroid %in% two_locations) %>% 
    select(data_type, centroid, City, Location, sample) %>% 
    group_by(centroid) %>% 
    summarise(N = n_distinct(City)) %>% 
    ungroup() %>% 
    group_by(N) %>% 
    summarise(n = n()) %>% 
    mutate(gr = ifelse(N == 1, "Single city", "Multi city"),
           lev = "City"),
  
  genes %>% 
    filter(centroid %in% two_samples) %>% 
    select(data_type, centroid, City, Location, sample) %>% 
    group_by(centroid) %>% 
    summarise(N = n_distinct(Location)) %>% 
    ungroup() %>% 
    group_by(N) %>% 
    summarise(n = n()) %>% 
    mutate(gr = ifelse(N == 1, "Single location", "Multi location"),
           lev = "location"),
  
  
  genes %>% select(data_type, centroid, City, Location, sample) %>% 
    group_by(centroid) %>% 
    summarise(N = n_distinct(sample)) %>% 
    ungroup() %>% 
    group_by(N) %>% 
    summarise(n = n()) %>% 
    mutate(gr = ifelse(N == 1, "Single sample", "Multi sample"),
           lev = "sample"))


count_locations <- count_locations %>%
  ungroup() %>% 
  group_by(gr, lev) %>% 
  summarise(n = sum(n))



two_mags_mags <- mags %>% select(data_type, centroid, City, Location, sample, id) %>% 
  group_by(centroid) %>% 
  mutate(N = n_distinct(id)) %>% 
  filter(N > 1) %>% ungroup() %>% select(centroid) %>% distinct() %>% 
  pull() 

two_samples_mags <- mags %>% select(data_type, centroid, City, Location, sample) %>% 
  group_by(centroid) %>% 
  mutate(N = n_distinct(sample)) %>% 
  filter(N > 1) %>% ungroup() %>% select(centroid) %>% distinct() %>% 
  pull() 

two_locations_mags <- mags %>% select(data_type, centroid, City, Location, sample) %>% 
  group_by(centroid) %>% 
  mutate(N = n_distinct(Location)) %>% 
  filter(N > 1) %>% ungroup() %>% select(centroid) %>% distinct() %>% 
  pull() 

count_locations_mag <- rbind(
  mags %>% 
    filter(centroid %in% two_locations_mags) %>% 
    select(data_type, centroid, City, Location, sample) %>% 
    group_by(centroid) %>% 
    summarise(N = n_distinct(City)) %>% 
    ungroup() %>% 
    group_by(N) %>% 
    summarise(n = n()) %>% 
    mutate(gr = ifelse(N == 1, "Single city", "Multi city"),
           lev = "City"),
  
  mags %>% 
    filter(centroid %in% two_samples_mags) %>% 
    select(data_type, centroid, City, Location, sample) %>% 
    group_by(centroid) %>% 
    summarise(N = n_distinct(Location)) %>% 
    ungroup() %>% 
    group_by(N) %>% 
    summarise(n = n()) %>% 
    mutate(gr = ifelse(N == 1, "Single location", "Multi location"),
           lev = "location"),
  
  
  mags %>% select(data_type, centroid, City, Location, sample) %>% 
    filter(centroid %in% two_mags_mags) %>% 
    group_by(centroid) %>% 
    summarise(N = n_distinct(sample)) %>% 
    ungroup() %>% 
    group_by(N) %>% 
    summarise(n = n()) %>% 
    mutate(gr = ifelse(N == 1, "Single sample", "Multi sample"),
           lev = "sample"),
  
  mags %>% select(data_type, centroid, City, Location, sample, id) %>% 
    group_by(centroid) %>% 
    summarise(N = n_distinct(id)) %>% 
    ungroup() %>% 
    group_by(N) %>% 
    summarise(n = n()) %>% 
    mutate(gr = ifelse(N == 1, "Single MAG", "Multi MAGs"),
           lev = "mag"))


count_locations_mag <- count_locations_mag %>%
  ungroup() %>% 
  group_by(gr, lev) %>% 
  summarise(n = sum(n))


count_locations_mag <- count_locations_mag %>%
  mutate(lev = factor(lev, levels = c("mag", "sample", "location", "City"))) %>%
  ungroup() %>% 
  group_by(lev) %>%
  mutate(p = n / sum(n))


count_locations <- count_locations %>%
  mutate(lev = factor(lev, levels = c("mag", "sample", "location", "City"))) %>% 
  ungroup() %>% 
  group_by(lev) %>%
  mutate(p = n / sum(n))

count_locations <- count_locations %>% 
  mutate(gr2 = ifelse(grepl("Multi", gr), "2", "1"))

count_locations_mag <- count_locations_mag %>% 
  mutate(gr2 = ifelse(grepl("Multi", gr), "2", "1"))


p_circle1 <- ggplot(count_locations, 
                    aes(x0 = 1 , y0 = 0, fill = gr2)) +
  geom_circle(aes(r = p), alpha = 1, show.legend = F) +
  scale_fill_manual(values = pal_8[3:8]) +
  geom_text(aes(x = 1 - p, y = p - 0.1, label = paste0(round((p) * 100), "%")), 
            size = (general_size-2) / ggplot2::.pt, vjust = 0, hjust = 0) +
  geom_text(aes(x = 1 - p, y = p - 0.1, label = paste0(n)), 
            size = (general_size-2) / ggplot2::.pt, vjust = 2.5, hjust = -0.5) +
  coord_equal() +
  theme_void() +
  labs(fill = "") + 
  guides(fill = guide_legend(ncol = 1)) +
  facet_wrap(lev ~ gr, ncol = 2) +
  theme(
    strip.text = element_text(size = general_size-2),
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y  = element_blank(),
    legend.text = element_text(size = general_size)) 



p_circle2 <- ggplot(count_locations_mag, 
                    aes(x0 = 1 , y0 = 0, fill = gr2)) +
  geom_circle(aes(r = p), alpha = 1, show.legend = F) +
  scale_fill_manual(values = pal_8[3:8]) +
  geom_text(aes(x = 1 - p, y = p - 0.1, label = paste0(round((p) * 100), "%")), 
            size = (general_size-2) / ggplot2::.pt, vjust = 0, hjust = 0) +
  geom_text(aes(x = 1 - p, y = p - 0.1, label = paste0(n)), 
            size = (general_size-2) / ggplot2::.pt, vjust = 2.5, hjust = -0.5) +
  coord_equal() +
  theme_void() +
  labs(fill = "") + 
  guides(fill = guide_legend(ncol = 1)) +
  facet_wrap(lev ~ gr, ncol = 2) +
  theme(
    strip.text = element_text(size = general_size-2),
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y  = element_blank(),
    legend.text = element_text(size = general_size)) 



counts_genes_geo_loc <- genes %>% select(data_type, centroid, City, Location, sample) %>% 
  group_by(data_type, City, centroid, Location, sample) %>%
  slice_head(n = 1) %>% 
  ungroup() %>%
  group_by(centroid) %>% 
  mutate(n_city = n_distinct(City)) %>% 
  ungroup() %>% 
  group_by(centroid, City) %>% 
  mutate(n_location = n_distinct(Location)) %>% 
  ungroup() %>% 
  group_by(centroid, City, Location) %>% 
  mutate(n_sample = n_distinct(sample)) %>% 
  ungroup() %>% 
  mutate(label = ifelse(n_city > 1, "Multiple cities",
                        ifelse(n_location > 1, "Multiple locations",
                               ifelse(n_sample > 1, "Multiple samples", "Single sample")))) 



genes_per_sample <- genes %>% 
  group_by(City, Location, sample) %>% 
  summarise(n = n_distinct(centroid)) 

wilcox.test(genes_per_sample$n[genes_per_sample$City == "Nanjing"],
            genes_per_sample$n[genes_per_sample$City != "Nanjing"], "two.sided", paired = F)


genes_per_sample_mags <- mags %>% 
  group_by(City, Location, sample) %>% 
  summarise(n = n_distinct(centroid)) 

wilcox.test(genes_per_sample_mags$n[genes_per_sample$City == "Nanjing"],
            genes_per_sample_mags$n[genes_per_sample$City != "Nanjing"], "two.sided", paired = F)


genes_per_city <- genes %>% group_by(City) %>% summarise(n = n_distinct(centroid)) %>% 
  mutate(City = paste(" ", City)) 

qantiles <- genes %>% 
  group_by(City, Location, sample) %>% 
  summarise(n = n_distinct(centroid)) %>%
  ungroup() %>% 
  group_by(City) %>% 
  summarise(q25 = quantile(n, 0.25), q75 = quantile(n, 0.75))


box2<- genes %>% 
  group_by(City, Location, sample) %>% 
  summarise(n = n_distinct(centroid)) %>% 
  ggplot(aes(x = City, y = n*10, fill = City)) +
  geom_rect(data = genes_per_city,
            aes(xmin = 0, xmax = as.numeric(factor(City)) + 0.4, 
                ymin = 0, ymax = n, fill = City), alpha = 0.3,
            inherit.aes = FALSE, color = "black", linewidth =  0.2,
            show.legend = F) +
  geom_rect(data = qantiles,
            aes(xmin = as.numeric(factor(City)) + 2, 
                xmax = 5, 
                ymin = q25*10, ymax = q75*10, fill = City), alpha = 0.3,
            inherit.aes = FALSE, color = "black", linewidth =  0.2,
            show.legend = F) +
  geom_col(data = genes_per_city, aes(x = City, y = n, fill = City),  
           width = 0.8, color = "black", show.legend = F) +
  geom_boxplot(show.legend = F) +
  geom_jitter(color = "black", height = 0, width = 0.3, show.legend = F) + 
  xlab("") + 
  scale_fill_manual(values = pal_8[c(1,2,1,2)], labels = function(x) gsub(" ", "\n", x))  +
  theme_minimal() + 
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x}) + 
  scale_y_continuous(name = "ARGs per city",
                     sec.axis = sec_axis(~ . / 10, name = "ARGs per sample")) +
  theme(
    legend.position = "bottom",
    axis.text.x.top = element_blank(),
    strip.text = element_text(size = general_size, face = "bold"),
    legend.text = element_text(size = general_size),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = general_size, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank())



ggsave("figure6/main_args_1b.svg", p_circle1 +
         theme(plot.margin = margin(0, 0, 0, 0), panel.spacing = unit(5, "mm")) , width = 90, height = 100, unit = "mm")

ggsave("figure6/main_args_1a.svg", box2 + 
         theme(axis.text.x = element_blank(),
               plot.margin = margin(0, 0, 0, 0), panel.spacing = unit(5, "mm")) , width = 90, height = 100, unit = "mm")



dev.off()
svg("figure6/main_args_1c.svg", width = 70/25.4, height = 100/25.4)
draw(pheat11, annotation_legend_side = "bottom",  merge_legends = FALSE)
dev.off()

#+
#  theme_minimal()

##########################################################################################
##########################################################################################
##########################################################################################
