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
                                "macrolide", "tetracycline", "qnr")))



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


set.seed(2026)
JI <- JI[sample(1:dim(JI)[1], replace = F),]

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
  col = list(City = setNames(pal_8[1:2], cities)),
  annotation_name_side = NULL,
  show_annotation_name = FALSE)

row_ann = rowAnnotation(
  City     = metadata$City[match(rownames(mat), metadata$Sample_id)],
  col = list(City = setNames(pal_8[1:2], cities)),
  show_legend = FALSE,
  annotation_name_side = NULL,
  show_annotation_name = FALSE)



col_ann_nolegend = HeatmapAnnotation(
  City = metadata$City[match(colnames(mat), metadata$Sample_id) ],
  col = list(City = setNames(pal_8[1:2], cities)),
  annotation_name_side = NULL,
  show_annotation_name = FALSE,
  show_legend = FALSE)

row_ann_nolegend = rowAnnotation(
  City     = metadata$City[match(rownames(mat), metadata$Sample_id)],
  col = list(City = setNames(pal_8[1:2], cities)),
  show_legend = FALSE,
  annotation_name_side = NULL,
  show_annotation_name = FALSE)

# heatmap Jaccard's index
pheat11 <- ComplexHeatmap::Heatmap(
  mat, cluster_rows = T, cluster_columns = T,
  rect_gp = gpar(type = "none"), column_dend_side = "top",
  cell_fun = function(j, i, x, y, w, h, fill) {
    if(as.numeric(x) >= 1 - as.numeric(y)) {
      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
    }},
  show_row_dend = FALSE,
  name = "J.I.", right_annotation = row_ann,
  top_annotation = col_ann, col = pal_seq, na_col = "white",
  show_row_names = FALSE, show_column_names = FALSE,
  heatmap_legend_param = list(direction = "horizontal"))


pheat12 <- ComplexHeatmap::Heatmap(
  mat, cluster_rows = T, cluster_columns = T,
  rect_gp = gpar(type = "none"), column_dend_side = "top",
  cell_fun = function(j, i, x, y, w, h, fill) {
    if(as.numeric(x) >= 1 - as.numeric(y)) {
      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
    }},
  show_row_dend = FALSE,
  name = "J.I.", right_annotation = row_ann_nolegend,
  top_annotation = col_ann_nolegend, col = pal_seq, na_col = "white",
  show_row_names = FALSE, show_column_names = FALSE,
  heatmap_legend_param = list(direction = "horizontal"),
  show_heatmap_legend = FALSE)



##########################################################################################
##########################################################################################
##########################################################################################v


genes <- genes %>% 
  mutate(
    hclass = 
      factor(hclass, levels = c("aminoglycoside", "beta-lactamase", 
                                "macrolide", "tetracycline", "qnr")))


##########################################################################################
##########################################################################################
##########################################################################################


two_samples <- genes %>% select(centroid, City, Location, sample) %>% 
  group_by(centroid) %>% 
  mutate(N = n_distinct(sample)) %>% 
  filter(N > 1) %>% ungroup() %>% select(centroid) %>% distinct() %>% 
  pull() 

two_locations <- genes %>% select(centroid, City, Location, sample) %>% 
  group_by(centroid) %>% 
  mutate(N = n_distinct(Location)) %>% 
  filter(N > 1) %>% ungroup() %>% select(centroid) %>% distinct() %>% 
  pull() 

count_locations <- rbind(
  genes %>% 
    filter(centroid %in% two_locations) %>% 
    select(centroid, City, Location, sample) %>% 
    group_by(centroid) %>% 
    summarise(N = n_distinct(City)) %>% 
    ungroup() %>% 
    group_by(N) %>% 
    summarise(n = n()) %>% 
    mutate(gr = ifelse(N == 1, "Single city", "Multi city"),
           lev = "City"),
  
  genes %>% 
    filter(centroid %in% two_samples) %>% 
    select(centroid, City, Location, sample) %>% 
    group_by(centroid) %>% 
    summarise(N = n_distinct(Location)) %>% 
    ungroup() %>% 
    group_by(N) %>% 
    summarise(n = n()) %>% 
    mutate(gr = ifelse(N == 1, "Single location", "Multi location"),
           lev = "location"),
  
  
  genes %>% select(centroid, City, Location, sample) %>% 
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


count_locations <- count_locations %>%
  mutate(lev = factor(lev, levels = c("mag", "sample", "location", "City"))) %>% 
  ungroup() %>% 
  group_by(lev) %>%
  mutate(p = n / sum(n))

count_locations <- count_locations %>% 
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

count_locations$gr2 <- factor(count_locations$gr2, levels = c(1,2))

p_circle1 <- ggplot(count_locations, aes(x = "", y = p, fill = gr2)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y", start = pi/1.5) +
  geom_text(aes(label = comma(n)),
            position = position_stack(vjust = 0.5),
            size = (general_size) / ggplot2::.pt) +
  facet_wrap(lev ~ ., ncol = 1) +
  theme_void() +
  scale_fill_manual(values = c(pal_8[4:3],pal_8[4:3],pal_8[4:3])) +
  theme(
    strip.text = element_text(size = general_size, color = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.text.y  = element_blank(),
    legend.text = element_text(size = general_size),
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.background = element_blank())


genes_per_sample <- genes %>% 
  group_by(City, Location, sample) %>% 
  summarise(n = n_distinct(centroid)) 

wilcox.test(genes_per_sample$n[genes_per_sample$City == "Nanjing"],
            genes_per_sample$n[genes_per_sample$City != "Nanjing"], "two.sided", paired = F)


genes_per_city <- genes %>% group_by(City) %>% summarise(n = n_distinct(centroid)) %>% 
  mutate(City = paste(" ", City)) 

qantiles <- genes %>% 
  group_by(City, Location, sample) %>% 
  summarise(n = n_distinct(centroid)) %>%
  ungroup() %>% 
  group_by(City) %>% 
  summarise(q25 = quantile(n, 0.25), q75 = quantile(n, 0.75))


box2 <- genes %>% 
  group_by(City, Location, sample) %>% 
  summarise(n = n_distinct(centroid)) %>% 
  ggplot(aes(x = City, y = n, fill = City)) +
  #geom_rect(data = genes_per_city,
  #          aes(xmin = 0, xmax = as.numeric(factor(City)) + 0.4, 
  #              ymin = 0, ymax = n, fill = City), alpha = 0.3,
  #          inherit.aes = FALSE, color = "black", linewidth =  0.2,
  #          show.legend = F) +
  #geom_rect(data = qantiles,
  #          aes(xmin = as.numeric(factor(City)) + 2, 
  #              xmax = 5, 
  #              ymin = q25*10, ymax = q75*10, fill = City), alpha = 0.3,
  #          inherit.aes = FALSE, color = "black", linewidth =  0.2,
  #          show.legend = F) +
  #geom_col(data = genes_per_city, aes(x = City, y = n, fill = City),  
  #         width = 0.8, color = "black", show.legend = F) +
  geom_boxplot(show.legend = F) +
  scale_y_continuous(labels = label_comma()) +
  geom_jitter(color = "black", height = 0, width = 0.3, show.legend = F) + 
  xlab("") + 
  scale_fill_manual(values = pal_8[c(1,2,1,2)], labels = function(x) gsub(" ", "\n", x))  +
  theme_minimal() + 
  #scale_x_discrete(labels = function(x) {
  #  x <- gsub("-", "-\n", x)
  #  x <- gsub(" ", "\n", x)
  #  x}) + 
  #scale_y_continuous(name = "ARGs per city",
  #                   sec.axis = sec_axis(~ . / 10, name = "ARGs per sample")) +
  ylab("ARGs per sample") + 
  theme(
    legend.position = "bottom",
    axis.text.x.top = element_blank(),
    strip.text = element_text(size = general_size, face = "bold", color = "black"),
    legend.text = element_text(size = general_size),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size , color = "black", face = "bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = general_size, color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.background = element_blank())

box2

ggsave("figure6/6_b.svg", p_circle1 +
         theme(plot.margin = margin(0, 0, 0, 0), 
               panel.spacing = unit(5, "mm")) +
         ggtitle("B"), 
       width = 90, height = 100, unit = "mm")


# dev.off()
# CairoSVG("figure6/6_B.svg", 
#        width = 90/25.4, height = 100/25.4)
# print(  p_circle1 +
#           theme(plot.margin = margin(0, 0, 0, 0), 
#                 panel.spacing = unit(5, "mm")))
# dev.off()

ggsave("figure6/6_a.svg", box2 + 
         theme(axis.text.x = element_blank(),
               plot.margin = margin(0, 0, 0, 0), 
               panel.spacing = unit(5, "mm")) +
         ggtitle("A"), 
       width = 90, height = 100, unit = "mm")



dev.off()
svg("figure6/6_c.svg", width = 70/25.4, height = 100/25.4)
draw(draw(pheat11, heatmap_legend_side = "bottom", merge_legend = TRUE))
grid.text("C", x = 0.5, y = unit(1, "npc") + unit(2, "mm"), gp = gpar(fontsize=general_size + 2))
dev.off()




heatmap_grob <- grid.grabExpr(draw(pheat11, heatmap_legend_side = "bottom", merge_legend = TRUE))

heatmap_grob <- arrangeGrob(
  heatmap_grob,
  #top = textGrob("C", just = "left", gp = gpar(fontsize = general_size + 2, color = "black", fontface = "bold")),
  padding = unit(c(0, 0, 0, 0), "mm")
)



grid.arrange(box2 + 
               theme(axis.text.x = element_blank(),
                     plot.margin = margin(0, 0, 0, 0), 
                     panel.spacing = unit(5, "mm")) +
               ggtitle("A"),
             p_circle1 +
               theme(plot.margin = margin(0, 0, 0, 0), 
                     panel.spacing = unit(5, "mm")) +
               ggtitle("B"),
             heatmap_grob, 
             nrow = 1)



#+
#  theme_minimal()

##########################################################################################
##########################################################################################
##########################################################################################
