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
  c("aac2p",	"aminoglycoside acetyltransferase",	"aac2p", "aminoglycoside"),
  c("aac3_1",	"aminoglycoside acetyltransferase",	"aac3", "aminoglycoside"), c("aac3_2",	"aminoglycoside acetyltransferase",	"aac3", "aminoglycoside"),
  c("aac6p_1",	"aminoglycoside acetyltransferase",	"aac6p'", "aminoglycoside"), c("aac6p_2",	"aminoglycoside acetyltransferase",	"aac6p'", "aminoglycoside"),
  c("aac6p_3",	"aminoglycoside acetyltransferase",	"aac6p'", "aminoglycoside"), c("aph2b",	"aminoglycoside phosphotransferase", "aph2''", "aminoglycoside"),
  c("aph3p",	"aminoglycoside phosphotransferase", "aph3'", "aminoglycoside"), c("aph6p",	"aminoglycoside phosphotransferase", "aph6'", "aminoglycoside"),
  c("class_a",	"beta-lactamase A",	"beta-lactamase A", "beta-lactamase"), c("class_b1_b2",	"beta-lactamase B",	"beta-lactamase B1B2", "beta-lactamase"),
  c("class_b3",	"beta-lactamase B",	"beta-lactamase B3", "beta-lactamase"), c("class_c",	"beta-lactamase C",	"beta-lactamase C", "beta-lactamase"),
  c("class_d1",	"beta-lactamase D",	"beta-lactamase D", "beta-lactamase"), c("class_d2",	"beta-lactamase D",	"beta-lactamase D", "beta-lactamase"),
  c("erm_1",	"macrolide erm", "erm", "macrolide"), c("erm_2",	"macrolide erm", "erm",	'macrolide'),
  c("mph",	"macrolide mph", "mph",	"macrolide"), c("qnr",	"qrn", "qnr", "qnr"), c("tet_efflux",	"tetracycline efflux pump", "tet efflux", "tetracycline"),
  c("tet_enzyme",	"tetracycline inactivating enzymes", "tet enzyme", "tetracycline"), c("tet_rpg",	"tetracycline ribosomal protection genes", "tet rpg", "tetracycline")))

metadata <- readxl::read_excel("../resource_generation/12_ARGs_fARGene/supp1_sample_metadata.xlsx")
metadata <- metadata %>% mutate(Sample_id = gsub("ample", "", Sample_id)) 

genes <- read.table("../resource_generation/12_ARGs_fARGene/Assembly/predicted_gene_names.txt", quote="\"", comment.char="")
clusters <- read.table("../resource_generation/12_ARGs_fARGene/Assembly/clusters.uc", quote="\"", comment.char="")

genes <- genes %>% mutate(centroid = clusters$V10[match(V1, clusters$V9)])
genes <- genes %>% mutate(centroid = ifelse(centroid == "*", V1, centroid))



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



rm( split1, split2, split3, split4)



######################################
######################################

## add metadata 


genes <- genes %>% mutate(description = gene_classes$X2[match(c_hmm, gene_classes$X1)],
                          class = gene_classes$X3[match(c_hmm, gene_classes$X1)],
                          hclass = gene_classes$X4[match(c_hmm, gene_classes$X1)],
                          City = metadata$City[match(sample, metadata$Sample_id)],
                          Location = metadata$Location[match(sample, metadata$Sample_id)],
                          Date = metadata$Date[match(sample, metadata$Sample_id)],
                          Longitude = metadata$Date[match(sample, metadata$Longitude)],
                          Latitude = metadata$Date[match(sample, metadata$Latitude)])

genes <- genes %>% mutate(read = paste(contig, sample, sep = "_")) 

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



amrfinder_genes <- read.delim("../resource_generation/12_ARGs_fARGene/amrfinder_on_fargene_centroids95_prodigal.tsv")
amrfinder_genes <- amrfinder_genes %>% 
  mutate(centroid =  Protein.id) %>% 
  mutate(centroid = sub("_1$", "", centroid)) %>% 
  mutate(centroid = sub("_2$", "", centroid)) %>% 
  mutate(centroid = sub("_3$", "", centroid)) %>% 
  mutate(centroid = sub("_4$", "", centroid)) %>% 
  mutate(centroid = sub("_5$", "", centroid)) %>% 
  mutate(centroid = sub("_6$", "", centroid)) 

amrfinder_genes <- amrfinder_genes %>% 
  mutate(centroid = ifelse(centroid %in% genes$centroid, centroid,
                           ifelse(Protein.id %in% genes$centroid, Protein.id, NA)))

genes <- genes %>% mutate(amrfinder = ifelse(centroid %in% amrfinder_genes$centroid, 1, 0))
genes <- genes %>% mutate(amrfinder.class = amrfinder_genes$Class[match(centroid, amrfinder_genes$centroid)])
genes <- genes %>% mutate(amrfindersubclass = amrfinder_genes$Subclass[match(centroid, amrfinder_genes$centroid)])
genes <- genes %>% mutate(amrfinder.hmm = amrfinder_genes$HMM.description[match(centroid, amrfinder_genes$centroid)])

eggnog <- read.delim("~/Documents/GitHub/urban_soil/resource_generation/12_ARGs_fARGene/centroids.emapper.annotations", comment.char="#")
eggnog <- eggnog %>% 
  mutate(centroid =  query) %>% 
  mutate(centroid = sub("_1$", "", centroid)) %>% 
  mutate(centroid = sub("_2$", "", centroid)) %>% 
  mutate(centroid = sub("_3$", "", centroid)) %>% 
  mutate(centroid = sub("_4$", "", centroid)) %>% 
  mutate(centroid = sub("_5$", "", centroid)) %>% 
  mutate(centroid = sub("_6$", "", centroid)) 

genes <- genes %>% mutate(eggnog_description = eggnog$Description[match(centroid, eggnog$centroid)])
genes <- genes %>% mutate(eggnog_preferred_name = eggnog$Preferred_name[match(centroid, eggnog$centroid)])
genes <- genes %>% mutate(eggnog_pfam = eggnog$PFAMs[match(centroid, eggnog$centroid)])

genes <- genes %>% group_by(centroid) %>% slice_head(n = 1)

table(genes$hmm, genes$amrfinder.class)
table(genes$hmm, genes$eggnog_pfam)
      
genes <- genes %>% mutate(amrfinder.class = ifelse(is.na(amrfinder.class), "Not found", amrfinder.class))
genes <- genes %>% mutate(amrfindersubclass = ifelse(is.na(amrfindersubclass), "Not found", amrfindersubclass))

col_ann = HeatmapAnnotation(
  City = metadata$City[match(colnames(mat), metadata$Sample_id) ],
  col = list(City = setNames(pal_8[1:2], cities)))

row_ann = rowAnnotation(
  City     = metadata$City[match(rownames(mat), metadata$Sample_id)],
  col = list(City = setNames(pal_8[1:2], cities)),
  show_legend = FALSE)



genes$eggnog_pfam[is.na(genes$eggnog_pfam)] <- "Not found"
mat_pfam <- table(genes$eggnog_pfam, genes$class)
mat_pfam <- t(t(mat_pfam)/colSums(mat_pfam))
mat_pfam[mat_pfam == 0] <- NA
mat_pfam <- mat_pfam[c(rownames(mat_pfam)[!rownames(mat_pfam) %in% "Not found"],"Not found"),]

# heatmap Jaccard's index
pheat_pfam <- ComplexHeatmap::Heatmap(
  mat_pfam, cluster_rows = F, cluster_columns = F,
  col = pal_seq, na_col = "white",
  name = "% of gene class",
)
pheat_pfam

genes$eggnog_description[is.na(genes$eggnog_description)] <- "Not found"
mat_egg_desc <- table(genes$eggnog_description, genes$class)
mat_egg_desc <- t(t(mat_egg_desc)/colSums(mat_egg_desc))
mat_egg_desc[mat_egg_desc == 0] <- NA
mat_egg_desc <- mat_egg_desc[c(rownames(mat_egg_desc)[!rownames(mat_egg_desc) %in% "Not found"],"Not found"),]

# heatmap Jaccard's index
pheat_egg_desc <- ComplexHeatmap::Heatmap(
  mat_egg_desc, cluster_rows = F, cluster_columns = F,
  col = pal_seq, na_col = "white"
)
pheat_egg_desc


genes$amrfinder.class[is.na(genes$amrfinder.class)] <- "Not found"
mat_amrfinderclass <- table(genes$amrfinder.class, genes$class)
mat_amrfinderclass <- t(t(mat_amrfinderclass)/colSums(mat_amrfinderclass))
mat_amrfinderclass[mat_amrfinderclass == 0] <- NA
mat_amrfinderclass <- mat_amrfinderclass[c(rownames(mat_amrfinderclass)[!rownames(mat_amrfinderclass) %in% "Not found"],"Not found"),]

# heatmap Jaccard's index
pheat_amrfinderclass <- ComplexHeatmap::Heatmap(
  mat_amrfinderclass, cluster_rows = F, cluster_columns = F,
  col = pal_seq, na_col = "white"
)

pheat_amrfinderclass


genes %>% group_by(eggnog_pfam) %>% 
  summarise(n = n()) %>% 
  mutate(p = n / sum(n)) %>% 
  filter(eggnog_pfam %in% c("-", "Not found")) %>% 
  mutate(P = sum(p))

data.frame(genes %>% group_by(class, eggnog_pfam) %>% summarise(n = n()) %>% mutate(p = n / sum(n)))
weird_eggnog_pfam <- c("BPD_transp_2",
                       "FtsX,MacB_PCD",
                       "MTHFR",
                       "Ribosom_S12_S23",
                       "ZapA",
                       "Peptidase_M56,Transpeptidase",
                       "Transpeptidase",
                       "RrnaAD",
                       "GTP-bdg_M,GTP-bdg_N,MMR_HSR1",
                       "Pentapeptide,Pentapeptide_4",
                       #"MFS_1,Sugar_tr"
                       "FAD_binding_3",
                       "EFG_C,EFG_II,EFG_IV,GTP_EFTU",
                       "EFG_C,EFG_II,EFG_IV,GTP_EFTU,GTP_EFTU_D2", "-")

genes %>% filter(eggnog_pfam %in% weird_eggnog_pfam) %>% group_by(class, eggnog_description) %>% summarise(n = n()) %>% mutate(p = n / sum(n))



eggnog_description_ok <- c(
  "beta-lactamase",
  "Beta-lactamase",
  "COG2602 Beta-lactamase class D",
  "PFAM Penicillin binding protein transpeptidase domain",
  "PFAM penicillin-binding protein transpeptidase",
  "Penicillin binding protein transpeptidase domain",
  "penicillin binding",
  "fluoroquinolone resistance protein",
  "Tetracycline resistance protein"
)


weird_description <- c(
"Activator of cell division through the inhibition of FtsZ GTPase activity, therefore promoting FtsZ assembly into bundles of protofilaments necessary for the formation of the division Z ring. It is recruited early at mid-cell but it is not essential for cell division",
"cell division through the inhibition of FtsZ GTPase activity, therefore promoting FtsZ assembly into bundles of protofilaments necessary for the formation of the division Z ring. It is recruited early at mid-cell but it is not essential for cell division",
"Belongs to the binding-protein-dependent transport system permease family",
"Interacts with and stabilizes bases of the 16S rRNA that are involved in tRNA selection in the A site and with the mRNA backbone. Located at the interface of the 30S and 50S subunits, it traverses the body of the 30S subunit contacting proteins on the other side and probably holding the rRNA structure together. The combined cluster of proteins S8, S12 and S17 appears to hold together the shoulder and platform of the 30S subunit",
"MacB-like periplasmic core domain",
"Methylenetetrahydrofolate reductase",
"Belongs to the class I-like SAM-binding methyltransferase superfamily. rRNA adenine N(6)-methyltransferase family",
"Ribosomal RNA adenine dimethylase",
"rRNA (adenine-N6,N6-)-dimethyltransferase activity",
"GTPase that associates with the 50S ribosomal subunit and may have a role during protein synthesis or ribosome biogenesis",
"COG0654 2-polyprenyl-6-methoxyphenol hydroxylase and related FAD-dependent oxidoreductases",
"FAD binding domain",
"Catalyzes the GTP-dependent ribosomal translocation step during translation elongation. During this step, the ribosome changes from the pre-translocational (PRE) to the post- translocational (POST) state as the newly formed A-site-bound peptidyl-tRNA and P-site-bound deacylated tRNA move to the P and E sites, respectively. Catalyzes the coordinated movement of the two tRNA molecules, the mRNA and conformational changes in the ribosome",
"Elongation factor G, domain IV",
"GTP-binding protein",
"elongation factor G",
"-")



genes %>% filter(eggnog_pfam %in% weird_eggnog_pfam & eggnog_description %in% weird_description) %>% group_by(class, eggnog_preferred_name) %>% summarise(n = n()) %>% mutate(p = n / sum(n))

weird_names <-c(
  "metF",
  "rbsC",
  "rpsL", #unspecific/mutations confer AMR
  "-",
  "ksgA") #unspecific/mutations confer AMR
  
unspecific_names <- c(
  "rpsL", 
  "ksgA") 



genes <- genes %>% 
  mutate(eggnog_no_amr_info = 
           ifelse(eggnog_pfam %in% weird_eggnog_pfam & 
                    eggnog_description %in% weird_description & 
                    eggnog_preferred_name %in% weird_names, "no AMR info","yes AMR info"))
  
data.frame(genes %>% group_by(class, eggnog_no_amr_info) %>% summarise(n = n()) %>% mutate(p = n/sum(n)))
data.frame(genes %>% group_by( eggnog_no_amr_info) %>% summarise(n = n()) %>% mutate(p = n/sum(n)))



genes %>% filter(eggnog_no_amr_info %in% "no AMR info") %>% group_by(class, eggnog_description, eggnog_preferred_name, eggnog_pfam) %>% summarise(n = n()) %>% mutate(p = n / sum(n))

genes %>% mutate(in_pfam = ifelse(eggnog_pfam %in% weird_eggnog_pfam, "no amr", "yes amr")) %>% group_by(in_pfam) %>% summarise(n = n()) %>% mutate(p = n / sum(n))
genes %>% filter(eggnog_pfam %in% weird_eggnog_pfam) %>%  mutate(in_description = ifelse(eggnog_description %in% weird_description, "no amr", "yes amr")) %>% group_by(in_description) %>% summarise(n = n()) %>% mutate(p = n / sum(n))
genes %>% filter(eggnog_pfam %in% weird_eggnog_pfam, eggnog_description %in% weird_description) %>%  mutate(in_genename = ifelse(eggnog_preferred_name %in% weird_names, "no amr", "yes amr")) %>% group_by(in_genename) %>% summarise(n = n()) %>% mutate(p = n / sum(n))



amr_eggnog_summary <- rbind(genes %>% mutate(in_pfam = ifelse(eggnog_pfam %in% weird_eggnog_pfam, "no amr", "yes amr")) %>% group_by(in_pfam) %>% summarise(n = n()) %>% mutate(p = n / sum(n)) %>% filter(in_pfam %in% "yes amr") %>% mutate(level = "pfam") %>% rename(amr = in_pfam),
genes %>% filter(eggnog_pfam %in% weird_eggnog_pfam) %>%  mutate(in_description = ifelse(eggnog_description %in% weird_description, "no amr", "yes amr")) %>% group_by(in_description) %>% summarise(n = n()) %>% mutate(p = n / sum(n)) %>% filter(in_description %in% "yes amr") %>% mutate(level = "description") %>% rename(amr = in_description),
genes %>% filter(eggnog_pfam %in% weird_eggnog_pfam, eggnog_description %in% weird_description) %>%  mutate(in_genename = ifelse(eggnog_preferred_name %in% weird_names, "no amr", "yes amr")) %>% group_by(in_genename) %>% summarise(n = n()) %>% mutate(p = n / sum(n))  %>% mutate(level = "gene name") %>% rename(amr = in_genename))

amr_eggnog_summary %>% arrange(desc(n)) %>% mutate(p = n / sum(n))


# number of genes | percentage of the total | description
# 30851 | 96.4% found at the pfam level
# 760 | 2.38%  found at the description level
# 173 | 0.541% found at the preferred gene name level
# 214 | 0.669% no info on amr (see note below)
# 
# of the 214 with no amr info:
# 71 FAD binding - tetracycline inactivation enzymes (that is how some of them work https://pmc.ncbi.nlm.nih.gov/articles/PMC7229144/)
# 118 methyltransferases - that is the function of macrolide methyltransferases
# 7 out of 17766 beta lactamase b3 no info on resistance
# 4 tet rpg no info on resistance
# 14 mph no info on resistance





