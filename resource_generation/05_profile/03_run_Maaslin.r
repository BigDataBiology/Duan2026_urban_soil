#Maaslin3
data <- read.delim("genus_abundance_filter.tsv",row.names = 1, header = TRUE, sep = "\t")
property <- read.delim("property.csv", row.names = 1,header = TRUE, sep = ",")
property$City <-factor(property$City, levels = c('ShangHai', 'NanJing'))
fit_out <- maaslin3(input_data = data,
                    input_metadata = property,
                    output = 'genus_3_abundance_output',
                    formula = '~ Electric.conductivity + All.carbon +
                    Organic.carbon + N + P + K + City',
                    normalization = 'None',
                    transform = 'LOG',
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    min_prevalence=0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 100,
                    cores = 1,
                    save_models = TRUE)
#Maaslin2
fit_data = Maaslin2(input_data     = data, 
                    input_metadata = property, 
                    output         = "genus_2_abundance_noTSS_output", 
                    max_significance = 0.1,
                    normalization = "None",
                    fixed_effects  = c("Electric.conductivity","All.carbon","Organic.carbon","N","P","K","City"),
                    reference      = c("City,ShangHai"))