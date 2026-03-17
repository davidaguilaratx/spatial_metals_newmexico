setwd("C:/Users/david/Documents/Research/PhD/Spatial_Omics/Xenium")

library(tidyverse)

# load in xenuium panel and 100 custom add-on genes
xenium_panel_genes = read.csv('custom_panel_building/datasets/xenium_panel_genes.csv')
custom_genes = read.csv('Final_custom_gene_panel/gene_list.csv')

# format for merging
colnames(xenium_panel_genes)[1] = 'Gene'
custom_genes$Annotation = NA
custom_genes = custom_genes %>% select(Gene, Ensembl.ID, Annotation)

# merge xenium panel with custom genes
panel = rbind(xenium_panel_genes, custom_genes)

# Load DE genes from CELLxGENE between breast tissues of African American and European descent.
# Positive effect sizes pertain to African American (AA) women, while negative effect
# sizes pertain to women of European descent (ED).
# Data for AA come from Kumar et al. Nature 2023 and 
# Bhat-Nakshatri et al. Cell Reports Medicine. 2021.
# Data for ED come from the same source, but also Gray et al. Developmental Cell. 2022,
# and Reed et al. Nat Genet. 2024.
de_genes = read.csv('Final_custom_gene_panel/differential_expression_results_CELLxGENE_AAvsED.csv',
                    skip=2,header=T)

# merge DE genes from CELLxGENE data to get estimated effect size for
# for each gene in full customized xenium panel.
panel_de_genes = merge(panel, de_genes, by='Gene')
write.csv(panel_de_genes, 'Final_custom_gene_panel/panel_de_genes_and_effect_sizes_AAvsED.csv')

summary(panel_de_genes$Effect.Size)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.30200 -0.04100  0.00000  0.01286  0.04800  0.41100 

# positive, AA only
summary(panel_de_genes$Effect.Size[panel_de_genes$Effect.Size > 0])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00200 0.02000 0.05000 0.08308 0.11600 0.41100
length(panel_de_genes$Effect.Size[panel_de_genes$Effect.Size > 0])
# [1] 173


# Negative, ED only
summary(panel_de_genes$Effect.Size[panel_de_genes$Effect.Size < 0])
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.30200 -0.07650 -0.04100 -0.05619 -0.01500 -0.00100 
length(panel_de_genes$Effect.Size[panel_de_genes$Effect.Size < 0])
# [1] 175


