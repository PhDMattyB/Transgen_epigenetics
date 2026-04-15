##############################
## Mapk pathway gene locations
##
## Matt Brachmann (PhDMattyB)
##
## 15.04.2026
##
##############################


library(tidyverse)
library(data.table)
library(clusterProfiler)
library(org.Dr.eg.db)

run_go <- function(gene_vec, ontology) {
  
  enrichGO(
    gene          = gene_vec,
    OrgDb         = org.Dr.eg.db,
    keyType       = "SYMBOL",
    ont           = ontology,     # "BP", "MF", or "CC"
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
}

simplify_GO <- function(go_obj) {
  
  if (is.null(go_obj)) return(NULL)
  
  clusterProfiler::simplify(
    go_obj,
    cutoff      = 0.7,
    by          = "p.adjust",
    select_fun  = min,
    measure     = "Wang"
  )
}


gene_annotation = read_tsv('~/Parsons_Postdoc/Stickleback_Genomic/Stickleback_Annotation_features/stickleback_v5_ensembl_genes.gff3.gz', 
                          col_names = F, 
                          skip = 1) %>% 
  # filter(X3 %in% c('gene', 
  #                  'exon', 
  #                  'CDS')) %>% 
  group_by(X1) %>% 
  arrange(X4, 
          X5) %>% 
  ## arrange each gene by its start and end points on each chromosome
  mutate(mid = X4 + (X5-X4)/2) %>% 
  dplyr::select(X1, 
                X3:X5, 
                X9:mid) %>% 
  dplyr::rename(chromosome = X1, 
         feature = X3, 
         start = X4, 
         end = X5, 
         gene_id = X9, 
         position = mid) %>% 
  na.omit() 
# %>% 
#   filter(feature == 'gene')




methy_outliers = read_tsv('~/Parsons_Postdoc/Stickleback_Genomic/Methylation_outliers.txt')%>% 
  separate(col = TETWarm, 
           into = c('chromosome', 
                    'position'),
           sep = '-') %>% 
  group_by(chromosome) %>% 
  arrange(position)

methy_outliers$position = as.numeric(methy_outliers$position)

## create a 1Kb window around the snp identified
## settting the start and end range for the methylation data
## this number can be as big or small as you want depending
## on how liberal you want it

methy_out_window = methy_outliers %>%
  group_by(chromosome) %>% 
  mutate(start = position-100, 
         end = position+100)

library(data.table)

setDT(methy_out_window)
setDT(gene_annotation)

setkey(methy_out_window, 
       chromosome, 
       start, 
       end)

gene_overlap = foverlaps(gene_annotation, 
                         methy_out_window, 
                         type="any")


gene_overlap_tib = as_tibble(gene_overlap) %>% 
  na.omit() 

View(gene_overlap_tib)

gene_name_1 = gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_id, 
                F118, 
                F218, 
                TETWarm.F118, 
                TETWarm.F218, 
                TETWarm.F118.F218) %>% 
  separate(col = gene_id, 
           into = c('ensemble_id', 
                    'gene_name', 
                    'parent_code', 
                    'gene_name2'), 
           sep = ';') %>%
  separate(col = gene_name, 
           into = c('Garbage', 
                    'gene_name'), 
           sep = '=') %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_name) %>% 
  na.omit()

gene_name_2 = gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_id,
                F118, 
                F218, 
                TETWarm.F118, 
                TETWarm.F218, 
                TETWarm.F118.F218) %>% 
  separate(col = gene_id, 
           into = c('ensemble_id', 
                    'gene_name', 
                    'parent_code', 
                    'gene_name2'), 
           sep = ';') %>%
  separate(col = parent_code, 
           into = c('Garbage', 
                    'gene_name'), 
           sep = '=') %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_name) %>% 
  na.omit()


methy_genes = bind_rows(gene_name_1, 
                        gene_name_2) %>% 
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>% 
  filter(!grepl('ENSG', 
                gene_name))%>% 
  separate(gene_name, 
           c('gene_name', 
             'trash'), 
           sep = '-') %>% 
  dplyr::select(-trash) %>% 
  distinct(gene_name)


# methy_genes %>% 
#   write_tsv('Methylation_outlier_genes.txt')

methy_genes_go = methy_genes%>%
  summarise(
    genes = list(unique(gene_name)),
    .groups = "drop"
  ) %>%
  mutate(
    BP = map(genes, ~ tryCatch(run_go(.x, "BP"), error = function(e) NULL)),
    MF = map(genes, ~ tryCatch(run_go(.x, "MF"), error = function(e) NULL)),
    CC = map(genes, ~ tryCatch(run_go(.x, "CC"), error = function(e) NULL))
  )%>%
  mutate(
    BP = map(BP, simplify_GO),
    MF = map(MF, simplify_GO),
    CC = map(CC, simplify_GO)
  ) %>% 
  pivot_longer(
    cols = c(BP, MF, CC),
    names_to = "ontology",
    values_to = "go"
  ) %>%
  filter(!map_lgl(go, is.null)) %>%
  mutate(go_tbl = map(go, as_tibble)) %>%
  dplyr::select(ontology, go_tbl) %>%
  unnest(go_tbl)
# %>% 
#   dplyr::filter(ontology == 'BP')




# align with mapk genes --------------------------------------------------


mapk_pathway = read_csv('~/Parsons_Postdoc/Stickleback_Genomic/Stickleback_Annotation_features/mapk_pathway_genelist_fixed.csv')%>%
  mutate_all(., .funs = tolower) %>% 
  dplyr::rename(gene_name = Symbol) %>% 
  dplyr::select(gene_name)

methy_genes = methy_genes %>% 
  mutate_all(., .funs = tolower)

methy_genes %>% 
  filter(gene_name %in% mapk_pathway$gene_name)


inner_join(methy_genes, 
           mapk_pathway, 
           by = 'gene_name')

intersect(mapk_pathway, 
          methy_genes)


# Outlier top1% expression distribution -----------------------------------
# ecotype outliers --------------------------------------------------------


# eco_div_outliers = read_csv("~/Parsons_Postdoc/Methylation_data/GLM_results/Methylation_Eco_outlier_top5dist.csv") %>% 
#   mutate(chromosome = Chromosome, 
#          start = win_start, 
#          end = win_end)

eco_div_outliers = read_csv("~/Parsons_Postdoc/Methylation_data/GLM_results/Methylation_Eco_outlier_top1dist.csv") %>% 
  mutate(chromosome = Chromosome, 
         start = win_start, 
         end = win_end)


setDT(eco_div_outliers)
setDT(gene_annotation)

setkey(eco_div_outliers, 
       chromosome, 
       start, 
       end)

eco_gene_overlap = foverlaps(gene_annotation, 
                         eco_div_outliers, 
                         type="any")


eco_gene_overlap_tib = as_tibble(eco_gene_overlap) %>% 
  na.omit()

# View(gene_overlap_tib)


eco_gene_name_1 = eco_gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_id, 
                win_mid, 
                expression_mean, 
                expression_sd, 
                feature) %>% 
  separate(col = gene_id, 
           into = c('ensemble_id', 
                    'gene_name', 
                    'parent_code', 
                    'gene_name2'), 
           sep = ';') %>%
  separate(col = gene_name, 
           into = c('Garbage', 
                    'gene_name'), 
           sep = '=') %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_name) %>% 
  na.omit()

eco_gene_name_2 = eco_gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_id, 
                win_mid, 
                expression_mean, 
                expression_sd, 
                feature) %>% 
  separate(col = gene_id, 
           into = c('ensemble_id', 
                    'gene_name', 
                    'parent_code', 
                    'gene_name2'), 
           sep = ';') %>%
  separate(col = parent_code, 
           into = c('Garbage', 
                    'gene_name'), 
           sep = '=') %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_name) %>% 
  na.omit()


eco_methy_genes = bind_rows(eco_gene_name_1, 
                            eco_gene_name_2) %>% 
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>% 
  filter(!grepl('ENSG', 
                gene_name))%>% 
  separate(gene_name, 
           c('gene_name', 
             'trash'), 
           sep = '-') %>% 
  dplyr::select(-trash) %>% 
  distinct(gene_name)%>% 
  mutate_all(., .funs = tolower)


mapk_pathway = read_csv('~/Parsons_Postdoc/Stickleback_Genomic/Stickleback_Annotation_features/mapk_pathway_genelist_fixed.csv')%>%
  mutate_all(., .funs = tolower) %>% 
  dplyr::rename(gene_name = Symbol) %>% 
  dplyr::select(gene_name)


eco_methy_genes %>% 
  filter(gene_name %in% mapk_pathway$gene_name) %>% 
  write_tsv("~/Parsons_Postdoc/Methylation_data/GLM_results/Eco_top1_outlier_mapk_genes.txt")


eco_methy_genes %>% 
  write_tsv('~/Parsons_Postdoc/Methylation_data/GLM_results/Eco_top1_Outlier_genes.txt', 
            col_names = F)







# F1 eco outliers -------------------------------------------------------------
# f1_div_outliers = read_csv("~/Parsons_Postdoc/Methylation_data/GLM_results/Methylation_F1_eco_outlier_top5dist.csv") %>% 
#   mutate(chromosome = Chromosome, 
#          start = win_start, 
#          end = win_end)

f1_div_outliers = read_csv("~/Parsons_Postdoc/Methylation_data/GLM_results/Methylation_F1_eco_outlier_top1dist.csv") %>% 
  mutate(chromosome = Chromosome, 
         start = win_start, 
         end = win_end)


setDT(f1_div_outliers)
setDT(gene_annotation)

setkey(f1_div_outliers, 
       chromosome, 
       start, 
       end)

f1_gene_overlap = foverlaps(gene_annotation, 
                         f1_div_outliers, 
                         type="any")


f1_gene_overlap_tib = as_tibble(f1_gene_overlap) %>% 
  na.omit()

# View(gene_overlap_tib)


F1_eco_gene_name_1 = f1_gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_id, 
                win_mid, 
                expression_mean, 
                expression_sd, 
                feature) %>% 
  separate(col = gene_id, 
           into = c('ensemble_id', 
                    'gene_name', 
                    'parent_code', 
                    'gene_name2'), 
           sep = ';') %>%
  separate(col = gene_name, 
           into = c('Garbage', 
                    'gene_name'), 
           sep = '=') %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_name) %>% 
  na.omit()

F1_eco_gene_name_2 = f1_gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_id, 
                win_mid, 
                expression_mean, 
                expression_sd, 
                feature) %>% 
  separate(col = gene_id, 
           into = c('ensemble_id', 
                    'gene_name', 
                    'parent_code', 
                    'gene_name2'), 
           sep = ';') %>%
  separate(col = parent_code, 
           into = c('Garbage', 
                    'gene_name'), 
           sep = '=') %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_name) %>% 
  na.omit()


F1_eco_methy_genes = bind_rows(F1_eco_gene_name_1, 
                               F1_eco_gene_name_2) %>% 
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>% 
  filter(!grepl('ENSG', 
                gene_name))%>% 
  separate(gene_name, 
           c('gene_name', 
             'trash'), 
           sep = '-') %>% 
  dplyr::select(-trash) %>% 
  distinct(gene_name)%>% 
  mutate_all(., .funs = tolower)


mapk_pathway = read_csv('~/Parsons_Postdoc/Stickleback_Genomic/Stickleback_Annotation_features/mapk_pathway_genelist_fixed.csv')%>%
  mutate_all(., .funs = tolower) %>% 
  dplyr::rename(gene_name = Symbol) %>% 
  dplyr::select(gene_name)


F1_eco_methy_genes %>% 
  filter(gene_name %in% mapk_pathway$gene_name) %>% 
  write_tsv("~/Parsons_Postdoc/Methylation_data/GLM_results/F1_eco_top1_outlier_mapk_genes.txt")

## fgf13 related to ecotype*f1 effects


F1_eco_methy_genes %>% 
  write_tsv('~/Parsons_Postdoc/Methylation_data/GLM_results/F1_Eco_top1_Outlier_genes.txt', 
            col_names = F)






# F2 eco outliers -------------------------------------------------------------

# f2_div_outliers = read_csv("~/Parsons_Postdoc/Methylation_data/GLM_results/Methylation_F2_eco_outlier_top5dist.csv") %>% 
#   mutate(chromosome = Chromosome, 
#          start = win_start, 
#          end = win_end)
f2_div_outliers = read_csv("~/Parsons_Postdoc/Methylation_data/GLM_results/Methylation_F2_eco_outlier_top1dist.csv") %>% 
  mutate(chromosome = Chromosome, 
         start = win_start, 
         end = win_end)


setDT(f2_div_outliers)
setDT(gene_annotation)

setkey(f2_div_outliers, 
       chromosome, 
       start, 
       end)

f2_gene_overlap = foverlaps(gene_annotation, 
                         f2_div_outliers, 
                         type="any")


f2_gene_overlap_tib = as_tibble(f2_gene_overlap) %>% 
  na.omit()

# View(gene_overlap_tib)


F2_eco_gene_name_1 = f2_gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_id, 
                win_mid, 
                expression_mean, 
                expression_sd, 
                feature) %>% 
  separate(col = gene_id, 
           into = c('ensemble_id', 
                    'gene_name', 
                    'parent_code', 
                    'gene_name2'), 
           sep = ';') %>%
  separate(col = gene_name, 
           into = c('Garbage', 
                    'gene_name'), 
           sep = '=') %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_name) %>% 
  na.omit()

F2_eco_gene_name_2 = f2_gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_id, 
                win_mid, 
                expression_mean, 
                expression_sd, 
                feature) %>% 
  separate(col = gene_id, 
           into = c('ensemble_id', 
                    'gene_name', 
                    'parent_code', 
                    'gene_name2'), 
           sep = ';') %>%
  separate(col = parent_code, 
           into = c('Garbage', 
                    'gene_name'), 
           sep = '=') %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_name) %>% 
  na.omit()


F2_eco_methy_genes = bind_rows(F2_eco_gene_name_1, 
                               F2_eco_gene_name_2) %>% 
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>% 
  filter(!grepl('ENSG', 
                gene_name))%>% 
  separate(gene_name, 
           c('gene_name', 
             'trash'), 
           sep = '-') %>% 
  dplyr::select(-trash) %>% 
  distinct(gene_name)%>% 
  mutate_all(., .funs = tolower)


mapk_pathway = read_csv('~/Parsons_Postdoc/Stickleback_Genomic/Stickleback_Annotation_features/mapk_pathway_genelist_fixed.csv')%>%
  mutate_all(., .funs = tolower) %>% 
  dplyr::rename(gene_name = Symbol) %>% 
  dplyr::select(gene_name)


F2_eco_methy_genes %>% 
  filter(gene_name %in% mapk_pathway$gene_name) %>% 
  write_tsv("~/Parsons_Postdoc/Methylation_data/GLM_results/F2_eco_top1_outlier_mapk_genes.txt")


F2_eco_methy_genes %>% 
  write_tsv('~/Parsons_Postdoc/Methylation_data/GLM_results/F2_Eco_top1_Outlier_genes.txt', 
            col_names = F)


# f1 f2 eco outliers ------------------------------------------------------
# F1F2Eco_div_outliers = read_csv("~/Parsons_Postdoc/Methylation_data/GLM_results/Methylation_F1_F2_eco_outlier_top5dist.csv") %>% 
#   mutate(chromosome = Chromosome, 
#          start = win_start, 
#          end = win_end)
F1F2Eco_div_outliers = read_csv("~/Parsons_Postdoc/Methylation_data/GLM_results/Methylation_F1_F2_eco_outlier_top1dist.csv") %>% 
  mutate(chromosome = Chromosome, 
         start = win_start, 
         end = win_end)


setDT(F1F2Eco_div_outliers)
setDT(gene_annotation)

setkey(F1F2Eco_div_outliers, 
       chromosome, 
       start, 
       end)

f1_f2_gene_overlap = foverlaps(gene_annotation, 
                         F1F2Eco_div_outliers, 
                         type="any")


f1_f2_gene_overlap_tib = as_tibble(f1_f2_gene_overlap) %>% na.omit()

# View(gene_overlap_tib)


F1F2_eco_gene_name_1 = f1_f2_gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_id, 
                win_mid, 
                expression_mean, 
                expression_sd, 
                feature) %>% 
  separate(col = gene_id, 
           into = c('ensemble_id', 
                    'gene_name', 
                    'parent_code', 
                    'gene_name2'), 
           sep = ';') %>%
  separate(col = gene_name, 
           into = c('Garbage', 
                    'gene_name'), 
           sep = '=') %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_name) %>% 
  na.omit()

F1F2_eco_gene_name_2 = f1_f2_gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_id, 
                win_mid, 
                expression_mean, 
                expression_sd, 
                feature) %>% 
  separate(col = gene_id, 
           into = c('ensemble_id', 
                    'gene_name', 
                    'parent_code', 
                    'gene_name2'), 
           sep = ';') %>%
  separate(col = parent_code, 
           into = c('Garbage', 
                    'gene_name'), 
           sep = '=') %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_name) %>% 
  na.omit()


F1F2_eco_methy_genes = bind_rows(F1F2_eco_gene_name_1, 
                                 F1F2_eco_gene_name_2) %>% 
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>% 
  filter(!grepl('ENSG', 
                gene_name))%>% 
  separate(gene_name, 
           c('gene_name', 
             'trash'), 
           sep = '-') %>% 
  dplyr::select(-trash) %>% 
  distinct(gene_name)%>% 
  mutate_all(., .funs = tolower)


mapk_pathway = read_csv('~/Parsons_Postdoc/Stickleback_Genomic/Stickleback_Annotation_features/mapk_pathway_genelist_fixed.csv')%>%
  mutate_all(., .funs = tolower) %>% 
  dplyr::rename(gene_name = Symbol) %>% 
  dplyr::select(gene_name)


F1F2_eco_methy_genes %>% 
  filter(gene_name %in% mapk_pathway$gene_name)

F1F2_eco_methy_genes %>% 
  write_tsv('~/Parsons_Postdoc/Methylation_data/GLM_results/F1F2_Eco_top5_Outlier_genes.txt', 
            col_names = F)



# quantify overlap between model res top1 abs outliers --------------------
F1_eco_methy_genes
F2_eco_methy_genes
eco_methy_genes

intersect(eco_methy_genes, 
          F1_eco_methy_genes)

intersect(eco_methy_genes, 
          F2_eco_methy_genes)

intersect(F1_eco_methy_genes, 
          F2_eco_methy_genes)

inner_join(eco_methy_genes, 
           F1_eco_methy_genes) %>% 
  inner_join(., 
             F2_eco_methy_genes)





# enrichR -----------------------------------------------------------------

F2_eco_methy_genes_go = F2_eco_methy_genes%>%
  summarise(
    genes = list(unique(gene_name)),
    .groups = "drop"
  ) %>%
  mutate(
    BP = map(genes, ~ tryCatch(run_go(.x, "BP"), error = function(e) NULL)),
    MF = map(genes, ~ tryCatch(run_go(.x, "MF"), error = function(e) NULL)),
    CC = map(genes, ~ tryCatch(run_go(.x, "CC"), error = function(e) NULL))
  )%>%
  # mutate(
  #   BP = map(BP, simplify_GO),
  #   MF = map(MF, simplify_GO),
  #   CC = map(CC, simplify_GO)
  # ) %>% 
  pivot_longer(
    cols = c(BP, MF, CC),
    names_to = "ontology",
    values_to = "go"
  ) %>%
  filter(!map_lgl(go, is.null)) %>%
  mutate(go_tbl = map(go, as_tibble)) %>%
  dplyr::select(ontology, go_tbl) %>%
  unnest(go_tbl)


F1_eco_methy_genes_go = F1_eco_methy_genes%>%
  summarise(
    genes = list(unique(gene_name)),
    .groups = "drop"
  ) %>%
  mutate(
    BP = map(genes, ~ tryCatch(run_go(.x, "BP"), error = function(e) NULL)),
    MF = map(genes, ~ tryCatch(run_go(.x, "MF"), error = function(e) NULL)),
    CC = map(genes, ~ tryCatch(run_go(.x, "CC"), error = function(e) NULL))
  )%>%
  # mutate(
  #   BP = map(BP, simplify_GO),
  #   MF = map(MF, simplify_GO),
  #   CC = map(CC, simplify_GO)
  # ) %>% 
  pivot_longer(
    cols = c(BP, MF, CC),
    names_to = "ontology",
    values_to = "go"
  ) %>%
  filter(!map_lgl(go, is.null)) %>%
  mutate(go_tbl = map(go, as_tibble)) %>%
  dplyr::select(ontology, go_tbl) %>%
  unnest(go_tbl)


eco_methy_genes_go = eco_methy_genes%>%
  summarise(
    genes = list(unique(gene_name)),
    .groups = "drop"
  ) %>%
  mutate(
    BP = map(genes, ~ tryCatch(run_go(.x, "BP"), error = function(e) NULL)),
    MF = map(genes, ~ tryCatch(run_go(.x, "MF"), error = function(e) NULL)),
    CC = map(genes, ~ tryCatch(run_go(.x, "CC"), error = function(e) NULL))
  )%>%
  # mutate(
  #   BP = map(BP, simplify_GO),
  #   MF = map(MF, simplify_GO),
  #   CC = map(CC, simplify_GO)
  # ) %>% 
  pivot_longer(
    cols = c(BP, MF, CC),
    names_to = "ontology",
    values_to = "go"
  ) %>%
  filter(!map_lgl(go, is.null)) %>%
  mutate(go_tbl = map(go, as_tibble)) %>%
  dplyr::select(ontology, go_tbl) %>%
  unnest(go_tbl)
