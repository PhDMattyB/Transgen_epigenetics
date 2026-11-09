## RDA results for TGP and WGP

library(tidyverse)
library(data.table)

setwd('~/Parsons_Postdoc/Methylation_data/GLM_results')

# gene annotation ---------------------------------------------------------

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
  dplyr::rename(Chromosome = X1, 
                feature = X3, 
                win_start = X4, 
                win_end = X5, 
                gene_id = X9, 
                win_mid = mid) %>% 
  na.omit() %>% 
  mutate(win_start_loc = win_start, 
         win_end_loc = win_end)


# Head WGP ----------------------------------------------------------------


head_wgp = read_csv('~/Methylation_data/WGP_HEAD_clean_RDA_CAND_corr.csv') %>% 
  separate(col = loc, 
           into = c('Chromosome', 
                    'POS'), 
           sep = '-') %>% 
  mutate(wgp_HEAD_clean_great = abs(WGP_HEAD_clean) > abs(eco_num) & abs(WGP_HEAD_clean) > abs(csize)) %>% 
  filter(wgp_HEAD_clean_great == 'TRUE') %>% 
  # filter(Association == 'WGP_HEAD_clean') %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  mutate(POS = as.numeric(POS))


F2_eco_res = read_tsv('Methylation_F2_eco_effect_win1000_step1000.txt') %>% 
  mutate(Effect = 'F2', 
         abs_effect = abs(expression_mean)) %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  filter(expression_n > 3) %>%
  group_by(Chromosome) %>%
  mutate(threshold = quantile(abs_effect,
                              0.99,
                              na.rm = TRUE)) %>%
  filter(abs_effect >= threshold) %>% 
  mutate(direction = case_when(
    expression_mean > threshold ~ "Hyper",
    expression_mean < -threshold ~ "Hypo",
    TRUE ~ "None")) %>%
  # Order by genomic position
  arrange(Chromosome, win_mid)


head_RDA_WGP_Ecotype = inner_join(head_wgp,
           F2_eco_res,
           by = join_by(Chromosome,
                        POS <= win_end,
                        POS >= win_start),
           suffix = c("_RDA", "_Meth"),
           relationship = "many-to-many") %>% 
  mutate(win_start_loc = POS-1, 
         win_end_loc = POS+1)



setDT(head_RDA_WGP_Ecotype)
setDT(gene_annotation)

setkey(head_RDA_WGP_Ecotype, 
       Chromosome, 
       win_start_loc, 
       win_end_loc)


head_RDA_WGP_Ecotype_gene_overlap = foverlaps(gene_annotation,
                                head_RDA_WGP_Ecotype,
                                # by.x = start,
                                # by.y = end,
                                type="any")


head_RDA_WGP_Ecotype_gene_overlap = as_tibble(head_RDA_WGP_Ecotype_gene_overlap) %>% 
  na.omit() %>%
  # filter(Chromosome != 'chrUn') %>% 
  arrange(Chromosome, 
          i.win_mid)

head_RDA_WGP_Ecotype_gene_overlap %>% 
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
  dplyr::select(gene_name) %>% 
  unique() %>% 
  write_tsv('RDA_HEAD_WGP_ECO_GENES.txt')


# Head tgp ----------------------------------------------------------------



head_tgp = read_csv('~/Methylation_data/TGP_HEAD_clean_RDA_CAND_corr.csv') %>% 
  separate(col = loc, 
           into = c('Chromosome', 
                    'POS'), 
           sep = '-')%>% 
  mutate(tgp_HEAD_clean_great = abs(TGP_HEAD_clean) > abs(eco_num) & abs(TGP_HEAD_clean) > abs(csize)) %>% 
  filter(tgp_HEAD_clean_great == 'TRUE') %>% 
  # filter(Association == 'TGP_HEAD_clean') %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  mutate(POS = as.numeric(POS))


F1_eco_res = read_tsv('Methylation_F1_eco_effect_win1000_step1000.txt') %>% 
  mutate(Effect = 'F1', 
         abs_effect = abs(expression_mean)) %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  filter(expression_n > 3) %>%
  group_by(Chromosome) %>%
  mutate(threshold = quantile(abs_effect,
                              0.99,
                              na.rm = TRUE)) %>%
  filter(abs_effect >= threshold) %>% 
  mutate(direction = case_when(
    expression_mean > threshold ~ "Hyper",
    expression_mean < -threshold ~ "Hypo",
    TRUE ~ "None")) %>%
  # Order by genomic position
  arrange(Chromosome, win_mid)


head_RDA_TGP_Ecotype = inner_join(head_tgp,
                                  F1_eco_res,
                                  by = join_by(Chromosome,
                                               POS <= win_end,
                                               POS >= win_start),
                                  suffix = c("_RDA", "_Meth"),
                                  relationship = "many-to-many") 




# Head RDA F1 F2 ecotype interactions -------------------------------------



F1_F2_eco_res = read_tsv('Methylation_F1_F2_eco_effect_win1000_step1000.txt') %>% 
  mutate(Effect = 'F1_F2', 
         abs_effect = abs(expression_mean)) %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  filter(expression_n > 3) %>%
  group_by(Chromosome) %>%
  mutate(threshold = quantile(abs_effect,
                              0.99,
                              na.rm = TRUE)) %>%
  filter(abs_effect >= threshold) %>% 
  mutate(direction = case_when(
    expression_mean > threshold ~ "Hyper",
    expression_mean < -threshold ~ "Hypo",
    TRUE ~ "None")) %>%
  # Order by genomic position
  arrange(Chromosome, win_mid)




head_RDA_WGP_f1_f2 = inner_join(head_wgp,
           F1_F2_eco_res,
           by = join_by(Chromosome,
                        POS <= win_end,
                        POS >= win_start),
           suffix = c("_RDA", "_Meth"),
           relationship = "many-to-many") %>% 
  mutate(win_start_loc = POS-1, 
         win_end_loc = POS+1)


head_RDA_WGP_Ecotype %>% 
  group_by(Chromosome, 
           POS) %>% 
  intersect(., head_RDA_WGP_f1_f2)


setDT(head_RDA_WGP_f1_f2)
setDT(gene_annotation)

setkey(head_RDA_WGP_f1_f2, 
       Chromosome, 
       win_start_loc, 
       win_end_loc)


head_RDA_WGP_f1_f2_gene_overlap = foverlaps(gene_annotation,
                                head_RDA_WGP_f1_f2,
                                # by.x = start,
                                # by.y = end,
                                type="any")


head_RDA_WGP_f1_f2_gene_overlap = as_tibble(head_RDA_WGP_f1_f2_gene_overlap) %>% 
  na.omit() %>%
  # filter(Chromosome != 'chrUn') %>% 
  arrange(Chromosome, 
          i.win_mid)

head_RDA_WGP_f1_f2_gene_overlap %>% 
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
  dplyr::select(gene_name) %>% 
  unique() %>% 
  write_tsv('RDA_WGP_HEAD_F1_F2_ECO_GENES.txt')


head_RDA_TGP_f1_f2 = inner_join(head_tgp,
                                F1_F2_eco_res,
                                by = join_by(Chromosome,
                                             POS <= win_end,
                                             POS >= win_start),
                                suffix = c("_RDA", "_Meth"),
                                relationship = "many-to-many") %>% 
  mutate(win_start_loc = POS-1, 
         win_end_loc = POS+1)




setDT(head_RDA_TGP_f1_f2)
setDT(gene_annotation)

setkey(head_RDA_TGP_f1_f2, 
       Chromosome, 
       win_start_loc, 
       win_end_loc)


head_RDA_TGP_f1_f2_gene_overlap = foverlaps(gene_annotation,
                                            head_RDA_TGP_f1_f2,
                                            # by.x = start,
                                            # by.y = end,
                                            type="any")


head_RDA_TGP_f1_f2_gene_overlap = as_tibble(head_RDA_TGP_f1_f2_gene_overlap) %>% 
  na.omit() %>%
  # filter(Chromosome != 'chrUn') %>% 
  arrange(Chromosome, 
          i.win_mid)

head_RDA_TGP_f1_f2_gene_overlap %>% 
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
  dplyr::select(gene_name) %>% 
  unique() %>% 
  write_tsv('RDA_TGP_HEAD_F1_F2_ECO_GENES.txt')

# Body wgp ----------------------------------------------------------------


body_wgp = read_csv('~/Methylation_data/BODY_WGP_clean_RDA_CAND_corr.csv') %>% 
  separate(col = loc, 
           into = c('Chromosome', 
                    'POS'), 
           sep = '-') %>% 
  mutate(wgp_HEAD_clean_great = abs(BODY_WGP_clean) > abs(eco_num) & abs(BODY_WGP_clean) > abs(csize)) %>% 
  filter(wgp_HEAD_clean_great == 'TRUE') %>% 
  # filter(Association == 'BODY_WGP_clean') %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  mutate(POS = as.numeric(POS))



F2_eco_res = read_tsv('Methylation_F2_eco_effect_win1000_step1000.txt') %>% 
  mutate(Effect = 'F2', 
         abs_effect = abs(expression_mean)) %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  filter(expression_n > 3) %>%
  group_by(Chromosome) %>%
  mutate(threshold = quantile(abs_effect,
                              0.99,
                              na.rm = TRUE)) %>%
  filter(abs_effect >= threshold) %>% 
  mutate(direction = case_when(
    expression_mean > threshold ~ "Hyper",
    expression_mean < -threshold ~ "Hypo",
    TRUE ~ "None")) %>%
  # Order by genomic position
  arrange(Chromosome, win_mid)


body_RDA_WGP_Ecotype = inner_join(body_wgp,
                                  F2_eco_res,
                                  by = join_by(Chromosome,
                                               POS <= win_end,
                                               POS >= win_start),
                                  suffix = c("_RDA", "_Meth"),
                                  relationship = "many-to-many") %>% 
  mutate(win_start_loc = POS-1, 
         win_end_loc = POS+1)



setDT(body_RDA_WGP_Ecotype)
setDT(gene_annotation)

setkey(body_RDA_WGP_Ecotype, 
       Chromosome, 
       win_start_loc, 
       win_end_loc)


body_RDA_WGP_Ecotype_gene_overlap = foverlaps(gene_annotation,
                                              body_RDA_WGP_Ecotype,
                                              # by.x = start,
                                              # by.y = end,
                                              type="any")


body_RDA_WGP_Ecotype_gene_overlap = as_tibble(body_RDA_WGP_Ecotype_gene_overlap) %>% 
  na.omit() %>%
  # filter(Chromosome != 'chrUn') %>% 
  arrange(Chromosome, 
          i.win_mid)

body_RDA_WGP_Ecotype_gene_overlap %>% 
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
  dplyr::select(gene_name) %>% 
  unique() %>% 
write_tsv('RDA_body_wgp_ECO_GENES.txt')


# Body tgp ecotype----------------------------------------------------------------


body_tgp = read_csv('~/Methylation_data/BODY_TGP_clean_RDA_CAND_corr.csv') %>% 
  separate(col = loc, 
           into = c('Chromosome', 
                    'POS'), 
           sep = '-') %>% 
  mutate(tgp_HEAD_clean_great = abs(BODY_TGP_clean) > abs(eco_num) & abs(BODY_TGP_clean) > abs(csize)) %>% 
  filter(tgp_HEAD_clean_great == 'TRUE') %>% 
  # filter(Association == 'BODY_TGP_clean') %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  mutate(POS = as.numeric(POS))



F1_eco_res = read_tsv('Methylation_F1_eco_effect_win1000_step1000.txt') %>% 
  mutate(Effect = 'F1', 
         abs_effect = abs(expression_mean)) %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  filter(expression_n > 3) %>%
  group_by(Chromosome) %>%
  mutate(threshold = quantile(abs_effect,
                              0.99,
                              na.rm = TRUE)) %>%
  filter(abs_effect >= threshold) %>% 
  mutate(direction = case_when(
    expression_mean > threshold ~ "Hyper",
    expression_mean < -threshold ~ "Hypo",
    TRUE ~ "None")) %>%
  # Order by genomic position
  arrange(Chromosome, win_mid)


body_RDA_tgp_Ecotype = inner_join(body_tgp,
                                  F1_eco_res,
                                  by = join_by(Chromosome,
                                               POS <= win_end,
                                               POS >= win_start),
                                  suffix = c("_RDA", "_Meth"),
                                  relationship = "many-to-many") %>% 
  mutate(win_start_loc = POS-1, 
         win_end_loc = POS+1)


# WGP whole fish ----------------------------------------------------------
wgp_whole_fish = read_csv('~/Methylation_data/WGP2_clean_RDA_CAND_corr.csv') %>% 
  separate(col = loc, 
           into = c('Chromosome', 
                    'POS'), 
           sep = '-') %>% 
  mutate(wgp_whole_clean_great = abs(WGP_clean) > abs(eco_num)) %>% 
  filter(wgp_whole_clean_great == 'TRUE') %>% 
  # filter(Association == 'WGP_clean') %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  mutate(POS = as.numeric(POS))



F2_eco_res = read_tsv('Methylation_F2_eco_effect_win1000_step1000.txt') %>% 
  mutate(Effect = 'F2', 
         abs_effect = abs(expression_mean)) %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  filter(expression_n > 3) %>%
  group_by(Chromosome) %>%
  mutate(threshold = quantile(abs_effect,
                              0.99,
                              na.rm = TRUE)) %>%
  filter(abs_effect >= threshold) %>% 
  mutate(direction = case_when(
    expression_mean > threshold ~ "Hyper",
    expression_mean < -threshold ~ "Hypo",
    TRUE ~ "None")) %>%
  # Order by genomic position
  arrange(Chromosome, win_mid)


whole_RDA_wgp_Ecotype = inner_join(wgp_whole_fish,
                                  F2_eco_res,
                                  by = join_by(Chromosome,
                                               POS <= win_end,
                                               POS >= win_start),
                                  suffix = c("_RDA", "_Meth"),
                                  relationship = "many-to-many") %>% 
  mutate(win_start_loc = POS-1, 
         win_end_loc = POS+1)


setDT(whole_RDA_wgp_Ecotype)
setDT(gene_annotation)

setkey(whole_RDA_wgp_Ecotype, 
       Chromosome, 
       win_start_loc, 
       win_end_loc)


whole_RDA_WGP_Ecotype_gene_overlap = foverlaps(gene_annotation,
                                               whole_RDA_wgp_Ecotype,
                                              # by.x = start,
                                              # by.y = end,
                                              type="any")


whole_RDA_WGP_Ecotype_gene_overlap = as_tibble(whole_RDA_WGP_Ecotype_gene_overlap) %>% 
  na.omit() %>%
  # filter(Chromosome != 'chrUn') %>% 
  arrange(Chromosome, 
          i.win_mid)

whole_RDA_WGP_Ecotype_gene_overlap %>% 
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
  dplyr::select(gene_name) %>% 
  unique() %>% 
  write_tsv('RDA_wholefish_wgp_ECO_GENES.txt')



# TGP whole fish ----------------------------------------------------------


tgp_whole_fish = read_csv('~/Methylation_data/TGP2_clean_RDA_CAND_corr.csv') %>% 
  separate(col = loc, 
           into = c('Chromosome', 
                    'POS'), 
           sep = '-') %>% 
  mutate(tgp_whole_clean_great = abs(TGP_clean) > abs(eco_num)) %>% 
  filter(tgp_whole_clean_great == 'TRUE') %>% 
  # filter(Association == 'TGP_clean') %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  mutate(POS = as.numeric(POS))



F1_eco_res = read_tsv('Methylation_F1_eco_effect_win1000_step1000.txt') %>% 
  mutate(Effect = 'F1', 
         abs_effect = abs(expression_mean)) %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  filter(expression_n > 3) %>%
  group_by(Chromosome) %>%
  mutate(threshold = quantile(abs_effect,
                              0.99,
                              na.rm = TRUE)) %>%
  filter(abs_effect >= threshold) %>% 
  mutate(direction = case_when(
    expression_mean > threshold ~ "Hyper",
    expression_mean < -threshold ~ "Hypo",
    TRUE ~ "None")) %>%
  # Order by genomic position
  arrange(Chromosome, win_mid)


whole_RDA_tgp_Ecotype = inner_join(tgp_whole_fish,
                                   F2_eco_res,
                                   by = join_by(Chromosome,
                                                POS <= win_end,
                                                POS >= win_start),
                                   suffix = c("_RDA", "_Meth"),
                                   relationship = "many-to-many") %>% 
  mutate(win_start_loc = POS-1, 
         win_end_loc = POS+1)


setDT(whole_RDA_tgp_Ecotype)
setDT(gene_annotation)

setkey(whole_RDA_tgp_Ecotype, 
       Chromosome, 
       win_start_loc, 
       win_end_loc)


whole_RDA_tgp_Ecotype_gene_overlap = foverlaps(gene_annotation,
                                               whole_RDA_tgp_Ecotype,
                                               # by.x = start,
                                               # by.y = end,
                                               type="any")


whole_RDA_tgp_Ecotype_gene_overlap = as_tibble(whole_RDA_tgp_Ecotype_gene_overlap) %>% 
  na.omit() %>%
  # filter(Chromosome != 'chrUn') %>% 
  arrange(Chromosome, 
          i.win_mid)

whole_RDA_tgp_Ecotype_gene_overlap %>% 
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
  dplyr::select(gene_name) %>% 
  unique() %>% 
  write_tsv('RDA_wholefish_tgp_ECO_GENES.txt')



# Opercular Functional trait - unmodified ---------------------------------
orig_opercKT_fish = read_csv('~/Methylation_data/func_orig_operc_RDA_CAND_corr.csv') %>% 
  separate(col = loc, 
           into = c('Chromosome', 
                    'POS'), 
           sep = '-') %>% 
  mutate(orig_operc_great = abs(orig_opercular_kt) > abs(eco_num) & abs(orig_opercular_kt) > abs(csize)) %>% 
  filter(orig_operc_great == 'TRUE') %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  mutate(POS = as.numeric(POS))


eco_res  = read_tsv('Methylation_eco_effect_win1000_step1000.txt') %>% 
  mutate(Effect = 'Eco', 
         abs_effect = abs(expression_mean)) %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  filter(expression_n > 3) %>%
  group_by(Chromosome) %>%
  mutate(threshold = quantile(abs_effect,
                              0.99,
                              na.rm = TRUE)) %>%
  filter(abs_effect >= threshold) %>% 
  mutate(direction = case_when(
    expression_mean > threshold ~ "Hyper",
    expression_mean < -threshold ~ "Hypo",
    TRUE ~ "None")) %>%
  # Order by genomic position
  arrange(Chromosome, win_mid)


opercKT_RDA_orig_Ecotype = inner_join(orig_opercKT_fish,
                                   eco_res,
                                   by = join_by(Chromosome,
                                                POS <= win_end,
                                                POS >= win_start),
                                   suffix = c("_RDA", "_Meth"),
                                   relationship = "many-to-many") %>% 
  mutate(win_start_loc = POS-1, 
         win_end_loc = POS+1)




# Opercular Functional trait - TGp ----------------------------------------

tgp_opercKT_fish = read_csv('~/Methylation_data/func_TGP_operc_RDA_CAND_corr.csv') %>% 
  separate(col = loc, 
           into = c('Chromosome', 
                    'POS'), 
           sep = '-') %>% 
  # dplyr::select(-Interaction) %>% 
  mutate(TGP_operc_kt_great = TGP_operc_kt > eco_num & TGP_operc_kt > csize) %>% 
  filter(TGP_operc_kt_great == 'TRUE') %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  mutate(POS = as.numeric(POS))

F1_eco_res = read_tsv('Methylation_F1_eco_effect_win1000_step1000.txt') %>% 
  mutate(Effect = 'F1', 
         abs_effect = abs(expression_mean)) %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  filter(expression_n > 3) %>%
  group_by(Chromosome) %>%
  mutate(threshold = quantile(abs_effect,
                              0.99,
                              na.rm = TRUE)) %>%
  filter(abs_effect >= threshold) %>% 
  mutate(direction = case_when(
    expression_mean > threshold ~ "Hyper",
    expression_mean < -threshold ~ "Hypo",
    TRUE ~ "None")) %>%
  # Order by genomic position
  arrange(Chromosome, win_mid)



opercKT_RDA_tgp_Ecotype = inner_join(tgp_opercKT_fish,
                                      F1_eco_res,
                                      by = join_by(Chromosome,
                                                   POS <= win_end,
                                                   POS >= win_start),
                                      suffix = c("_RDA", "_Meth"),
                                      relationship = "many-to-many") %>% 
  mutate(win_start_loc = POS-1, 
         win_end_loc = POS+1)


# Opercular Functional traits - WGP ---------------------------------------


wgp_opercKT_fish = read_csv('~/Methylation_data/func_WGP_operc_RDA_CAND_corr.csv') %>% 
  separate(col = loc, 
           into = c('Chromosome', 
                    'POS'), 
           sep = '-') %>% 
  # dplyr::select(-Interaction) %>% 
  mutate(wgp_operc_kt_great = WGP_operc_kt > eco_num & WGP_operc_kt > csize) %>% 
  filter(wgp_operc_kt_great == 'TRUE') %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  mutate(POS = as.numeric(POS))

F2_eco_res = read_tsv('Methylation_F2_eco_effect_win1000_step1000.txt') %>% 
  mutate(Effect = 'F2', 
         abs_effect = abs(expression_mean)) %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  filter(expression_n > 3) %>%
  group_by(Chromosome) %>%
  mutate(threshold = quantile(abs_effect,
                              0.99,
                              na.rm = TRUE)) %>%
  filter(abs_effect >= threshold) %>% 
  mutate(direction = case_when(
    expression_mean > threshold ~ "Hyper",
    expression_mean < -threshold ~ "Hypo",
    TRUE ~ "None")) %>%
  # Order by genomic position
  arrange(Chromosome, win_mid)



opercKT_RDA_wgp_Ecotype = inner_join(wgp_opercKT_fish,
                                     F2_eco_res,
                                     by = join_by(Chromosome,
                                                  POS <= win_end,
                                                  POS >= win_start),
                                     suffix = c("_RDA", "_Meth"),
                                     relationship = "many-to-many") %>% 
  mutate(win_start_loc = POS-1, 
         win_end_loc = POS+1)

