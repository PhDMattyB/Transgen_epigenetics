## RDA results for TGP and WGP

library(tidyverse)




# Head WGP ----------------------------------------------------------------


head_wgp = read_csv('~/Methylation_data/WGP_HEAD_clean_RDA_CAND_corr.csv') %>% 
  separate(col = loc, 
           into = c('Chromosome', 
                    'POS'), 
           sep = '-') %>% 
  filter(Association == 'WGP_HEAD_clean') %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  mutate(POS = as.numeric(POS))


# head_wgp %>%
#   group_by(Association) %>% 
#   summarize(n = n())

# head_wgp %>% 
#   filter(Chromosome == 'chrXXI') %>% 
#   arrange(POS) %>% View()

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
           relationship = "many-to-many") 

# Head tgp ----------------------------------------------------------------



head_tgp = read_csv('~/Methylation_data/TGP_HEAD_clean_RDA_CAND_corr.csv') %>% 
  separate(col = loc, 
           into = c('Chromosome', 
                    'POS'), 
           sep = '-')%>% 
  filter(Association == 'TGP_HEAD_clean') %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  mutate(POS = as.numeric(POS))

# head_tgp %>%
#   group_by(Association) %>% 
#   summarize(n = n())
# 
# head_tgp %>% 
#   filter(Chromosome == 'chrXXI') %>% 
#   arrange(POS) 

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



head_RDA_TGP_f1_f2 = inner_join(head_tgp,
           F1_F2_eco_res,
           by = join_by(Chromosome,
                        POS <= win_end,
                        POS >= win_start),
           suffix = c("_RDA", "_Meth"),
           relationship = "many-to-many") 

head_RDA_WGP_f1_f2 = inner_join(head_wgp,
           F1_F2_eco_res,
           by = join_by(Chromosome,
                        POS <= win_end,
                        POS >= win_start),
           suffix = c("_RDA", "_Meth"),
           relationship = "many-to-many") 


head_RDA_WGP_Ecotype %>% 
  group_by(Chromosome, 
           POS) %>% 
  intersect(., head_RDA_WGP_f1_f2)
