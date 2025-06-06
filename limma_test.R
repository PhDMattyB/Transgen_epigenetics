### LIMMA test run
## Full model matrix


library(tidyverse)
library(limma)
library(edgeR)
library(qvalue)
library(glmmTMB)

setwd('~/Methylation_data/')


mvalues = read_csv('MVALUES_methylation_cleaned_data.csv') %>% 
  select(-Location_data)
meta_data = read_csv('Methylation_meta_data.csv')%>% 
  mutate(ecotype = as.character(case_when(
    poppair == 'GTS' ~ 'W',
    poppair == 'CSWY' ~ 'C', 
    poppair == 'SKRC' ~ 'W', 
    poppair == 'SKRW' ~ 'C', 
    poppair == 'MYVC' ~ 'C', 
    poppair == 'MYVW' ~ 'W', 
    poppair ==  'ASHNC' ~ 'C', 
    poppair == 'ASHNW' ~ 'W'))) %>% 
  mutate(poppair = as.character(case_when(
    poppair == 'GTS' ~ 'GTS',
    poppair == 'CSWY' ~ 'CSWY', 
    poppair == 'SKRC' ~ 'SKR', 
    poppair == 'SKRW' ~ 'SKR', 
    poppair == 'MYVC' ~ 'MYV', 
    poppair == 'MYVW' ~ 'MYV', 
    poppair ==  'ASHNC' ~ 'ASHN', 
    poppair == 'ASHNW' ~ 'ASHN'
  )))


meta_data$F1_temp = as.character(meta_data$F1_temp)
meta_data$F2_temp = as.character(meta_data$F2_temp)



## making the metadata
# 
# mval_small_meta %>%
#   separate(col = Location_data, 
#            into = c('group1', 
#                     'group2', 
#                     'group3'), 
#            sep = '_', 
#            remove = F) %>% 
#   dplyr::select(-group1) %>% 
#   separate(col = group2, 
#            into = c('individualID', 
#                     'pop_temp'), 
#            sep = '-', 
#            remove = F) %>% 
#   separate(col = group3, 
#            into = c('EU_ID', 
#                     'EU_individualID', 
#                     'garbarge'), 
#            sep = '-') %>% 
#   select(-garbarge, 
#          -group2) %>%
#   extract(pop_temp, 
#           into = c('poppair', 
#                    'temps'), 
#           '^([A-Z]+)(\\d+)$') %>% 
#   mutate(len = str_length(temps), 
#          half = len %/% 2, 
#          F1_temp = str_sub(temps, 1, half), 
#          F2_temp = str_sub(temps, 
#                            half + 1, 
#                            len)) %>% 
#   dplyr::select(Location_data, 
#                 individualID, 
#                 poppair, 
#                 temps, 
#                 EU_ID, 
#                 EU_individualID, 
#                 F1_temp, 
#                 F2_temp) %>% 
#   write_csv('Methylation_meta_data.csv')
# 
# 
# meta_fish_ID = meta_data %>%
#   separate_wider_regex(fish,
#                        c(var1 = ".*?",
#                          "_",
#                          var2 = ".*")) %>%
#   separate_wider_regex(var2,
#                        c(f1_temp = ".*?",
#                          "@",
#                          f2_temp = ".*")) %>%
#    separate(f2_temp,
#            into = c('f2_temp',
#                     'trash'),
#            sep = '[.]') %>%
#   select(var1,
#          f1_temp,
#          f2_temp)
# 
# meta_fish_ID$var1 = gsub("'", '', meta_fish_ID$var1)
# 
# meta_fish_ID = meta_fish_ID %>%
#   unite(col = 'Fish_ID',
#         sep = '')%>%
#   mutate(Fish_ID = gsub("Myvat",
#                         "MYV",
#                         Fish_ID)) 
# 




## THIS SHOULD NOT BE NEEDED. 
## THIS WAS TO TRY SOMETHING ELSE
## WHERE I HAD TO ROTATE THE DATA SET
# trans_methy = mvals_cleaned %>%
#   select(-Fish_ID) %>% 
#   # remove_rownames() %>% 
#   # select(-fish, 
#   #        -ID, 
#   #        -Full_ID, 
#   #        -F1, 
#   #        -F2, 
#   #        -poppair, 
#   #        -ecotype, 
#   #        -csize_real) %>% 
#   # row.names(.)
#   # rownames_to_column(var = 'Fish_ID') %>% 
#   t() %>% 
#   as_tibble()
  
# mvals_cleaned %>% 
#   select(1) %>% 
#   write_tsv('mvals_small_colnames.txt')
# 
# names(mval_small) %>% 
#   as_tibble() %>% 
#   slice(-1) %>% 
#   write_tsv('mvals_small_chr_rownames.txt')
# 



# analysis  -----------------------------------------------------------------
set.seed(1738)



## make sure that the rows of the model results table
## are the same as the rows in the mvalues data frame

## the number of columns does not change with the output
model_results_table = as.data.frame(matrix(nrow = 2828701, 
                                          ncol = 34))


for(i in 1:ncol(mvalues)){
  gene = mvalues[,i]
  # all_others = rowSums(test_mvals[,-i])
  # Y = cbind(gene, 
            # all_others)
  
  Model = glmmTMB(gene ~ F1_temp * F2_temp * ecotype + (1|poppair), 
                # family = gaussian(), 
                # control = glmerControl(
                  # optimizer = 'optimx', optCtrl = list(method = 'nlminb')), 
                data = meta_data)
  
  model_results = summary(Model)
  
  
  model_results_table[i, 1] = colnames(test_mvals)[i]
  model_results_table[i, 2] = sum(test_mvals[,i]/sum(test_mvals))
  model_results_table[i, 3] = model_results$coefficients$cond[1,1]
  model_results_table[i, 4] = model_results$coefficients$cond[2,1]
  model_results_table[i, 5] = model_results$coefficients$cond[3,1]
  model_results_table[i, 6] = model_results$coefficients$cond[4,1]
  model_results_table[i, 7] = model_results$coefficients$cond[5,1]
  model_results_table[i, 8] = model_results$coefficients$cond[6,1]
  model_results_table[i, 9] = model_results$coefficients$cond[7,1]
  model_results_table[i, 10] = model_results$coefficients$cond[8,1]
  model_results_table[i, 11] = model_results$coefficients$cond[1,2]
  model_results_table[i, 12] = model_results$coefficients$cond[2,2]
  model_results_table[i, 13] = model_results$coefficients$cond[3,2]
  model_results_table[i, 14] = model_results$coefficients$cond[4,2]
  model_results_table[i, 15] = model_results$coefficients$cond[5,2]
  model_results_table[i, 16] = model_results$coefficients$cond[6,2]
  model_results_table[i, 17] = model_results$coefficients$cond[7,2]
  model_results_table[i, 18] = model_results$coefficients$cond[8,2]
  model_results_table[i, 19] = model_results$coefficients$cond[1,3]
  model_results_table[i, 20] = model_results$coefficients$cond[2,3]
  model_results_table[i, 21] = model_results$coefficients$cond[3,3]
  model_results_table[i, 22] = model_results$coefficients$cond[4,3]
  model_results_table[i, 23] = model_results$coefficients$cond[5,3]
  model_results_table[i, 24] = model_results$coefficients$cond[6,3]
  model_results_table[i, 25] = model_results$coefficients$cond[7,3]
  model_results_table[i, 26] = model_results$coefficients$cond[8,3]
  model_results_table[i, 27] = model_results$coefficients$cond[1,4]
  model_results_table[i, 28] = model_results$coefficients$cond[2,4]
  model_results_table[i, 29] = model_results$coefficients$cond[3,4]
  model_results_table[i, 30] = model_results$coefficients$cond[4,4]
  model_results_table[i, 31] = model_results$coefficients$cond[5,4]
  model_results_table[i, 32] = model_results$coefficients$cond[6,4]
  model_results_table[i, 33] = model_results$coefficients$cond[7,4]
  model_results_table[i, 34] = model_results$coefficients$cond[8,4]
  
  
}



# model results table -----------------------------------------------------

model_results_table = read_csv('model_results_glm_methylation.csv') %>% 
  rename(Methy_loc = V1, 
         expression = V2,
         estimate_intercept = V3,
         estimate_F1 = V4, 
         estimate_F2 = V5, 
         estimate_ecotype = V6, 
         estimate_F1_F2 = V7, 
         estimate_F1_eco = V8, 
         estimate_F2_eco = V9, 
         estimate_F1_F2_eco = V10,
         std_err_intercept = V11,
         std_err_F1 = V12, 
         std_err_F2 = V13, 
         std_err_eco = V14, 
         std_err_F1_F2 = V15,
         std_err_F1_eco = V16, 
         std_err_F2_eco = V17, 
         std_err_F1_F2_eco = V18,
         zval_intercept = V19,
         zval_F1 = V20, 
         zval_F2 = V21, 
         zval_eco = V22, 
         zval_F1_F2 = V23, 
         zval_F1_eco = V24, 
         zval_F2_eco = V25, 
         zval_F1_F2_eco = V26, 
         pval_intercept = V27, 
         pval_F1 = V28, 
         pval_F2 = V29, 
         pval_eco = V30, 
         pval_F1_F2 = V31, 
         pval_F1_eco = V32, 
         pval_F2_eco = V33, 
         pval_F1_F2_eco = V34) %>% 
  separate(col = Methy_loc, 
           into = c('Chromosome', 
                    'Pos'), 
           sep = '-', 
           remove = F)


### F1 temperature FDR
model_F1 = model_results_table %>% 
  select(Chromosome, 
         Pos, 
         expression,
         pval_F1)
F1_qvalues = qvalue(p = model_F1$pval_F1)  

F1_qvalues = F1_qvalues$qvalues

bind_cols(model_F1, 
          F1_qvalues) %>% 
  rename(F1_qvalues = 5) %>% 
  arrange(F1_qvalues) %>% 
  filter(F1_qvalues <= 0.05) %>% 
  arrange(Chromosome, 
          Pos)

## F2 temperature FDR
model_F2 = model_results_table %>% 
  select(Chromosome, 
         Pos, 
         expression,
         pval_F2)
F2_qvalues = qvalue(p = model_F2$pval_F2)  

F2_qvalues = F2_qvalues$qvalues

bind_cols(model_F2, 
          F2_qvalues) %>% 
  rename(F2_qvalues = 5) %>% 
  arrange(F2_qvalues) %>% 
  filter(F2_qvalues <= 0.05) %>% 
  arrange(Chromosome, 
          Pos)


## Ecotype FDR
model_eco = model_results_table %>% 
  select(Chromosome, 
         Pos, 
         expression,
         pval_eco)
qvalues_eco = qvalue(p = model_eco$pval_eco)  

qvalues_eco = qvalues_eco$qvalues

bind_cols(model_eco, 
          qvalues_eco) %>% 
  rename(qvalues_eco = 5) %>% 
  arrange(qvalues_eco) %>% 
  filter(qvalues_eco <= 0.05) %>% 
  arrange(Chromosome, 
          Pos)


### F1 * F2 FDR
model_F1_F2 = model_results_table %>% 
  select(Chromosome, 
         Pos, 
         expression,
         pval_F1_F2)
qvalues_F1_F2 = qvalue(p = model_F1_F2$pval_F1_F2)  

qvalues_F1_F2 = qvalues_F1_F2$qvalues

bind_cols(model_F1_F2, 
          qvalues_F1_F2) %>% 
  rename(qvalues_F1_F2 = 5) %>% 
  arrange(qvalues_F1_F2) %>% 
  filter(qvalues_F1_F2 <= 0.05) %>% 
  arrange(Chromosome, 
          Pos)

## F1 * eco FDR

model_F1_eco = model_results_table %>% 
  select(Chromosome, 
         Pos, 
         expression,
         pval_F1_eco)
qvalues_F1_eco = qvalue(p = model_F1_eco$pval_F1_eco)  

qvalues_F1_eco = qvalues_F1_eco$qvalues

bind_cols(model_F1_eco, 
          qvalues_F1_eco) %>% 
  rename(qvalues_F1_eco = 5) %>% 
  arrange(qvalues_F1_eco) %>% 
  filter(qvalues_F1_eco <= 0.05) %>% 
  arrange(Chromosome, 
          Pos)


## F2 * eco FDR

model_F2_eco = model_results_table %>% 
  select(Chromosome, 
         Pos, 
         expression,
         pval_F2_eco)
qvalues_F2_eco = qvalue(p = model_F2_eco$pval_F2_eco)  

qvalues_F2_eco = qvalues_F2_eco$qvalues

bind_cols(model_F2_eco, 
          qvalues_F2_eco) %>% 
  rename(qvalues_F2_eco = 5) %>% 
  arrange(qvalues_F2_eco) %>% 
  filter(qvalues_F2_eco <= 0.05) %>% 
  arrange(Chromosome, 
          Pos)


## F1 * F2 * eco FDR

model_F1_F2_eco = model_results_table %>% 
  select(Chromosome, 
         Pos, 
         expression,
         pval_F1_F2_eco)
qvalues_F1_F2_eco = qvalue(p = model_F1_F2_eco$pval_F1_F2_eco)  

qvalues_F1_F2_eco = qvalues_F1_F2_eco$qvalues

bind_cols(model_F1_F2_eco, 
          qvalues_F1_F2_eco) %>% 
  rename(qvalues_F1_F2_eco = 5) %>% 
  arrange(qvalues_F1_F2_eco) %>% 
  filter(qvalues_F1_F2_eco <= 0.05) %>% 
  arrange(Chromosome, 
          Pos)

