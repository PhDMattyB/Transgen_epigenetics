
# RDA FUNCTIONAL TRAITS ---------------------------------------------------

## 07.03.2025

## Matt Brachmann

## PhDMattyb

# Functions ---------------------------------------------------------------

outliers = function(x,z){
  lims = mean(x) + c(-1,1)*z*sd(x)
  
  x[x<lims[1] | x>lims[2]]
}



# start up ----------------------------------------------------------------

setwd('~/Methylation_data/')

library(tidyverse)
library(sjmisc)
library(vegan)
library(patchwork)

theme_set(theme_bw())

meta_data = read_csv('formattedDataEU.csv') %>% 
  select(fish, 
         ID, 
         Full_ID, 
         F1, 
         F2, 
         poppair, 
         ecotype...13, 
         csize_real) %>% 
  rename(ecotype = ecotype...13) 
# %>% 


# TGP clean traits ---------------------------------------------------
TGP_clean_data = read_csv('tgp_wgp_CandiscCleaned_data.csv') %>% 
  inner_join(.,
             meta_data, 
           by = 'fish', 
           relationship = 'many-to-many')


TGP_pheno_fish = TGP_clean_data %>% 
  filter(str_detect(fish,
                    '_#G')) %>% 
  distinct(fish, 
           .keep_all = T) %>% 
  select(fish)

TGP_pheno_fish_ID = TGP_pheno_fish %>% 
  separate_wider_regex(fish, 
                       c(var1 = ".*?", 
                         "_", 
                         var2 = ".*")) %>% 
  separate_wider_regex(var2, 
                       c(f1_temp = ".*?", 
                         "@", 
                         f2_temp = ".*")) %>% 
  separate(f2_temp, 
           into = c('f2_temp', 
                    'trash'), 
           sep = '[.]') %>% 
  select(var1, 
         f1_temp, 
         f2_temp) 

TGP_pheno_fish_ID$var1 = gsub("'", '', TGP_pheno_fish_ID$var1)

TGP_pheno_fish_ID = TGP_pheno_fish_ID %>% 
  unite(col = 'Fish_ID', 
        sep = '')%>% 
  mutate(Fish_ID = gsub("Myvat", 
                        "MYV", 
                        Fish_ID)) 


TGP_pheno_fish = TGP_clean_data %>%
  filter(str_detect(fish,
                    '_#G')) %>% 
  distinct(fish, 
           .keep_all = T) %>%  
           # dplyr::select(Full_ID) %>%
           bind_cols(.,
                     TGP_pheno_fish_ID)
           


mvalues = read_csv('MVALUES_methylation_cleaned_data.csv')

meth_fish = mvalues %>%
  select(1)
# 
# meth_fish = mvalues %>% 
#   select(1)
# 
meth_fish_ID = meth_fish %>% 
  # separate(Location_data, 
  #         into = c('garbage', 
  #            'ID'), 
  #          sep = '-') %>% 
  separate_wider_regex(Location_data, 
                       c(var1 = ".*?", 
                         "-", 
                         var2 = ".*")) %>% 
  select(var2) %>% 
  separate(var2, 
           into = c('ID', 
                    'ID2'), 
           sep = '-') %>% 
  mutate(new_vals = paste0('#G', 
                           str_pad(ID2, 
                                   0, 
                                   side = 'left'))) %>% 
  select(ID, 
         new_vals) %>% 
  unite(col = Fish_ID, 
        sep = '_')




## identifying potential issues with phenotypic data
## 5 individuals that were sequenced did not have the #G identifier
## Four from GTS18@12 and one from MYVC12@12
anti_join(meth_fish_ID,
          TGP_pheno_fish_ID)


TGP_clean_pheno_fish_final = inner_join(meth_fish_ID, 
                                        TGP_pheno_fish) %>% 
  arrange(Fish_ID)

## now we need to order the phenotypic data and methylation data
## they have to be in the same order otherwise this will all
## be fucked. Right now they aren't. 
## order the data by the Fish_ID column

## Add the new id column to the methylation data

TGP_clean_mvalues_final = bind_cols(meth_fish_ID, 
                                    mvalues) %>% 
  arrange(Fish_ID)  %>% 
  filter(Fish_ID != 'GTS1812_EU1_#G1', 
         Fish_ID != 'GTS1812_EU1_#G2',
         Fish_ID != 'GTS1812_EU1_#G3',
         Fish_ID != 'GTS1812_EU1_#G4', 
         Fish_ID != 'MYVC1212_EU4_#G1')

TGP_clean_mvalues_final$Fish_ID == TGP_clean_pheno_fish_final$Fish_ID


## need to check that everythings in order for the analyses
## If there are any FALSE we're fucked. 
## Shooting for all TRUES

# setdiff(TGP_clean_mvalues_final$Fish_ID, 
#         TGP_clean_pheno_fish_final$Fish_ID)
## TGP RDA 
TGP_clean_mvalues_only = TGP_clean_mvalues_final %>% 
  # select(-1) %>% 
  select(-Location_data, 
         -Fish_ID)

TGP_clean_RDA_full = rda(TGP_clean_mvalues_only ~ TGPclean + ecotype + TGPclean*ecotype + csize_real, 
                        data = TGP_clean_pheno_fish_final, 
                        scale = T)

TGP_clean_eco_RDA2 = rda(TGP_clean_mvalues_only ~ TGPclean + F1text,
                         data = TGP_clean_pheno_fish_final,
                         scale = T)

RsquareAdj(TGP_clean_eco_RDA)
summary(eigenvals(TGP_clean_eco_RDA, 
                  model = 'constrained'))

screeplot(TGP_clean_eco_RDA)

## Run after all other coding is finished
## This will take a while 
TGP_eco_signif_full = anova.cca(TGP_clean_eco_RDA, 
                                parallel = getOption('mc.cores'))

TGP_eco_signif_full2 = anova.cca(TGP_clean_eco_RDA2,
                                 parallel = getOption('mc.cores'))

vif.cca(TGP_clean_eco_RDA2)

TGP_clean_sum_rda2 = summary(TGP_clean_eco_RDA2)

TGP_clean_sum_rda2$species %>%
  as_tibble() %>%
  write_csv('TGP_clean_RDA_Uncorrected_PCA_locations.csv')

TGP_clean_sum_rda2$sites %>%
  as_tibble() %>%
  write_csv('TGP_clean_RDA_Uncorrected_PCA_individuals.csv')

TGP_clean_sum_rda2$biplot %>%
  as_tibble() %>%
  write_csv('TGP_clean_RDA_Uncorrected_PCA_biplot.csv')


TGP_clean_rda_scores2 = scores(TGP_clean_eco_RDA2,
                               choices = c(1:5),
                               display = 'species')
# 
hist(TGP_clean_rda_scores2[,1])
hist(TGP_rda_scores[,2])
hist(TGP_rda_scores[,3])
hist(TGP_rda_scores[,4])
hist(TGP_rda_scores[,5])
# 
TGP2_clean_rda_outliers_axis1 = outliers(TGP_clean_rda_scores2[,1], 3)
# 
# 
TGP2_clean_rda_out_axis1 = cbind.data.frame(rep(1,
                                                times = length(TGP2_clean_rda_outliers_axis1)),
                                            names(TGP2_clean_rda_outliers_axis1),
                                            unname(TGP2_clean_rda_outliers_axis1))
# 
TGP2_clean_rda_out_axis1 = TGP2_clean_rda_out_axis1 %>%
  as_tibble() %>%
  dplyr::rename(axis = 1,
                loc = 2,
                scores = 3)
# 
# 
TGP_all_loc2 = TGP_clean_rda_scores2[,1]
# 
TGP2_rda_normal = cbind.data.frame(rep(2,
                                       times = length(TGP_all_loc2)),
                                   names(TGP_all_loc2),
                                   unname(TGP_all_loc2))
TGP2_rda_normal = TGP2_rda_normal %>%
  as_tibble() %>%
  dplyr::rename(axis = 1,
                loc = 2,
                scores = 3)
# 
TGP2_rda_normal = TGP2_rda_normal[!TGP2_rda_normal$loc %in% TGP2_clean_rda_out_axis1$loc,]
# 
write_csv(TGP2_rda_normal,
          'TGP2_clean_RDA_PCaxes_nonoutliers_methylation.csv')

write_csv(TGP2_clean_rda_out_axis1,
          'TGP2_clean_RDA_outliers_AXIS1_RAW_PCaxes_methylation.csv')
# 
TGP2_clean_rda_out = as.data.frame(TGP2_clean_rda_out_axis1)
TGP2_all_loc = as.data.frame(TGP_all_loc2)
TGP_clean_pheno_fish_final2 = as.data.frame(TGP_clean_pheno_fish_final)
# 
TGP_clean_phenotypes = TGP_clean_pheno_fish_final2 %>%
  as_tibble() %>%
  mutate(eco_num = as.numeric(case_when(
    F1text == 'cold' ~ '1',
    F1text == 'warm' ~ '2'))) %>%
  # dplyr::select(-Population) %>%
  as.data.frame()

# # nam = rda_out[1:45, 2]
# # out_loc = all_loc[nam,]
# # out_cor = apply(test_pheno,
# #                 2, 
# #                 function(x)cor(x, out_loc))
# 
foo = matrix(nrow=(948),
             ncol = 2)
# colnames(foo) = c('TGP_clean',
#                   'eco_num', 
#                   'interaction')
colnames(foo) = c('TGP_clean',
                  'eco_num')

TGP_clean_phenotypes = TGP_clean_phenotypes %>% 
  dplyr::select(TGPclean, 
                eco_num) 
# %>% 
#   mutate(interaction = TGPclean*eco_num)

for (i in 1:length(TGP2_clean_rda_out$loc)){
  nam = TGP2_clean_rda_out[i,2]
  loc.gen = TGP_clean_mvalues_only[,nam]
  foo[i,] = apply(TGP_clean_phenotypes,2,function(x)cor(x,loc.gen))
}

TGP2_clean_candidates = cbind.data.frame(TGP2_clean_rda_out,
                                         foo)

TGP2_clean_candidates %>%
  as_tibble() %>%
  write_csv('TGP2_clean_RDA_outliers_methylation_correlations.csv')
# 
##check for duplicates
length(TGP2_clean_candidates$loc[duplicated(TGP2_clean_candidates$loc)])
# 
for(i in 1:length(TGP2_clean_candidates$loc)){
  bar = TGP2_clean_candidates[i,]
  TGP2_clean_candidates[i,6] = names(which.max(abs(bar[4:5])))
  TGP2_clean_candidates[i,7] = max(abs(bar[4:5]))
}
# 
TGP2_clean_candidates %>% 
  as_tibble() %>%
  rename(Association = V6, 
         Correlation = V7) %>%
  write_csv("TGP2_clean_RDA_CAND_corr.csv") %>% 
  group_by(Association) %>% 
  summarize(n = n())

