### LIMMA test run
## Full model matrix


library(tidyverse)
library(limma)
library(edgeR)

setwd('~/Methylation_data/')

setwd('C:/Users/mkb6d/Documents/')

# mvalues = read_csv('MVALUES_methylation_cleaned_data.csv')
# 
# meta_data = read_csv('formattedDataEU.csv') %>%
#   select(fish,
#          ID,
#          Full_ID,
#          F1,
#          F2,
#          poppair,
#          ecotype...13,
#          csize_real) %>%
#   rename(ecotype = ecotype...13)
# 
mval_small = read_csv('mvalues_chrI.csv')
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

  
  # separate(Fish_ID,
  #          into = c('ID',
  #                   'ID2',
  #                   'ID3'),
  #          sep = '_',
  #          remove = F) %>%
  # separate(ID3,
  #          into = c('symbol',
  #                   'ID4',
  #                   'num'),
  #          sep = '',
  #          remove = F) %>%
  # select(-symbol) %>%
  # mutate(new_vals = paste0('#G',
  #                          str_pad(num,
  #                                  0,
  #                                  side = 'left'))) %>%
  # select(-ID3,
  #        -ID4,
  #        -num) %>%
  # unite(Fish_ID2,
  #       c('ID',
  #         'ID2',
  #         'new_vals'),
  #       sep = '_') %>%
  # select(Fish_ID2) %>%
  # rename(Fish_ID = Fish_ID2)

# meta_data_cleaned = bind_cols(meta_fish_ID,
#                               meta_data) %>%
#   arrange(Fish_ID)%>%
#   filter(grepl('#G', Fish_ID))


# 
meta_data_cleaned = read_csv('Cleaned_meta_data.csv') %>%
  arrange(Fish_ID)

mval_small_ID = mval_small %>% 
  select(Location_data)

methy_fish_ID = mval_small_ID %>% 
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
        sep = '_') %>% 
  arrange(Fish_ID) 

mvals_cleaned = bind_cols(methy_fish_ID, 
                          mval_small) %>% 
  select(-Location_data) %>% 
  arrange(Fish_ID)


meta_data_cleaned = meta_data_cleaned %>% 
  separate(Fish_ID, 
           into = c('G', 
                    'num'), 
           sep = '#') %>% 
  filter(num %in% c('G1', 
                    'G2', 
                    'G3', 
                    'G4', 
                    'G5')) %>% 
  unite(Fish_ID, 
        c('G', 
          'num'), 
        sep = '#') %>% 
  arrange(Fish_ID)

methy_fish_ID = inner_join(meta_data_cleaned, 
           methy_fish_ID, 
           by = 'Fish_ID')

## organized data. 

## Finally cleaned
methy_fish_ID = read_csv("Methylation_meta_data_issue.csv") %>% 
  arrange(Fish_ID)


methy_fish_ID$F1 = factor(methy_fish_ID$F1)
methy_fish_ID$F2 = factor(methy_fish_ID$F2)
methy_fish_ID$ecotype = factor(methy_fish_ID$ecotype)
methy_fish_ID$poppair = factor(methy_fish_ID$poppair)

design_formula = ~0+F1*F2*ecotype*poppair

design = model.matrix(design_formula, data = methy_fish_ID)


design_formula2 = ~0+F1*F2*ecotype

design2 = model.matrix(design_formula2, data = methy_fish_ID)

# print(design)

set.seed(1738)

trans_methy = mvals_cleaned %>%
  select(-Fish_ID) %>% 
  # remove_rownames() %>% 
  # select(-fish, 
  #        -ID, 
  #        -Full_ID, 
  #        -F1, 
  #        -F2, 
  #        -poppair, 
  #        -ecotype, 
  #        -csize_real) %>% 
  # row.names(.)
  # rownames_to_column(var = 'Fish_ID') %>% 
  t() %>% 
  as_tibble()
  
mvals_cleaned %>% 
  select(1) %>% 
  write_tsv('mvals_small_colnames.txt')

names(mval_small) %>% 
  as_tibble() %>% 
  slice(-1) %>% 
  write_tsv('mvals_small_chr_rownames.txt')

fit = lmFit(trans_methy, design2)

ebayes_fit = eBayes(fit)

topTable(ebayes_fit) %>%
  rownames_to_column() %>% 
  as_tibble() %>% View()
  write_csv('Refined_model_Top_table_example_df_CHRI_methylation.csv')

  
ordinary.t = fit$coefficients/fit$stdev.unscaled/fit$sigma 


head(ordinary.t)

# ebayes_fit$p.value %>% 
#   as.data.frame() %>% 
#   as_tibble() %>% 
#   arrange(`F118:F218:ecotypew`)
