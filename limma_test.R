### LIMMA test run
## Full model matrix


library(tidyverse)
library(limma)
setwd('~/Methylation_data/')

# mvalues = read_csv('MVALUES_methylation_cleaned_data.csv')

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

mval_small = read_csv('mvalues_chrI.csv')


meta_fish_ID = meta_data %>% 
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

meta_fish_ID$var1 = gsub("'", '', meta_fish_ID$var1)

meta_fish_ID = meta_fish_ID %>% 
  unite(col = 'Fish_ID', 
        sep = '')%>% 
  mutate(Fish_ID = gsub("Myvat", 
                        "MYV", 
                        Fish_ID)) %>% 
  separate(Fish_ID, 
           into = c('ID', 
                    'ID2', 
                    'ID3'), 
           sep = '_', 
           remove = F) %>%
  separate(ID3, 
           into = c('symbol', 
                    'ID4', 
                    'num'), 
           sep = '', 
           remove = F) %>% 
  select(-symbol) %>% 
  mutate(new_vals = paste0('#G', 
                           str_pad(num, 
                                   0, 
                                   side = 'left'))) %>% 
  select(-ID3, 
         -ID4, 
         -num) %>% 
  unite(Fish_ID2, 
        c('ID', 
          'ID2', 
          'new_vals'), 
        sep = '_') %>% 
  select(Fish_ID2) %>% 
  rename(Fish_ID = Fish_ID2)

meta_data_cleaned = bind_cols(meta_fish_ID, 
                              meta_data)


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
        sep = '_')

mvals_cleaned = bind_cols(methy_fish_ID, 
                          mval_small) %>% 
  select(-Location_data) %>% 
  arrange(Fish_ID)


methy_meta_data = inner_join(meta_data_cleaned, 
           methy_fish_ID, 
           by = 'Fish_ID', 
           relationship = 'many-to-many') %>% 
  arrange(Fish_ID)

mvals_cleaned$Fish_ID == methy_meta_data$Fish_ID


# Final_IDs = inner_join(meta_fish_ID,
#            methy_fish_ID, 
#            by = 'Fish_ID', 
#            relationship = 'many-to-many') %>% 
