##############################
## Transgenerational epigenetics
##
## Matt Brachmann (PhDMattyB)
##
## 07.10.2024
##
##############################


# setwd('C:/Users/phdma/OneDrive/OneDrive - University of Glasgow/Documents/Parsons_Postdoc/Stickleback_Genomic/Transgen_epigenetics/')
setwd('~/Methylation_data/')

library(tidyverse)
library(sjmisc)
# methy_rds = read_rds(file = "Methylated.cov.noZero.rds")
# head(methy_rds)
# 
# methy_rds %>% 
#   # as.matrix() %>% 
#   # rownames_to_column(var = 'Location') %>% 
#   as.data.frame() %>%
#   as_tibble() %>% 
#   # rowid_to_column(var = 'rowid')
#   write_csv('Methylated_cov_nozero.csv')
# 
# 
# row.names(methy_rds) %>% 
#   as.data.frame() %>% 
#   as.tibble() %>% 
#   rename(Location_data = 1) %>% 
#   write_csv('Methylation_location_data.csv')
# 
#   
# unmethy_rds = read_rds(file = 'UnMethylated.cov.noZero.rds')
# 
# unmethy %>% 
#   as.data.frame() %>% 
#   write_csv('Unmethylated_cov_nozero.csv')
# 
# row.names(unmethy_rds) %>% 
#   as.data.frame() %>% 
#   as.tibble() %>% 
#   rename(Location_data = 1) %>% 
#   write_csv('UnMethylation_location_data.csv')


## read new data format to make sure this actually worked

# methy = read_csv('Methylated_cov_nozero.csv')
# Loc_data = read_csv('Methylation_location_data.csv')
# methy_data = bind_cols(Loc_data, 
#                        methy)
# 
# methy_data %>% 
#   rotate_df() %>% 
#   rownames_to_column() %>% 
#   write_csv('Methylated_data_clean.csv', 
#             col_names = F)

# unmethy = read_csv('Unmethylated_cov_nozero.csv')
# unmethy_loc = read_csv('UnMethylation_location_data.csv')
# 
# unmethy_data = bind_cols(unmethy_loc, 
#                           unmethy)
# unmethy_data %>% 
#   rotate_df() %>% 
#   rownames_to_column() %>% 
#   write_csv('UnMethylated_data_clean.csv', 
#             col_names = F)


methy_clean = read_csv('Methylated_data_clean.csv')

unmethy_clean = read_csv('UnMethylated_data_clean.csv')

methy = methy_clean %>% 
  select(-Location_data)

unmethy = unmethy_clean %>% 
  select(-Location_data)

# methy_test = methy_clean %>% 
#   slice(1:10) %>% 
#   select(1:15) %>% 
#   select(-Location_data)
# 
# unmethy_test = unmethy_clean %>% 
#   slice(1:10) %>% 
#   select(1:15) %>% 
#   select(-Location_data)
# 
# beta_denom = map2_df(methy_test, 
#         unmethy_test, 
#         `+`) %>%  
#   map2_df(., 
#           100, 
#           `+`)
# 
# beta_test = map2_df(methy_test, 
#                     beta_denom, 
#                     `/`)

beta_denom = map2_df(methy,
                     unmethy,
                             `+`) %>%
                       map2_df(.,
                               100,
                               `+`)
beta_values = map2_df(methy,
                    beta_denom,
                    `/`)
