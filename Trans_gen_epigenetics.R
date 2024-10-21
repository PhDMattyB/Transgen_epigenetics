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

methy = read_csv('Methylated_cov_nozero.csv')
Loc_data = read_csv('Methylation_location_data.csv')
methy_data = bind_cols(Loc_data,
                       methy)
# 
# methy_data %>% 
#   rotate_df() %>% 
#   rownames_to_column() %>% 
#   write_csv('Methylated_data_clean.csv', 
#             col_names = F)

unmethy = read_csv('Unmethylated_cov_nozero.csv')
unmethy_loc = read_csv('UnMethylation_location_data.csv')
# 
unmethy_data = bind_cols(unmethy_loc,
                          unmethy)
# unmethy_data %>% 
#   rotate_df() %>% 
#   rownames_to_column() %>% 
#   write_csv('UnMethylated_data_clean.csv', 
#             col_names = F)

methy = read_csv('Methylated_data_clean.csv', 
                 col_names = F)

unmethy = read_csv('UnMethylated_data_clean.csv', 
                   col_names = F)
 

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


# calculate beta values ---------------------------------------------------


beta_denom = map2_df(methy,
                     unmethy,
                             `+`) %>%
                       map2_df(.,
                               100,
                               `+`)
beta_values = map2_df(methy,
                    beta_denom,
                    `/`)
# 
# beta_values %>% 
#   slice(1:10) %>% 
#   select(1:10)

meta = methy_clean %>% 
  select(Location_data)


full_data = bind_cols(meta, 
                      beta_values)

sum(is.na(full_data))




# calculate m-values ------------------------------------------------------

Methy_m = map2_df(methy,
        100,
        `+`) 
Unmethy_m = map2_df(unmethy, 
                    100, 
                    `+`)

divde = map2_df(Methy_m, 
                Unmethy_m, 
                `/`)

# test_df = divde %>% 
#   slice(1:10) %>% 
#   select(1:10)
# 
# test_df %>% 
#   mutate(across(everything(),
#                 ~log2(.)))

mvalues = divde %>% 
  mutate(across(everything(),
                ~log2(.)))

## checks to see if the data is the way it should be
sum(is.na(mvalues))
mvalues %>% 
  purrr::keep(~any(. < 0))

## combine with individual id data
mvalues_final = bind_cols(meta, 
                          mvalues)  

## write a copy so we don't have to do that again
mvalues_final %>% 
  write_csv('MVALUES_methylation_cleaned_data.csv')
# mvalues = log2(divide)

## read in mvalue data
mvalues = read_csv('MVALUES_methylation_cleaned_data.csv')

## test set for down stream data wrangling
mval_test = mvalues %>% 
  select(1:10) 

meta = mvalues %>%
  select(1)

meta %>%
  separate(col = Location_data,
           into = c('SampleID',
                    'Other',
                    'individual'),
           sep = '-') %>%
  separate(col = Other,
           into = c('Pop_data',
                    'id'),
           sep = '_') %>%
  unite(col = 'SampleID',
        c('SampleID',
          'id',
          'individual'),
        sep = '_') %>%
  separate(col = Pop_data,
           into = c('Population',
                    'temps'),
           sep = '(?<=[A-Za-z])(?=[0-9])') %>%
  separate(col = temps,
           into = c('F1_temp',
                    'F2_temp'),
           sep = 2, 
           remove = F) %>%
  mutate(Ecotype = as.factor(case_when(
    Population == 'GTS' ~ 'Geothermal',
    Population == 'CSWY' ~ 'Ambient',
    Population == 'ASHNW' ~ 'Geothermal',
    Population == 'ASHNC' ~ 'Ambient',
    Population == 'MYVW' ~ 'Geothermal',
    Population == 'MYVC' ~ 'Ambient',
    Population == 'SKRW' ~ 'Geothermal',
    Population == 'SKRC' ~ 'Ambient'))) %>%
  write_csv('Methylation_metadata.csv')

meta_data = read_csv('Methylation_metadata.csv')

test_df = bind_cols(meta_data, 
          mval_test)


# Phenotypic traits -------------------------------------------------------



# RDA ---------------------------------------------------------------------
methy_test = mval_test %>% 
  select(-1)

test_pheno = meta_data %>% 
  select(temps, 
         Population)
mvalues = mvalues %>% 
  select(-1)

library(vegan)

test_rda = rda(mvalues ~ temps * Population, 
               data = test_pheno, 
               scale = T)

RsquareAdj(test_rda)
summary(eigenvals(test_rda, 
                  model = 'constrained'))

screeplot(test_rda)

signif_full = anova.cca(test_rda, 
                        parallel = )

