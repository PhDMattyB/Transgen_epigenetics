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

mval_small %>% 
  select(Location_data)
