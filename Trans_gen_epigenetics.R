##############################
## Transgenerational epigenetics
##
## Matt Brachmann (PhDMattyB)
##
## 07.10.2024
##
##############################


setwd('C:/Users/phdma/OneDrive/OneDrive - University of Glasgow/Documents/Parsons_Postdoc/Stickleback_Genomic/Transgen_epigenetics/')

library(tidyverse)

methy = read_rds(file = "Methylated.cov.noZero.rds")

methy %>% 
  as.data.frame() %>% 
  write_csv('Methylated_cov_nozero.csv')

unmethy = read_rds(file = 'UnMethylated.cov.noZero.rds')

unmethy %>% 
  as.data.frame() %>% 
  write_csv('Unmethylated_cov_nozero.csv')

## read new data format to make sure this actually worked

methy = read_csv('Methylated_cov_nozero.csv')
