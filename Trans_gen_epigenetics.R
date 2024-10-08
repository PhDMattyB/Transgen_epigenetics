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
library(devtools)
install.packages('freeR')
install_github('cfree14/freeR')

library(freeR)

methy_rds = read_rds(file = "Methylated.cov.noZero.rds")
head(methy_rds)

methy_rds %>% 
  # as.matrix() %>% 
  # rownames_to_column(var = 'Location') %>% 
  as.data.frame() %>%
  as_tibble() %>% 
  # rowid_to_column(var = 'rowid')
  write_csv('Methylated_cov_nozero.csv')


row.names(methy_rds) %>% 
  as.data.frame() %>% 
  as.tibble() %>% 
  rename(Location_data = 1) %>% 
  write_csv('Methylation_location_data.csv')

  
unmethy_rds = read_rds(file = 'UnMethylated.cov.noZero.rds')

unmethy %>% 
  as.data.frame() %>% 
  write_csv('Unmethylated_cov_nozero.csv')

row.names(unmethy_rds) %>% 
  as.data.frame() %>% 
  as.tibble() %>% 
  rename(Location_data = 1) %>% 
  write_csv('UnMethylation_location_data.csv')


## read new data format to make sure this actually worked

methy = read_csv('Methylated_cov_nozero.csv')

Loc_data = read_csv('Methylation_location_data.csv')

methy_data = bind_cols(Loc_data, 
                       methy)

