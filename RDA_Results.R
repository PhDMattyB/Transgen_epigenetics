## RDA results for TGP and WGP

library(tidyverse)

head_wgp = read_csv('~/Methylation_data/WGP_HEAD_clean_RDA_CAND_corr.csv') %>% 
  separate(col = loc, 
           into = c('Chromosome', 
                    'POS'), 
           sep = '-')


# head_wgp %>% 
#   arrange(desc(Correlation)) %>% View()

head_wgp %>% 
  filter(Chromosome == 'chrXXI') %>% 
  arrange(POS) %>% View()
