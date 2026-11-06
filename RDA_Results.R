## RDA results for TGP and WGP

library(tidyverse)

head_wgp = read_csv('~/Methylation_data/WGP_HEAD_clean_RDA_CAND_corr.csv') %>% 
  separate(col = loc, 
           into = c('Chromosome', 
                    'POS'), 
           sep = '-')


head_wgp %>%
  group_by(Association) %>% 
  summarize(n = n())

head_wgp %>% 
  filter(Chromosome == 'chrXXI') %>% 
  arrange(POS) %>% View()


head_wgp %>% 
  arrange(desc(Interaction))


head_wgp %>% 
  filter(WGP_HEAD_clean > abs(0.1))


head_tgp = read_csv('~/Methylation_data/TGP_HEAD_clean_RDA_CAND_corr.csv') %>% 
  separate(col = loc, 
           into = c('Chromosome', 
                    'POS'), 
           sep = '-')

head_tgp %>%
  group_by(Association) %>% 
  summarize(n = n())

head_tgp %>% 
  filter(Chromosome == 'chrXXI') %>% 
  arrange(POS) 