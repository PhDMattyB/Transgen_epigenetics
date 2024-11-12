


methy_rds = read_rds(file = 'edgeR.data.gt10.Coverage.noZero.rds')
head(methy_rds)


chrI_data = methy_rds %>%
  # as.matrix() %>%
  # rownames_to_column(var = 'Location') %>%
  as.data.frame() %>%
  as_tibble() %>% 
  select(contains('SKR'))

chrI_names = row.names(methy_rds) %>%
  as.data.frame() %>%
  as.tibble() %>%
  rename(Location_data = 1) 

bind_cols(chrI_names, 
                 chrI_data) %>% 
  filter(str_detect(Location_data, 'chrI-')) %>% 
  write_csv('ChrI_data_Corin.csv')


test = read_csv("ChrI_data_Corin.csv")
