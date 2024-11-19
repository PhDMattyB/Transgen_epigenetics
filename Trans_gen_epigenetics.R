##############################
## Transgenerational epigenetics
##
## Matt Brachmann (PhDMattyB)
##
## 07.10.2024
##
##############################


# Functions ---------------------------------------------------------------

outliers = function(x,z){
  lims = mean(x) + c(-1,1)*z*sd(x)
  
  x[x<lims[1] | x>lims[2]]
}




# Start up ----------------------------------------------------------------


# setwd('C:/Users/phdma/OneDrive/OneDrive - University of Glasgow/Documents/Parsons_Postdoc/Stickleback_Genomic/Transgen_epigenetics/')
setwd('~/Methylation_data/')

library(tidyverse)
library(sjmisc)
library(vegan)
library(patchwork)

theme_set(theme_bw())

# Analysis start -------------------------------------------------------------------

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
# mval_test = mvalues %>% 
#   select(1:10) 

# meta = mvalues %>%
#   select(1)
# 
# meta %>%
#   separate(col = Location_data,
#            into = c('SampleID',
#                     'Other',
#                     'individual'),
#            sep = '-') %>%
#   separate(col = Other,
#            into = c('Pop_data',
#                     'id'),
#            sep = '_') %>%
#   unite(col = 'SampleID',
#         c('SampleID',
#           'id',
#           'individual'),
#         sep = '_') %>%
#   separate(col = Pop_data,
#            into = c('Population',
#                     'temps'),
#            sep = '(?<=[A-Za-z])(?=[0-9])') %>%
#   separate(col = temps,
#            into = c('F1_temp',
#                     'F2_temp'),
#            sep = 2, 
#            remove = F) %>%
#   mutate(Ecotype = as.factor(case_when(
#     Population == 'GTS' ~ 'Geothermal',
#     Population == 'CSWY' ~ 'Ambient',
#     Population == 'ASHNW' ~ 'Geothermal',
#     Population == 'ASHNC' ~ 'Ambient',
#     Population == 'MYVW' ~ 'Geothermal',
#     Population == 'MYVC' ~ 'Ambient',
#     Population == 'SKRW' ~ 'Geothermal',
#     Population == 'SKRC' ~ 'Ambient'))) %>%
#   write_csv('Methylation_metadata.csv')

meta_data = read_csv('Methylation_metadata.csv')

test_df = bind_cols(meta_data, 
          mval_test)


# Phenotypic traits -------------------------------------------------------


# percent methylation -----------------------------------------------------

percent_df = bind_cols(meta_data, 
                       mvalues)

df = percent_df %>% 
  select(-SampleID,
         -F1_temp, 
         -F2_temp, 
         -Ecotype) %>% 
  group_by(Population, 
           temps) %>% 
  pivot_longer(cols = starts_with('chr'),
               names_to = 'methy_loc', 
               values_to = 'Methylation') %>% 
  separate(col = methy_loc, 
           into = c('Chromosome', 
                    'BP'), 
           sep = '-') %>% 
  arrange(Chromosome, 
          BP)


# RDA ---------------------------------------------------------------------
# methy_test = mval_test %>% 
#   select(-1)

test_pheno = meta_data %>% 
  select(temps, 
         Population)
mvalues = mvalues %>% 
  select(-1)

RDA_treatment = rda(mvalues ~ temps * Population, 
               data = test_pheno, 
               scale = T)

RsquareAdj(RDA_treatment)
summary(eigenvals(RDA_treatment, 
                  model = 'constrained'))

screeplot(RDA_treatment)

signif_full = anova.cca(RDA_treatment, 
                        parallel = getOption('mc.cores'))

# signif_axis = anova.cca(RDA_treatment, 
#                         by = 'axis',
#                         parallel = getOption('mc.cores'))


vif.cca(RDA_treatment)

sum_rda = summary(RDA_treatment)

sum_rda$species %>% 
  as_tibble() %>% 
  write_csv('RDA_treatment_pops_methy_locations.csv')

sum_rda$sites %>% 
  as_tibble() %>% 
  write_csv('RDA_treatment_pops_individuals.csv')

sum_rda$biplot %>% 
  as_tibble() %>% 
  write_csv('RDA_treatment_pops_biplot.csv')


rda_scores = scores(RDA_treatment, 
                    choices = c(1:2), 
                    display = 'species')

hist(rda_scores[,1])
hist(rda_scores[,2])

rda_outliers = outliers(rda_scores[,1], 3)
# rda_outliers_axis2 = outliers(rda_scores[,2], 3)


rda_out = cbind.data.frame(rep(1, 
                               times = length(rda_outliers)), 
                           names(rda_outliers), 
                           unname(rda_outliers))

rda_out = rda_out %>% 
  as_tibble() %>%  
dplyr::rename(axis = 1, 
              loc = 2, 
              scores = 3)


all_loc = rda_scores[,1]

rda_normal = cbind.data.frame(rep(1, 
                                  times = length(all_loc)), 
                              names(all_loc),
                              unname(all_loc))
rda_normal = rda_normal %>% 
  as_tibble() %>% 
  dplyr::rename(axis = 1, 
                loc = 2, 
                scores = 3)

rda_normal = rda_normal[!rda_normal$loc %in% rda_out$loc,]

# write_csv(rda_normal, 
#           'RDA_nonoutliers_methylation.csv')
# 
# write_csv(rda_out, 
#           'RDA_outliers_methylation.csv')

rda_out = as.data.frame(rda_out)
all_loc = as.data.frame(all_loc)
# test_pheno = as.data.frame(test_pheno)

test_pheno2 = test_pheno %>% 
  as_tibble() %>% 
mutate(pop_num = as.numeric(case_when(
      Population == 'GTS' ~ '1',
      Population == 'CSWY' ~ '2',
      Population == 'ASHNW' ~ '3',
      Population == 'ASHNC' ~ '4',
      Population == 'MYVW' ~ '5',
      Population == 'MYVC' ~ '6',
      Population == 'SKRW' ~ '7',
      Population == 'SKRC' ~ '8'))) %>% 
  # dplyr::select(-Population) %>% 
  as.data.frame()

# nam = rda_out[1:45, 2]
# out_loc = all_loc[nam,]
# out_cor = apply(test_pheno,
#                 2, 
#                 function(x)cor(x, out_loc))

foo = matrix(nrow=(45), 
             ncol = 2)
colnames(foo) = c('temps', 
                  'pop_num')

for (i in 1:length(rda_out$loc)){
  nam = rda_out[i,2]
  loc.gen = mvalues[,nam]
  foo[i,] = apply(test_pheno,2,function(x)cor(x,loc.gen))
}

candidates = cbind.data.frame(rda_out, 
                              foo)

# candidates %>% 
#   as_tibble() %>% 
#   write_csv('RDA_outliers_methylation_correlations.csv')
            
##check for duplicates
length(candidates$loc[duplicated(candidates$loc)])

for(i in 1:length(candidates$loc)){
  bar = candidates[i,]
  candidates[i,6] = names(which.max(abs(bar[4:5])))
  candidates[i,7] = max(abs(bar[4:5]))
  }

candidates



# Hyper and hypo methylation candidate loci -------------------------------


cand_loc = candidates$loc

cand_methy = mvalues %>% 
  select(any_of(cand_loc))

cand_methy = bind_cols(meta_data, 
                       cand_methy)

cand_methy %>% 
  arrange(Population,
          temps) %>% 
  View()

cand_methy_pivot = cand_methy %>% 
  group_by(SampleID, 
           Population, 
           temps, 
           F1_temp, 
           F2_temp, 
           Ecotype) %>% 
  pivot_longer(cols = starts_with('chr'),
               names_to = 'methy_loc', 
               values_to = 'Methylation') %>% 
  separate(col = methy_loc, 
           into = c('Chromosome', 
                    'BP'), 
           sep = '-') %>% 
  arrange(Chromosome, 
          BP)

cand_methy_pivot$temps = as.character(cand_methy_pivot$temps)


# cand_methy_pivot %>% 
#   filter(Population %in% c('ASHNC', 
#                            'ASHNW')) %>% 
#   ggplot()+
#   geom_point(aes(x = BP, 
#                  y = Methylation, 
#                  col = temps))+
#   facet_grid(~Chromosome)

# cand_methy_pivot %>%
#   filter(Population == 'ASHNC') %>%
#   group_by(Population,
#            temps) %>%
#   distinct(Chromosome,
#            BP,
#            .keep_all = T) %>%
#   arrange(Chromosome,
#           BP,
#           temps) %>%
#   View()

F1_temps_pal = c('#0077b6', 
              '#a2d2ff', 
              '#ef959c',
              '#ef233c')

cand_methy_pivot$Chromosome = factor(cand_methy_pivot$Chromosome, 
                                     levels = c('chrIII', 
                                                'chrIV', 
                                                'chrVI', 
                                                'chrXI', 
                                                'chrXII', 
                                                'chrXIII', 
                                                'chrXVI', 
                                                'chrXX', 
                                                'chrUn'))


ASHNC_outlier_plot = cand_methy_pivot %>% 
  filter(Population == 'ASHNC') %>% 
  group_by(Population, 
           temps) %>% 
  distinct(Chromosome, 
           BP, 
           .keep_all = T) %>% 
  ggplot()+
  # geom_point(aes(x = BP, 
  #                y = Methylation, 
  #                col = temps))+
  geom_jitter(aes(x = BP, 
                 y = Methylation, 
                 fill = temps), 
              col = 'black',
              pch = 21,
              size = 2,
              width = 0, 
              height = 0.05)+
  scale_y_continuous(expand = c(0,0), 
                     limits = c(-6.0, 1.0), 
                     breaks = c(-6.0, 
                                -5.0, 
                                -4.0, 
                                -3.0, 
                                -2.0, 
                                -1.0, 
                                0.0, 
                                1.0))+
  geom_hline(yintercept = 0.0)+
  labs(title = 'ASHN Cold')+
  scale_fill_manual(values = F1_temps_pal)+
  facet_grid(~Chromosome)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        # axis.title.y = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(size = 12, 
                                  face = 'bold'),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(), 
        legend.title = element_blank(), 
        legend.position = 'none')


ASHNW_outlier_plot = cand_methy_pivot %>% 
  filter(Population == 'ASHNW') %>% 
  group_by(Population, 
           temps) %>% 
  distinct(Chromosome, 
           BP, 
           .keep_all = T) %>% 
  ggplot()+
  # geom_point(aes(x = BP, 
  #                y = Methylation, 
  #                col = temps))+
  geom_jitter(aes(x = BP, 
                 y = Methylation, 
                 fill = temps),
              col = 'black', 
              pch = 21,
              size = 2,
              width = 0, 
              height = 0.05)+
  scale_y_continuous(expand = c(0,0), 
                     limits = c(-6.0, 1.0), 
                     breaks = c(-6.0, 
                                -5.0, 
                                -4.0, 
                                -3.0, 
                                -2.0, 
                                -1.0, 
                                0.0, 
                                1.0))+
  geom_hline(yintercept = 0.0)+
  labs(title = 'ASHN Warm', 
       x = 'Base pair position')+
  scale_fill_manual(values = F1_temps_pal)+
  facet_grid(~Chromosome)+
  theme(axis.title.y = element_blank(), 
        # axis.title.y = element_text(size = 14), 
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(size = 12, 
                                  face = 'bold'),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(), 
        legend.title = element_blank(), 
        legend.position = 'none')


ASHN_meth_outliers = ASHNC_outlier_plot/ASHNW_outlier_plot


MYVC_outlier_plot = cand_methy_pivot %>% 
  filter(Population == 'MYVC') %>% 
  group_by(Population, 
           temps) %>% 
  distinct(Chromosome, 
           BP, 
           .keep_all = T) %>% 
  ggplot()+
  # geom_point(aes(x = BP, 
  #                y = Methylation, 
  #                col = temps))+
  geom_jitter(aes(x = BP, 
                  y = Methylation, 
                  fill = temps),
              col = 'black', 
              pch = 21,
              size = 2,
              width = 0, 
              height = 0.05)+
  scale_y_continuous(expand = c(0,0), 
                     limits = c(-6.0, 1.0), 
                     breaks = c(-6.0, 
                                -5.0, 
                                -4.0, 
                                -3.0, 
                                -2.0, 
                                -1.0, 
                                0.0, 
                                1.0))+
  geom_hline(yintercept = 0.0)+
  labs(title = 'MYV Cold', 
       x = 'Base pair position')+
  scale_fill_manual(values = F1_temps_pal)+
  facet_grid(~Chromosome)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        # axis.title.y = element_text(size = 14), 
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(size = 12, 
                                  face = 'bold'),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(), 
        legend.title = element_blank(), 
        legend.position = 'none')

MYVW_outlier_plot = cand_methy_pivot %>% 
  filter(Population == 'MYVW') %>% 
  group_by(Population, 
           temps) %>% 
  distinct(Chromosome, 
           BP, 
           .keep_all = T) %>% 
  ggplot()+
  # geom_point(aes(x = BP, 
  #                y = Methylation, 
  #                col = temps))+
  geom_jitter(aes(x = BP, 
                  y = Methylation, 
                  fill = temps),
              col = 'black', 
              pch = 21,
              size = 2,
              width = 0, 
              height = 0.05)+
  scale_y_continuous(expand = c(0,0), 
                     limits = c(-6.0, 1.0), 
                     breaks = c(-6.0, 
                                -5.0, 
                                -4.0, 
                                -3.0, 
                                -2.0, 
                                -1.0, 
                                0.0, 
                                1.0))+
  geom_hline(yintercept = 0.0)+
  labs(title = 'MYV Warm', 
       x = 'Base pair position')+
  scale_fill_manual(values = F1_temps_pal)+
  facet_grid(~Chromosome)+
  theme(axis.title = element_blank(), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(size = 12, 
                                  face = 'bold'),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(), 
        legend.title = element_blank(), 
        legend.position = 'none')

MYV_meth_outliers = MYVC_outlier_plot/MYVW_outlier_plot

SKRC_outlier_plot = cand_methy_pivot %>% 
  filter(Population == 'SKRC') %>% 
  group_by(Population, 
           temps) %>% 
  distinct(Chromosome, 
           BP, 
           .keep_all = T) %>% 
  ggplot()+
  # geom_point(aes(x = BP, 
  #                y = Methylation, 
  #                col = temps))+
  geom_jitter(aes(x = BP, 
                  y = Methylation, 
                  fill = temps),
              col = 'black', 
              pch = 21,
              size = 2,
              width = 0, 
              height = 0.05)+
  scale_y_continuous(expand = c(0,0), 
                     limits = c(-6.0, 1.0), 
                     breaks = c(-6.0, 
                                -5.0, 
                                -4.0, 
                                -3.0, 
                                -2.0, 
                                -1.0, 
                                0.0, 
                                1.0))+
  geom_hline(yintercept = 0.0)+
  labs(title = 'SKR Cold', 
       x = 'Base pair position')+
  scale_fill_manual(values = F1_temps_pal)+
  facet_grid(~Chromosome)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        # axis.title.y = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(size = 12, 
                                  face = 'bold'),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(), 
        legend.title = element_blank(), 
        legend.position = 'none')
  

SKRW_outlier_plot = cand_methy_pivot %>% 
  filter(Population == 'SKRW') %>% 
  group_by(Population, 
           temps) %>% 
  distinct(Chromosome, 
           BP, 
           .keep_all = T) %>% 
  ggplot()+
  # geom_point(aes(x = BP, 
  #                y = Methylation, 
  #                col = temps))+
  geom_jitter(aes(x = BP, 
                  y = Methylation, 
                  fill = temps),
              col = 'black', 
              pch = 21, 
              size = 2,
              width = 0, 
              height = 0.05)+
  scale_y_continuous(expand = c(0,0), 
                     limits = c(-6.0, 1.0), 
                     breaks = c(-6.0, 
                                -5.0, 
                                -4.0, 
                                -3.0, 
                                -2.0, 
                                -1.0, 
                                0.0, 
                                1.0))+
  geom_hline(yintercept = 0.0)+
  labs(title = 'SKR Warm', 
       x = 'Base pair position')+
  scale_fill_manual(values = F1_temps_pal)+
  facet_grid(~Chromosome)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(size = 12, 
                                  face = 'bold'),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(), 
        legend.title = element_blank(), 
        legend.position = 'none')

SKR_meth_outlier = SKRC_outlier_plot/SKRW_outlier_plot


CSWY_outlier_plot = cand_methy_pivot %>% 
  filter(Population == 'CSWY') %>% 
  group_by(Population, 
           temps) %>% 
  distinct(Chromosome, 
           BP, 
           .keep_all = T) %>% 
  ggplot()+
  # geom_point(aes(x = BP, 
  #                y = Methylation, 
  #                col = temps))+
  geom_jitter(aes(x = BP, 
                  y = Methylation, 
                  fill = temps),
              col = 'black', 
              pch = 21,
              size = 2,
              width = 0, 
              height = 0.05)+
  scale_y_continuous(expand = c(0,0), 
                     limits = c(-6.0, 1.0), 
                     breaks = c(-6.0, 
                                -5.0, 
                                -4.0, 
                                -3.0, 
                                -2.0, 
                                -1.0, 
                                0.0, 
                                1.0))+
  geom_hline(yintercept = 0.0)+
  labs(title = 'CSWY Cold', 
       x = 'Base pair position')+
  scale_fill_manual(values = F1_temps_pal)+
  facet_grid(~Chromosome)+
  theme(axis.title.x = element_text(size = 14), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        # axis.title.y = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(size = 12, 
                                  face = 'bold'),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(), 
        legend.title = element_blank(), 
        legend.position = 'bottom')

GTS_outlier_plot = cand_methy_pivot %>% 
  filter(Population == 'GTS') %>% 
  group_by(Population, 
           temps) %>% 
  distinct(Chromosome, 
           BP, 
           .keep_all = T) %>% 
  ggplot()+
  # geom_point(aes(x = BP, 
  #                y = Methylation, 
  #                col = temps))+
  geom_jitter(aes(x = BP, 
                  y = Methylation, 
                  fill = temps),
              col = 'black', 
              pch = 21,
              size = 2,
              width = 0, 
              height = 0.05)+
  scale_y_continuous(expand = c(0,0), 
                     limits = c(-6.0, 1.0), 
                     breaks = c(-6.0, 
                                -5.0, 
                                -4.0, 
                                -3.0, 
                                -2.0, 
                                -1.0, 
                                0.0, 
                                1.0))+
  geom_hline(yintercept = 0.0)+
  labs(title = 'GTS Warm', 
       x = 'Base pair position')+
  scale_fill_manual(values = F1_temps_pal)+
  facet_grid(~Chromosome)+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(size = 12, 
                                  face = 'bold'),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(), 
        legend.title = element_blank(), 
        legend.position = 'bottom')

CSWY_GTS_meth_outlier = CSWY_outlier_plot/GTS_outlier_plot


BIG_METHY = (ASHNC_outlier_plot|ASHNW_outlier_plot)/(MYVC_outlier_plot|MYVW_outlier_plot)/(SKRC_outlier_plot|SKRW_outlier_plot)/(CSWY_outlier_plot|GTS_outlier_plot)

ggsave('Big_methylaton_treatment_plot.tiff', 
       plot = BIG_METHY, 
       dpi = 'retina', 
       units = 'cm', 
       height = 25, 
       width = 40)

# GRAPHS! -----------------------------------------------------------------


candidates = read_csv('RDA_outliers_methylation_correlations.csv')
normal_loc = read_csv('RDA_nonoutliers_methylation.csv')
loc_data = read_csv('RDA_treatment_pops_methy_locations.csv')
individuals = read_csv('RDA_treatment_pops_individuals.csv')
biplot = read_csv('RDA_treatment_pops_biplot.csv')

outlier_label = rep('Outlier', 
                    length(candidates$loc)) %>% 
  as_tibble() %>% 
  rename(label = value)

candidates = bind_cols(candidates, 
                       outlier_label)

normal_label = rep('Normal', 
                   length(normal_loc$loc)) %>% 
  as_tibble() %>% 
  rename(label = value)

normal_loc = bind_cols(normal_loc, 
                       normal_label)


vars_rda = bind_cols(test_pheno, 
          individuals) 

vars_rda$temps = as.character(vars_rda$temps)

# all_scores = rda_scores %>%
#   as.data.frame() %>%
#   rownames_to_column() %>%
#   as_tibble() %>%
#   rename(loc = rowname)
# 
# inner_join(candidates,
#            all_scores,
#            by = 'loc') %>%
#   write_csv('RDA_outliers_methylation_correlations.csv')


theme_set(theme_bw())

ggplot()+
  geom_point(data = loc_data,
             aes(x = RDA1,
                 y = RDA2),
             col = '#ADA597',
             size = 2)+
  geom_point(data = candidates,
             aes(x = RDA1,
                 y = RDA2),
             col = '#0AB33A',
             size = 2)+
  labs(x = 'RDA 1 (0.65% variance explained)',
       y = 'RDA 2 (0.13% variance explained)')+
  theme(#legend.position = "none", 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.title = element_text(size = 15), 
    axis.text = element_text(size = 15), 
    # axis.ticks = element_line(size = 1), 
    plot.title = element_text(size = 15, 
                              hjust = 0), 
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12))
  # geom_point(data = vars_rda, 
  #            aes(x = RDA1, 
  #                y = RDA2, 
  #                col = temps), 
  #            size = 2)+
  # geom_segment(aes(xend = RDA_treatment$CCA$biplot[,1], 
  #                  yend = RDA_treatment$CCA$biplot[,2], 
  #                  x = 0, 
  #                  y = 0), 
  #              colour = 'black', 
  #              size = 1, 
  #              linetype = 1, 
  #              arrow = arrow(length = unit(0.1, 
  #                                          'npc')))+
  # geom_text(aes(x = 1.5*RDA_treatment$CCA$biplot[,1], 
  #               y = 1.2*RDA_treatment$CCA$biplot[,2], 
  #               label = colnames(vars_rda[,1:2])))



# RAW whole body RDA ------------------------------------------------------

mvalues = read_csv('MVALUES_methylation_cleaned_data.csv')
meta_data = read_csv('Methylation_metadata.csv')
raw_data = read_csv('Raw_RDA_phenotypes.csv')


meth_fish = mvalues %>% 
  select(1)

meth_fish_ID = meth_fish %>% 
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

pheno_fish = raw_data %>% 
  filter(str_detect(fish,
                    '_#G')) %>% 
  distinct(fish, 
           .keep_all = T) %>% 
  select(fish)

pheno_fish_ID = pheno_fish %>% 
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

pheno_fish_ID$var1 = gsub("'", '', pheno_fish_ID$var1)

pheno_fish_ID = pheno_fish_ID %>% 
  unite(col = 'Fish_ID', 
        sep = '')%>% 
  mutate(Fish_ID = gsub("Myvat", 
                        "MYV", 
                        Fish_ID)) 


pheno_fish = raw_data %>% 
  filter(str_detect(fish,
                    '_#G')) %>% 
  distinct(fish, 
           .keep_all = T) %>% 
  # select(Full_ID) %>% 
  bind_cols(pheno_fish_ID, 
                       .)
## identifying potential issues with phenotypic data
## 5 individuals that were sequenced did not have the #G identifier
## Four from GTS18@12 and one from MYVC12@12
anti_join(meth_fish_ID,
          pheno_fish_ID)


pheno_fish_final = inner_join(meth_fish_ID, 
          pheno_fish) %>% 
  arrange(Fish_ID)

## now we need to order the phenotypic data and methylation data
## they have to be in the same order otherwise this will all
## be fucked. Right now they aren't. 
## order the data by the Fish_ID column

## Add the new id column to the methylation data

mvalues_final = bind_cols(meth_fish_ID, 
          mvalues) %>% 
arrange(Fish_ID)  

## need to check that everythings in order for the analyses
## If there are any FALSE we're fucked. 
## Shooting for all TRUES
mvalues_final$Fish_ID == pheno_fish_final$Fish_ID

mvalues = mvalues %>% 
  select(-1)

RDA_treatment = rda(mvalues ~ Comp1 + Comp2 + Comp3 + Comp4 + Comp5, 
                    data = pheno_fish_final, 
                    scale = T)

RsquareAdj(RDA_treatment)
summary(eigenvals(RDA_treatment, 
                  model = 'constrained'))

screeplot(RDA_treatment)

signif_full = anova.cca(RDA_treatment, 
                        parallel = getOption('mc.cores'))

# signif_axis = anova.cca(RDA_treatment, 
#                         by = 'axis',
#                         parallel = getOption('mc.cores'))


vif.cca(RDA_treatment)

sum_rda = summary(RDA_treatment)

sum_rda$species %>% 
  as_tibble() %>% 
  write_csv('RAW_Uncorrected_PCA_locations.csv')

sum_rda$sites %>% 
  as_tibble() %>% 
  write_csv('RDA_Uncorrected_PCA_individuals.csv')

sum_rda$biplot %>% 
  as_tibble() %>% 
  write_csv('RDA_Uncorrected_PCA_biplot.csv')


rda_scores = scores(RDA_treatment, 
                    choices = c(1:5), 
                    display = 'species')

hist(rda_scores[,1])
hist(rda_scores[,2])
hist(rda_scores[,3])
hist(rda_scores[,4])
hist(rda_scores[,5])

rda_outliers_axis1 = outliers(rda_scores[,1], 3)
# rda_outliers_axis2 = outliers(rda_scores[,2], 3)


rda_out = cbind.data.frame(rep(1, 
                               times = length(rda_outliers)), 
                           names(rda_outliers), 
                           unname(rda_outliers))

rda_out = rda_out %>% 
  as_tibble() %>%  
  dplyr::rename(axis = 1, 
                loc = 2, 
                scores = 3)


all_loc = rda_scores[,1]

rda_normal = cbind.data.frame(rep(1, 
                                  times = length(all_loc)), 
                              names(all_loc),
                              unname(all_loc))
rda_normal = rda_normal %>% 
  as_tibble() %>% 
  dplyr::rename(axis = 1, 
                loc = 2, 
                scores = 3)

rda_normal = rda_normal[!rda_normal$loc %in% rda_out$loc,]

# write_csv(rda_normal, 
#           'RDA_nonoutliers_methylation.csv')
# 
# write_csv(rda_out, 
#           'RDA_outliers_methylation.csv')

rda_out = as.data.frame(rda_out)
all_loc = as.data.frame(all_loc)
# test_pheno = as.data.frame(test_pheno)

test_pheno2 = test_pheno %>% 
  as_tibble() %>% 
  mutate(pop_num = as.numeric(case_when(
    Population == 'GTS' ~ '1',
    Population == 'CSWY' ~ '2',
    Population == 'ASHNW' ~ '3',
    Population == 'ASHNC' ~ '4',
    Population == 'MYVW' ~ '5',
    Population == 'MYVC' ~ '6',
    Population == 'SKRW' ~ '7',
    Population == 'SKRC' ~ '8'))) %>% 
  # dplyr::select(-Population) %>% 
  as.data.frame()

# nam = rda_out[1:45, 2]
# out_loc = all_loc[nam,]
# out_cor = apply(test_pheno,
#                 2, 
#                 function(x)cor(x, out_loc))

foo = matrix(nrow=(45), 
             ncol = 2)
colnames(foo) = c('temps', 
                  'pop_num')

for (i in 1:length(rda_out$loc)){
  nam = rda_out[i,2]
  loc.gen = mvalues[,nam]
  foo[i,] = apply(test_pheno,2,function(x)cor(x,loc.gen))
}

candidates = cbind.data.frame(rda_out, 
                              foo)

# candidates %>% 
#   as_tibble() %>% 
#   write_csv('RDA_outliers_methylation_correlations.csv')

##check for duplicates
length(candidates$loc[duplicated(candidates$loc)])

for(i in 1:length(candidates$loc)){
  bar = candidates[i,]
  candidates[i,6] = names(which.max(abs(bar[4:5])))
  candidates[i,7] = max(abs(bar[4:5]))
}

candidates


