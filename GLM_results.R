setwd('~/Parsons_Postdoc/Methylation_data/GLM_results/')

library(tidyverse)
library(qvalue)
library(windowscanr)
library(future)
library(future.apply)
library(patchwork)

mod_results = read_table('model_results_table.tab')


F2_res = mod_results %>% 
  dplyr::select(Methy_loc,
         Chromosome, 
         Pos,
         expression,
         estimate_F1,
         pval_F1) %>% 
  filter(!pval_F1 == 'NaN') %>% 
  filter(pval_F1 <= 0.001)


F2_res = mod_results %>% 
  dplyr::select(Methy_loc, 
         Chromosome, 
         Pos,
         expression,
         estimate_F2,
         pval_F2) %>% 
  filter(!pval_F2 == 'NaN') %>% 
  filter(pval_F2 <= 0.001)


eco_res = mod_results %>% 
  dplyr::select(Methy_loc, 
         Chromosome, 
         Pos,
         expression, 
         estimate_ecotype,
         pval_eco) %>% 
  filter(!pval_eco == 'NaN') %>% 
  filter(pval_eco <= 0.001)

F2_eco_res = mod_results %>% 
  dplyr::select(Methy_loc, 
         Chromosome, 
         Pos,
         expression,
         estimate_F2_eco,
         pval_F2_eco) %>% 
  filter(!pval_F2_eco == 'NaN') %>% 
  filter(pval_F2_eco <= 0.001)

F2_eco_res = mod_results %>% 
  dplyr::select(Methy_loc, 
         Chromosome, 
         Pos,
         expression,
         estimate_F2_eco,
         pval_F2_eco) %>% 
  filter(!pval_F2_eco == 'NaN') %>% 
  filter(pval_F2_eco <= 0.001)

F2_F2_res = mod_results %>% 
  dplyr::select(Methy_loc, 
         Chromosome, 
         Pos,
         expression,
         estimate_F2_F2,
         pval_F2_F2) %>% 
  filter(!pval_F2_F2 == 'NaN') %>% 
  filter(pval_F2_F2 <= 0.001)

F2_F2_eco_res = mod_results %>% 
  dplyr::select(Methy_loc,
         Chromosome, 
         Pos,
         expression,
         estimate_F2_F2_eco,
         pval_F2_F2_eco) %>% 
  filter(!pval_F2_F2_eco == 'NaN') %>% 
  filter(pval_F2_F2_eco <= 0.001)



inner_join(F2_res, 
           F2_res, 
           by = c('Methy_loc'))




# sliding window ----------------------------------------------------------


F2_res$Pos = as.integer(F2_res$Pos)

F2_test = F2_res %>% 
  as.data.frame()

plan(multisession, workers = 8)

win_start = seq(1, 1, by = 100)

meth_F2_win = future_lapply(
  win_start,
  function(x){
    winScan(x = F2_res, 
            groups = 'Chromosome', 
            position = 'Pos', 
            values = c('expression'), 
            win_size = 1000, 
            win_step = 1000, 
            funs = c('mean', 'sd'))
  }
)

meth_F2_win[[1]] %>%
  as_tibble() %>% 
  filter(expression_n >= 2) %>% 
  write_tsv('Methylation_F2_effect_win1000_step1000.txt')


win_start = seq(1, 1, by = 100)

meth_F2_eco_win = future_lapply(
  win_start,
  function(x){
    winScan(x = F2_eco_res, 
            groups = 'Chromosome', 
            position = 'Pos', 
            values = c('expression'), 
            win_size = 1000, 
            win_step = 1000, 
            funs = c('mean', 'sd'))
  }
)

meth_F2_eco_win[[1]] %>%
  as_tibble() %>% 
  filter(expression_n >= 2) %>% 
  write_tsv('Methylation_F2_eco_effect_win1000_step1000.txt')



win_start = seq(1, 1, by = 100)

meth_f2_win = future_lapply(
  win_start,
  function(x){
    winScan(x = F2_res, 
            groups = 'Chromosome', 
            position = 'Pos', 
            values = c('expression'), 
            win_size = 1000, 
            win_step = 1000, 
            funs = c('mean', 'sd'))
  }
)

meth_f2_win[[1]] %>%
  as_tibble() %>% 
  filter(expression_n >= 2) %>% 
  write_tsv('Methylation_F2_effect_win1000_step1000.txt')


win_start = seq(1, 1, by = 100)

meth_f2_eco_win = future_lapply(
  win_start,
  function(x){
    winScan(x = F2_eco_res, 
            groups = 'Chromosome', 
            position = 'Pos', 
            values = c('expression'), 
            win_size = 1000, 
            win_step = 1000, 
            funs = c('mean', 'sd'))
  }
)

meth_f2_eco_win[[1]] %>%
  as_tibble() %>% 
  filter(expression_n >= 2) %>% 
  write_tsv('Methylation_F2_eco_effect_win1000_step1000.txt')


F2_F2_res

win_start = seq(1, 1, by = 100)

meth_F2_f2_win = future_lapply(
  win_start,
  function(x){
    winScan(x = F2_F2_res, 
            groups = 'Chromosome', 
            position = 'Pos', 
            values = c('expression'), 
            win_size = 1000, 
            win_step = 1000, 
            funs = c('mean', 'sd'))
  }
)

meth_F2_f2_win[[1]] %>%
  as_tibble() %>% 
  filter(expression_n >= 2) %>% 
  write_tsv('Methylation_F2_F2_effect_win1000_step1000.txt')

F2_F2_eco_res
win_start = seq(1, 1, by = 100)

meth_F2_f2_eco_win = future_lapply(
  win_start,
  function(x){
    winScan(x = F2_F2_eco_res, 
            groups = 'Chromosome', 
            position = 'Pos', 
            values = c('expression'), 
            win_size = 1000, 
            win_step = 1000, 
            funs = c('mean', 'sd'))
  }
)

meth_F2_f2_eco_win[[1]] %>%
  as_tibble() %>% 
  filter(expression_n >= 2) %>% 
  write_tsv('Methylation_F2_F2_eco_effect_win1000_step1000.txt')



eco_res
win_start = seq(1, 1, by = 100)

meth_eco_win = future_lapply(
  win_start,
  function(x){
    winScan(x = eco_res, 
            groups = 'Chromosome', 
            position = 'Pos', 
            values = c('expression'), 
            win_size = 1000, 
            win_step = 1000, 
            funs = c('mean', 'sd'))
  }
)

meth_eco_win[[1]] %>%
  as_tibble() %>% 
  filter(expression_n >= 2) %>% 
  write_tsv('Methylation_eco_effect_win1000_step1000.txt')


# sliding window neutral locations ----------------------------------------

mod_results = read_table('model_results_table.tab') %>% 
  filter(!pval_F1 == 'NaN') %>% 
  as.data.frame()

options(future.globals.maxSize= 2000 * 1024^2)

win_start = seq(1, 1, by = 100)

mod_neutral_win = future_lapply(
  win_start,
  function(x){
    winScan(x = mod_results, 
            groups = 'Chromosome', 
            position = 'Pos', 
            values = c('expression'), 
            win_size = 1000, 
            win_step = 1000, 
            funs = c('mean', 'sd'))
  }
)

mod_neutral_win[[1]] %>%
  as_tibble() %>% 
  filter(expression_n >= 2) %>% 
  write_tsv('Methylation_neutral_win1000_step1000.txt')


# sliding window results --------------------------------------------------

ecotype_win = read_tsv('Methylation_eco_effect_win1000_step1000.txt')

ecotype_win %>% 
  group_by(Chromosome, 
           win_mid) %>% 
  summarize(sum_n = sum(expression_n)) %>%
  filter(sum_n >= 40) %>% View()

F2_win = read_tsv('Methylation_F2_effect_win1000_step1000.txt')

F2_win %>% 
  group_by(Chromosome, 
           win_mid) %>% 
  summarize(sum_n = sum(expression_n)) %>% 
  filter(sum_n >= 40)

F2_win %>% 
  filter(expression_mean < 0) %>% 
  group_by(Chromosome) %>% 
  summarize(n = n())
F2_win %>% 
  filter(expression_mean > 0)


F2_win = read_tsv('Methylation_F2_eco_effect_win1000_step1000.txt')

F2_win %>% 
  group_by(Chromosome, 
           win_mid) %>% 
  summarize(sum_n = sum(expression_n)) %>% 
  filter(sum_n >= 40)

F2_win %>% 
  filter(expression_mean < 0)
F2_win %>% 
  filter(expression_mean > 0)

F2_eco_win = read_tsv('Methylation_F2_eco_effect_win1000_step1000.txt')
F2_eco_win %>% 
  group_by(Chromosome, 
           win_mid) %>% 
  summarize(sum_n = sum(expression_n)) %>% 
  filter(sum_n >= 40)

F2_eco_win = read_tsv('Methylation_F2_eco_effect_win1000_step1000.txt')
F2_eco_win %>% 
  group_by(Chromosome, 
           win_mid) %>% 
  summarize(sum_n = sum(expression_n)) %>% 
  filter(sum_n >= 40)

F2_F2_win = read_tsv('Methylation_F2_F2_effect_win1000_step1000.txt')
F2_F2_win %>% 
  group_by(Chromosome, 
           win_mid) %>% 
  summarize(sum_n = sum(expression_n)) %>% 
  filter(sum_n >= 40)

F2_F2_eco_win = read_tsv('Methylation_F2_F2_eco_effect_win1000_step1000.txt')
F2_F2_eco_win %>% 
  group_by(Chromosome, 
           win_mid) %>% 
  summarize(sum_n = sum(expression_n)) %>% 
  filter(sum_n >= 40) %>% View()


# plot methylation sliding window -----------------------------------------

neutral = read_tsv('Methylation_neutral_win1000_step1000.txt')
F2_res = read_tsv('Methylation_F2_effect_win1000_step1000.txt') %>% 
  mutate(Effect = 'F1')
F2_res = read_tsv('Methylation_F2_effect_win1000_step1000.txt') %>% 
  mutate(effect = 'F2')
F2_F2_res = read_tsv('Methylation_F2_F2_effect_win1000_step1000.txt') %>% 
  mutate(effect = 'F1*F2')


ggplot(data = neutral, 
       aes(win_mid, 
           y = expression_mean, 
           group = Chromosome))+
  geom_point(col = '#8C8C8C', 
             size = 2)

F2_plot = ggplot()+
  geom_point(data = F2_res, 
             aes(x = win_mid, 
                 y = expression_mean, 
                 group = Chromosome), 
             col = '#184E77', 
             size = 2)+ 
  facet_grid(~Chromosome, 
             scales = 'free')


f2_plot = ggplot()+
  geom_point(data = F2_res, 
             aes(x = win_mid, 
                 y = expression_mean, 
                 group = Chromosome), 
             col = '#D9ED92', 
             size = 2)+
  facet_grid(~Chromosome, 
             scales = 'free')

F2_f2_plot = ggplot()+
  geom_point(data = F2_F2_res, 
             aes(x = win_mid, 
                 y = expression_mean, 
                 group = Chromosome), 
             col = '#52B69A', 
             size = 2)+
  facet_grid(~Chromosome, 
             scales = 'free')


F2_plot/f2_plot/F2_f2_plot



# bring in the ecotype ----------------------------------------------------



F1_eco_res = read_tsv('Methylation_F1_eco_effect_win1000_step1000.txt') %>% 
  mutate(Effect = 'F1', 
         abs_effect = abs(expression_mean)) %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  filter(expression_n > 3) %>%
  group_by(Chromosome) %>%
  mutate(threshold = quantile(abs_effect,
                              0.99,
                              na.rm = TRUE)) %>%
  filter(abs_effect >= threshold) %>% 
  mutate(direction = case_when(
    expression_mean > threshold ~ "Hyper",
    expression_mean < -threshold ~ "Hypo",
    TRUE ~ "None")) %>%
   # Order by genomic position
  arrange(Chromosome, win_mid)

## Peak detection code
F1_peaks = F1_eco_res %>%
  group_by(Chromosome, direction) %>%
  mutate(consecutive = cumsum(
      row_number() == 1 |
        win_mid != lag(win_mid) + 1000)) %>%
  ungroup() %>% 
  filter(direction != "None") %>%
  group_by(Chromosome, direction, consecutive) %>% 
  filter(consecutive > 3)

F1_hyper <- F1_peaks %>%
  filter(expression_mean >= threshold)
F1_hypo <- F1_peaks %>%
  filter(expression_mean <= -threshold)

F1_hypo_regions <- F1_peaks %>%
  filter(expression_mean <= -threshold) %>%
  arrange(Chromosome, win_start) %>%
  group_by(Chromosome) %>%
  mutate(
    gap = win_start - lag(win_start),
    peak_id = cumsum(
      is.na(gap) | gap > 10000)) %>%
  ungroup() %>% 
  mutate(Methylation_state = 'Hypo-methylated')

F1_hyper_regions <- F1_peaks %>%
  filter(expression_mean >= threshold) %>%
  arrange(Chromosome, win_start) %>%
  group_by(Chromosome) %>%
  mutate(
    gap = win_start - lag(win_start),
    peak_id = cumsum(
      is.na(gap) | gap > 10000)) %>%
  ungroup() %>% 
  mutate(Methylation_state = 'Hyper-methylated')

F1_Peak_regions = bind_rows(F1_hypo_regions, 
                            F1_hyper_regions) %>% 
  group_by(Chromosome) %>% 
  arrange(Chromosome, win_mid)

F1_hypo_summits <- F1_hypo_regions %>%
  group_by(Chromosome, peak_id) %>%
  slice_min(expression_mean, n = 1, with_ties = FALSE) %>%
  ungroup()

F1_hyper_summits <- F1_hyper_regions %>%
  group_by(Chromosome, peak_id) %>%
  slice_min(expression_mean, 
            n = 1, 
            with_ties = FALSE) %>%
  ungroup()

F1_Peak_regions %>% 
  dplyr::select(Chromosome, 
                win_start, 
                win_end) %>% 
  unite(col = genome_annotation, 
        c('win_start', 
          'win_end'), 
        sep = '..') %>% 
  unite(col = genome_annotation, 
        c('Chromosome', 
          'genome_annotation'), 
        sep = ':')


F2_eco_res = read_tsv('Methylation_F2_eco_effect_win1000_step1000.txt') %>% 
  mutate(Effect = 'F2', 
         abs_effect = abs(expression_mean))   %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  filter(expression_n > 3) %>%
  group_by(Chromosome) %>%
  mutate(threshold = quantile(abs_effect,
                              0.99,
                              na.rm = TRUE)) %>%
  filter(abs_effect >= threshold) %>% 
  mutate(direction = case_when(
    expression_mean > threshold ~ "Hyper",
    expression_mean < -threshold ~ "Hypo",
    TRUE ~ "None")) %>%
  # Order by genomic position
  arrange(Chromosome, win_mid)

## Peak detection code
F2_peaks = F2_eco_res %>%
  group_by(Chromosome, direction) %>%
  mutate(consecutive = cumsum(
    row_number() == 1 |
      win_mid != lag(win_mid) + 1000)) %>%
  ungroup() %>% 
  filter(direction != "None") %>%
  group_by(Chromosome, direction, consecutive) %>% 
  filter(consecutive > 3)

F2_hyper <- F2_peaks %>%
  filter(expression_mean >= threshold)
F2_hypo <- F2_peaks %>%
  filter(expression_mean <= -threshold)

F2_hypo_regions <- F2_peaks %>%
  filter(expression_mean <= -threshold) %>%
  arrange(Chromosome, win_start) %>%
  group_by(Chromosome) %>%
  mutate(
    gap = win_start - lag(win_start),
    peak_id = cumsum(
      is.na(gap) | gap > 10000)) %>%
  ungroup() %>% 
  mutate(Methylation_state = 'Hypo-methylated')

F2_hyper_regions <- F2_peaks %>%
  filter(expression_mean >= threshold) %>%
  arrange(Chromosome, win_start) %>%
  group_by(Chromosome) %>%
  mutate(
    gap = win_start - lag(win_start),
    peak_id = cumsum(
      is.na(gap) | gap > 10000)) %>%
  ungroup() %>% 
  mutate(Methylation_state = 'Hyper-methylated')

F2_Peak_regions = bind_rows(F2_hypo_regions, 
                            F2_hyper_regions) %>% 
  group_by(Chromosome) %>% 
  arrange(Chromosome, win_mid)

F2_hypo_summits <- F2_hypo_regions %>%
  group_by(Chromosome, peak_id) %>%
  slice_min(expression_mean, n = 1, with_ties = FALSE) %>%
  ungroup()

F2_hyper_summits <- F2_hyper_regions %>%
  group_by(Chromosome, peak_id) %>%
  slice_min(expression_mean, 
            n = 1, 
            with_ties = FALSE) %>%
  ungroup()

F2_Peak_regions %>% 
  dplyr::select(Chromosome, 
                win_start, 
                win_end) %>% 
  unite(col = genome_annotation, 
        c('win_start', 
          'win_end'), 
        sep = '..') %>% 
  unite(col = genome_annotation, 
        c('Chromosome', 
          'genome_annotation'), 
        sep = ':') %>% View()


F1_F2_eco_res = read_tsv('Methylation_F1_F2_eco_effect_win1000_step1000.txt') %>% 
  mutate(Effect = 'F1*F2', 
         abs_effect = abs(expression_mean)) %>% 
  filter(! Chromosome %in% c('chrUn')) %>% 
  filter(expression_n > 3) %>%
  group_by(Chromosome) %>%
  mutate(threshold = quantile(abs_effect,
                              0.99,
                              na.rm = TRUE)) %>%
  filter(abs_effect >= threshold) %>% 
  mutate(direction = case_when(
    expression_mean > threshold ~ "Hyper",
    expression_mean < -threshold ~ "Hypo",
    TRUE ~ "None")) %>%
  # Order by genomic position
  arrange(Chromosome, win_mid)

## Peak detection code
F1_F2_peaks = F1_F2_eco_res %>%
  group_by(Chromosome, direction) %>%
  mutate(consecutive = cumsum(
    row_number() == 1 |
      win_mid != lag(win_mid) + 1000)) %>%
  ungroup() %>% 
  filter(direction != "None") %>%
  group_by(Chromosome, direction, consecutive) %>% 
  filter(consecutive > 3)

F1_F2_hyper <- F1_F2_peaks %>%
  filter(expression_mean >= threshold)
F1_F2_hypo <- F1_F2_peaks %>%
  filter(expression_mean <= -threshold)

F1_F2_hypo_regions <- F1_F2_peaks %>%
  filter(expression_mean <= -threshold) %>%
  arrange(Chromosome, win_start) %>%
  group_by(Chromosome) %>%
  mutate(
    gap = win_start - lag(win_start),
    peak_id = cumsum(
      is.na(gap) | gap > 10000)) %>%
  ungroup() %>% 
  mutate(Methylation_state = 'Hypo-methylated')

F1_F2_hyper_regions <- F1_F2_peaks %>%
  filter(expression_mean >= threshold) %>%
  arrange(Chromosome, win_start) %>%
  group_by(Chromosome) %>%
  mutate(
    gap = win_start - lag(win_start),
    peak_id = cumsum(
      is.na(gap) | gap > 10000)) %>%
  ungroup() %>% 
  mutate(Methylation_state = 'Hyper-methylated')

F1_F2_Peak_regions = bind_rows(F1_F2_hypo_regions, 
                            F1_F2_hyper_regions) %>% 
  group_by(Chromosome) %>% 
  arrange(Chromosome, win_mid)

F1_F2_hypo_summits <- F1_F2_hypo_regions %>%
  group_by(Chromosome, peak_id) %>%
  slice_min(expression_mean, n = 1, with_ties = FALSE) %>%
  ungroup()

F1_F2_hyper_summits <- F1_F2_hyper_regions %>%
  group_by(Chromosome, peak_id) %>%
  slice_min(expression_mean, 
            n = 1, 
            with_ties = FALSE) %>%
  ungroup()


ggplot(data = neutral, 
       aes(win_mid, 
           y = expression_mean, 
           group = Chromosome))+
  geom_point(col = '#8C8C8C', 
             size = 2) +
  facet_grid(~Chromosome, 
             scales = 'free')
  
F1_methylation_state_cols = c('#003049', 
                           '#780000')

F1_eco_plot = ggplot()+
  geom_point(data = F1_eco_res, 
             aes(x = win_mid, 
                 y = expression_mean, 
                 group = Chromosome, 
                 colour = direction), 
             # col = '#184E77',
             size = 2)+ 
  scale_color_manual(values = F1_methylation_state_cols)+
  geom_hline(data = distinct(F1_eco_res, Chromosome, threshold),
             aes(yintercept = threshold),
             linetype = "dashed") +
  geom_hline(data = distinct(F1_eco_res, Chromosome, threshold),
             aes(yintercept = -threshold),
             linetype = "dashed") +
  facet_grid(~Chromosome, 
             scales = 'free')+
  theme(axis.text.x = element_text(angle = 90))

F1_eco_plot_peak = ggplot()+
  geom_point(data = F1_Peak_regions, 
             aes(x = win_mid, 
                 y = expression_mean, 
                 group = Chromosome, 
                 col = Methylation_state), 
             # col = '#184E77', 
             size = 2)+ 
  scale_color_manual(values = F1_methylation_state_cols)+
  geom_hline(data = distinct(F1_Peak_regions, Chromosome, threshold),
    aes(yintercept = threshold),
    linetype = "dashed") +
  geom_hline(data = distinct(F1_Peak_regions, Chromosome, threshold),
    aes(yintercept = -threshold),
    linetype = "dashed") +
  facet_grid(~Chromosome, 
             scales = 'free')+
  theme(axis.text.x = element_text(angle = 90))


F2_methylation_state_cols = c('#669bbc', 
                              '#c1121f')
f2_eco_plot = ggplot()+
  geom_point(data = F2_eco_res, 
             aes(x = win_mid, 
                 y = expression_mean, 
                 group = Chromosome, 
                 col = direction), 
             size = 2)+
  scale_colour_manual(values = F2_methylation_state_cols)+
  geom_hline(data = distinct(F2_eco_res, Chromosome, threshold),
             aes(yintercept = threshold),
             linetype = "dashed") +
  geom_hline(data = distinct(F2_eco_res, Chromosome, threshold),
             aes(yintercept = -threshold),
             linetype = "dashed") +
  facet_grid(~Chromosome, 
             scales = 'free')+
  theme(axis.text.x = element_text(angle = 90))

F2_eco_plot_peak = ggplot()+
  geom_point(data = F2_Peak_regions, 
             aes(x = win_mid, 
                 y = expression_mean, 
                 group = Chromosome, 
                 col = Methylation_state), 
             # col = '#184E77', 
             size = 2)+ 
  scale_color_manual(values = F2_methylation_state_cols)+
  geom_hline(data = distinct(F2_Peak_regions, Chromosome, threshold),
             aes(yintercept = threshold),
             linetype = "dashed") +
  geom_hline(data = distinct(F2_Peak_regions, Chromosome, threshold),
             aes(yintercept = -threshold),
             linetype = "dashed") +
  facet_grid(~Chromosome, 
             scales = 'free')+
  theme(axis.text.x = element_text(angle = 90))


F1_F2_methylation_state_cols = c('#0f4c5c', 
                              '#e36414')

f1_f2_eco_plot = ggplot()+
  geom_point(data = F1_F2_eco_res, 
             aes(x = win_mid, 
                 y = expression_mean, 
                 group = Chromosome, 
                 col = direction), 
             size = 2)+
  scale_colour_manual(values = F1_F2_methylation_state_cols)+
  geom_hline(data = distinct(F1_F2_eco_res, Chromosome, threshold),
             aes(yintercept = threshold),
             linetype = "dashed") +
  geom_hline(data = distinct(F1_F2_eco_res, Chromosome, threshold),
             aes(yintercept = -threshold),
             linetype = "dashed") +
  facet_grid(~Chromosome, 
             scales = 'free')+
  theme(axis.text.x = element_text(angle = 90))

F1_F2_eco_plot_peak = ggplot()+
  geom_point(data = F1_F2_Peak_regions, 
             aes(x = win_mid, 
                 y = expression_mean, 
                 group = Chromosome, 
                 col = Methylation_state), 
             # col = '#184E77', 
             size = 2)+ 
  scale_color_manual(values = F1_F2_methylation_state_cols)+
  geom_hline(data = distinct(F1_F2_Peak_regions, Chromosome, threshold),
             aes(yintercept = threshold),
             linetype = "dashed") +
  geom_hline(data = distinct(F1_F2_Peak_regions, Chromosome, threshold),
             aes(yintercept = -threshold),
             linetype = "dashed") +
  facet_grid(~Chromosome, 
             scales = 'free')+
  theme(axis.text.x = element_text(angle = 90))


F2_eco_plot/f2_eco_plot/F2_f2_eco_plot

ggsave('F2_F2_interaction_sliding_window_results_22.09.2026.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       height = 25, 
       width = 30)



# outliers div -------------------------------------------------------------

eco_res = read_tsv('Methylation_eco_effect_win1000_step1000.txt') %>% 
  filter(expression_n > 3) 

# %>%
  
  # 
  # filter(abs(expression_mean) > 0.1) %>% 
  # write_csv('Methylation_Ecotype_outlier_win1000_step1000_abs0.1.csv')



F2_eco_res = read_tsv('Methylation_F2_eco_effect_win1000_step1000.txt') %>% 
  mutate(Effect = 'F1')%>% 
  filter(expression_n > 3) 
# %>%
  # filter(abs(expression_mean) > 0.0004) %>% 
  # write_csv('Methylation_F1EcoRes_outlier_win1000_step1000_abs0.1.csv')
# %>%
  # group_by(Chromosome, 
  #          win_start, 
  #          win_end) %>% 
  # summarize(mean_expression = mean(expression_mean), 
  #           max_expression = max(expression_mean))



F2_eco_res = read_tsv('Methylation_F2_eco_effect_win1000_step1000.txt') %>% 
  mutate(effect = 'F2')%>% 
  filter(expression_n > 3) 
# %>%
#   filter(abs(expression_mean) > 0.1) %>% 
#   write_csv('Methylation_F2EcoRes_outlier_win1000_step1000_abs0.1.csv')



F2_F2_eco_res = read_tsv('Methylation_F2_F2_eco_effect_win1000_step1000.txt') %>% 
  mutate(effect = 'F1*F2')%>% 
  filter(expression_n > 3) 

# %>%
#   filter(abs(expression_mean) > 0.1)%>% 
#   write_csv('Methylation_F1F2EcoRes_outlier_win1000_step1000_abs0.1.csv')



# top 5% of expression distrubtion ----------------------------------------

eco_res_top_dist = eco_res[abs(eco_res$expression_mean) > quantile(eco_res$expression_mean, 
                                                            prob = 1-5/100),]


F2_eco_res_top_dist = F2_eco_res[abs(F2_eco_res$expression_mean) > quantile(F2_eco_res$expression_mean, 
                                                              prob = 1-5/100),]


f2_eco_res_top_dist = F2_eco_res[abs(F2_eco_res$expression_mean) > quantile(F2_eco_res$expression_mean, 
                                                                       prob = 1-5/100),]

F2_f2_eco_res_top_dist = F2_F2_eco_res[abs(F2_F2_eco_res$expression_mean) > quantile(F2_F2_eco_res$expression_mean, 
                                                                       prob = 1-5/100),]


# eco_res[abs(eco_res$expression_mean) > quantile(eco_res$expression_mean, 
#                                            prob = 1-5/100),]%>%
#   summarize(mean_expression = mean(expression_mean), 
#             max_expression = max(expression_mean), 
#             min_expression = min(expression_mean))

eco_res_top_dist %>% 
  write_csv('Methylation_Eco_outlier_top5dist_abs.csv')

F2_eco_res_top_dist %>% 
  write_csv('Methylation_F2_eco_outlier_top5dist_abs.csv')

f2_eco_res_top_dist %>% 
  write_csv('Methylation_F2_eco_outlier_top5dist_abs.csv')

F2_f2_eco_res_top_dist %>% 
  write_csv("Methylation_F2_F2_eco_outlier_top5dist_abs.csv")



# Classify outliers -------------------------------------------------------

F2_eco_res_out = read_csv('Methylation_F2_eco_outlier_top5dist_abs.csv') %>% 
  mutate(Effect = 'F1') %>% 
  filter(! Chromosome %in% c('chrUn',
                             'chrM')) %>% 
  filter(expression_n > 3) %>% 
  filter(expression_mean > 1.000000e-05 |
           expression_mean < -1.000000e-05)

F2_eco_res_out = read_csv('Methylation_F2_eco_outlier_top5dist_abs.csv') %>% 
  mutate(effect = 'F2')  %>% 
  filter(! Chromosome %in% c('chrUn',
                             'chrM')) %>% 
  filter(expression_n > 3)
F2_F2_eco_res_out = read_csv('Methylation_F2_F2_eco_outlier_top5dist_abs.csv') %>% 
  mutate(effect = 'F1*F2') %>% 
  filter(! Chromosome %in% c('chrUn',
                             'chrM')) %>% 
  filter(expression_n > 3)

F2_eco_res %>% 
  filter(Chromosome == 'chrXVII') %>% View()


ggplot(data = neutral, 
       aes(win_mid, 
           y = expression_mean, 
           group = Chromosome))+
  geom_point(col = '#8C8C8C', 
             size = 2) +
  facet_grid(~Chromosome, 
             scales = 'free')

F2_eco_out_plot = ggplot()+
  geom_point(data = F2_eco_res_out, 
             aes(x = win_mid, 
                 y = expression_mean, 
                 group = Chromosome), 
             col = '#184E77', 
             size = 2)+ 
  facet_grid(~Chromosome, 
             scales = 'free')+
  theme(axis.text.x = element_text(angle = 90))


f2_eco_out_plot = ggplot()+
  geom_point(data = F2_eco_res_out, 
             aes(x = win_mid, 
                 y = expression_mean, 
                 group = Chromosome), 
             col = '#D9ED92', 
             size = 2)+
  facet_grid(~Chromosome, 
             scales = 'free')+
  theme(axis.text.x = element_text(angle = 90))


F2_f2_eco_out_plot = ggplot()+
  geom_point(data = F2_F2_eco_res_out, 
             aes(x = win_mid, 
                 y = expression_mean, 
                 group = Chromosome), 
             col = '#52B69A', 
             size = 2)+
  facet_grid(~Chromosome, 
             scales = 'free')+
  theme(axis.text.x = element_text(angle = 90))



F2_eco_out_plot/f2_eco_out_plot/F2_f2_eco_out_plot

ggsave('F2_F2_interaction_sliding_window_results.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       height = 25, 
       width = 30)

