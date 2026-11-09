
# install_github('tavareshugo/windowscanr')

library(tidyverse)
# library(dplyr)
# library(readr)
# library(future)
# library(future.apply)
library(windowscanr)

# Functions ---------------------------------------------------------------


dist_cal = function(data){
  data %>% 
    group_by(CHR) %>% 
    summarise(chr_len = max(POS)) %>% 
    mutate(total = cumsum(chr_len)-chr_len) %>% 
    dplyr::select(-chr_len) %>% 
    left_join(data, ., by = c('CHR'='CHR')) %>%
    arrange(CHR, 
            POS) %>% 
    mutate(BPcum = POS + total) 
}

stickle_CHR_reorder = function(data){
  data$CHR <- factor(data$CHR, 
                     levels=c("chr_I",
                              "chr_II", 
                              "chr_III", 
                              'chr_IV', 
                              'chr_V', 
                              'chr_VI', 
                              'chr_VII', 
                              'chr_VIII', 
                              'chr_IX', 
                              'chr_X', 
                              'chr_XI', 
                              'chr_XII', 
                              'chr_XIII', 
                              'chr_XIV', 
                              'chr_XV', 
                              'chr_XVI', 
                              'chr_XVII', 
                              'chr_XVIII', 
                              'chr_XIX', 
                              'chr_XX', 
                              'chr_XXI', 
                              'chr_Y', 
                              'chr_M', 
                              'chr_Un'))
  # data$CHR = as.character(data$CHR)
  return(data)
  
}


winScan_fixed <- function(
    x,
    groups = NULL,
    position = NULL,
    values,
    win_size,
    win_step = 0.5 * win_size,
    funs = "mean",
    cores = 1
) {
  
  funs <- as.list(funs)
  names(funs) <- sapply(funs, paste)
  
  if (is.null(groups)) {
    x$group <- 1
    groups <- "group"
  }
  
  if (is.null(position)) {
    x <- x %>%
      group_by(across(all_of(groups))) %>%
      mutate(pos = row_number()) %>%
      ungroup()
  } else {
    x <- x %>%
      rename(pos = all_of(position))
    
    stopifnot(is.numeric(x$pos))
  }
  
  if (win_size < win_step) {
    warning("Window size is smaller than step size!")
  }
  
  out <- x %>%
    group_by(across(all_of(groups))) %>%
    group_modify(~ windowscanr:::.winSlider(
      .x,
      values,
      win_size,
      win_step,
      funs,
      cores
    )) %>%
    as.data.frame()
  
  return(out)
}

# Warm cold Fst data ------------------------------------------------------


##Common FST outliers
WC_Fst_clean = read_csv('WC_Fst_clean.csv') %>%
  stickle_CHR_reorder() %>%
  dist_cal()


# Parallelizing winscan ---------------------------------------------------

# Set up 8 workers
plan(multisession, workers = 8)

# Split FST data by chromosome
chr_data <- split(WC_Fst_clean, WC_Fst_clean$CHR)

# Run sliding-window calculations in parallel
WC_fst_win <- future_lapply(
  chr_data,
  function(chr_df) {
    winScan_fixed(
      x = chr_df,
      groups = "CHR",
      position = "POS",
      values = "FST_zero",
      win_size = 1000,
      win_step = 1000,
      funs = c("mean", "sd")
    )
  },
  future.seed = TRUE
)

# Combine results
WC_fst_win <- bind_rows(WC_fst_win)

# Return to sequential processing
plan(sequential)

WC_fst_win %>%
  as_tibble() %>% 
  filter(FST_zero_n >= 3) %>% 
  write_tsv('WC_Fst_1kb_3obs_window.txt')

# Normal winscan ----------------------------------------------------------

WC_fst_win = winScan_fixed(
  x = WC_Fst_clean,
  groups = "CHR",
  position = "POS",
  values = "FST_zero",
  win_size = 1000,
  win_step = 1000,
  funs = c("mean", "sd")
)

WC_fst_win %>%
  as_tibble() %>% 
  filter(FST_n >= 3) %>% 
  write_tsv('WC_Fst_1kb_3obs_window.txt')



# windowed Fst data -------------------------------------------------------

WC_fst_win = read_tsv("WC_Fst_1kb_3obs_window.txt") %>% 
  dplyr::rename(Chromosome = CHR) %>% 
  separate(col = Chromosome, 
           into = c('chr', 
                    'num'), 
           sep = '_') %>% 
  unite(col = Chromosome, 
        c('chr', 
          'num'), 
        sep = '') %>% 
  mutate(fst_start = win_start,
    fst_end = win_end,
    buffer_start = pmax(1, win_start - 5000),
    buffer_end = win_end + 5000)

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


F2_eco_res = read_tsv('Methylation_F2_eco_effect_win1000_step1000.txt') %>% 
  mutate(Effect = 'F2', 
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


F1_F2_eco_res = read_tsv('Methylation_F1_F2_eco_effect_win1000_step1000.txt') %>% 
  mutate(Effect = 'F1_F2', 
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



## put a 5kb buffer around the fst values
## a methylation outlier window needs to be with 5kb away from an FST region
F1_Eco_res_FST = inner_join(WC_fst_win,
  F1_eco_res,
  by = join_by(Chromosome,
    buffer_start <= win_end,
    buffer_end >= win_start),
  suffix = c("_FST", "_Meth"),
  relationship = "many-to-many") %>% 
  filter(FST_zero_mean > 0.01)

## 4 methylation sites with high fst

F2_Eco_res_FST = inner_join(WC_fst_win,
                             F2_eco_res,
                             by = join_by(Chromosome,
                                          buffer_start <= win_end,
                                          buffer_end >= win_start),
                             suffix = c("_FST", "_Meth"),
                             relationship = "many-to-many")%>% 
  filter(FST_zero_mean > 0.01)

## 212 methylation sites with high Fst

F1_F2_Eco_res_FST = inner_join(WC_fst_win,
                             F1_F2_eco_res,
                             by = join_by(Chromosome,
                                          buffer_start <= win_end,
                                          buffer_end >= win_start),
                             suffix = c("_FST", "_Meth"),
                             relationship = "many-to-many")%>% 
  filter(FST_zero_mean > 0.01)

## 38 methylation sites with high Fst

# F1_Eco_res_FST = inner_join(WC_fst_win, 
#            F1_eco_res, 
#            by = c('Chromosome', 
#                   'win_start', 
#                   'win_end',
#                   'win_mid'))
# 
# F2_Eco_res_FST = inner_join(WC_fst_win, 
#                             F2_eco_res, 
#                             by = c('Chromosome', 
#                                    'win_start', 
#                                    'win_end',
#                                    'win_mid'))
# 
# F1_F2_Eco_res_FST = inner_join(WC_fst_win, 
#                             F1_F2_eco_res, 
#                             by = c('Chromosome', 
#                                    'win_start', 
#                                    'win_end',
#                                    'win_mid'))
