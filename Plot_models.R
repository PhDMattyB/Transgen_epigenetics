##############################
## Plotting RDA methylation model outputs
##
## Matt Brachmann (PhDMattyB)
##
## 01.07.2025
##
##############################

setwd('~/Methylation_data/')

library(tidyverse)
library(patchwork)
library(data.table)

# Genome annotation -------------------------------------------------------

## extract all of the stickle genome annotation data
gene_annotation = read_tsv('~/Parsons_Postdoc/Stickleback_Genomic/Stickleback_Annotation_features/stickleback_v5_ensembl_genes.gff3.gz', 
                           col_names = F, 
                           skip = 1) %>% 
  # filter(X3 %in% c('gene', 
  #                  'exon', 
  #                  'CDS')) %>% 
  group_by(X1) %>% 
  arrange(X4, 
          X5) %>% 
  ## arrange each gene by its start and end points on each chromosome
  mutate(mid = X4 + (X5-X4)/2) %>% 
  dplyr::select(X1, 
                X3:X5, 
                X9:mid) %>% 
  rename(chromosome = X1, 
         feature = X3, 
         start = X4, 
         end = X5, 
         gene_id = X9, 
         position = mid) %>% 
  na.omit()


gene_metadata = gene_annotation %>% 
  # filter(feature == 'gene') %>% 
  dplyr::select(position,
         chromosome,
         feature, 
         start, 
         end)

ensemlbe_annotation_data = gene_annotation %>% 
  # filter(feature == 'gene') %>% 
  pull(gene_id) %>% 
  as_tibble() %>% 
  separate(value, 
           into = c('ensemble_id', 
                    'gene_name',
                    'relationship'), 
           sep = ';')

ensemble_annotation_genes = ensemlbe_annotation_data %>% 
  separate(ensemble_id, 
           into = c('trash', 
                    'ensemble_id'), 
           sep = '=') %>%
  separate(gene_name, 
           into = c('trash', 
                    'gene_name'), 
           sep = '=') %>% 
  dplyr::select(-trash) %>% 
  separate(ensemble_id, 
           into = c('ensemble_name', 
                    'trash'), 
           sep = '.CDS') %>% 
  dplyr::select(-trash) 

annotation_data = bind_cols(gene_metadata, 
                            ensemble_annotation_genes) %>% 
  rename(CHR = chromosome, 
         BP = position)


# Body WGP genome data -------------------------------------------------------


Body_WGP_out = read_csv('BODY_WGP_clean_RDA_CAND_corr.csv')%>% 
  separate(col = loc, 
           into = c('CHR', 
                    'BP'), 
           sep = '-', 
           remove = F) %>% 
  arrange(CHR,
          BP) %>% 
  group_by(CHR)%>%
  mutate(BP = as.numeric(BP)) %>% 
  mutate(start = BP-100,   ## can change this to whatever window of interest you want around the site of interest
         end = BP+100) %>% 
  filter(Association == 'BODY_WGP_clean') %>% 
  select(loc,
         CHR, 
         BP, 
         scores, 
         start, 
         end)


setDT(Body_WGP_out)
setDT(annotation_data)

setkey(Body_WGP_out, 
       CHR, 
       start, 
       end)

## aligns the sites of interest with the annotated genome
Body_WGP_out_overlap = foverlaps(annotation_data,
                              Body_WGP_out,
                              # by.x = start,
                              # by.y = end,
                              type="any")

Body_WGP_out_overlap_tib = as_tibble(Body_WGP_out_overlap) %>% 
  na.omit() %>% 
  filter(CHR != 'chrUn') %>% 
  arrange(CHR, 
          BP)

Body_WGP_out_overlap_tib$gene_name %>%
  as_tibble() %>% 
  distinct() %>% 
  write_tsv('BODY_WGP_outlier_Genes_100bp_window.txt', 
            col_names = F)


Body_WGP_out_overlap_tib %>% 
  filter(CHR == 'chrXXI') %>% 
  distinct(gene_name)


# Body TGP genome data ----------------------------------------------------
Body_TGP_out = read_csv('BODY_TGP_clean_RDA_CAND_corr.csv')%>% 
  separate(col = loc, 
           into = c('CHR', 
                    'BP'), 
           sep = '-', 
           remove = F) %>% 
  arrange(CHR,
          BP) %>% 
  group_by(CHR)%>%
  mutate(BP = as.numeric(BP)) %>% 
  mutate(start = BP-100,   ## can change this to whatever window of interest you want around the site of interest
         end = BP+100) %>% 
  filter(Association == 'BODY_TGP_clean') %>% 
  select(loc,
         CHR, 
         BP, 
         scores, 
         start, 
         end)


setDT(Body_TGP_out)
setDT(annotation_data)

setkey(Body_TGP_out, 
       CHR, 
       start, 
       end)

## aligns the sites of interest with the annotated genome
Body_TGP_out_overlap = foverlaps(annotation_data,
                                 Body_TGP_out,
                                 # by.x = start,
                                 # by.y = end,
                                 type="any")

Body_TGP_out_overlap_tib = as_tibble(Body_TGP_out_overlap) %>% 
  na.omit() %>% 
  filter(CHR != 'chrUn') %>% 
  arrange(CHR, 
          BP)

Body_TGP_out_overlap_tib$gene_name %>%
  as_tibble() %>% 
  distinct() %>% 
  write_tsv('BODY_TGP_outlier_Genes_100bp_window.txt', 
            col_names = F)


Body_TGP_out_overlap_tib %>% 
  filter(CHR == 'chrXXI') %>% 
  distinct(gene_name)



# Body WGP out plot -------------------------------------------------------

Body_WGP_nonout = read_csv('BODY_WGP_clean_RDA_PCaxes_nonoutliers_methylation.csv')%>% 
  separate(col = loc, 
           into = c('CHR', 
                    'BP'), 
           sep = '-', 
           remove = F) %>% 
  arrange(CHR,
          BP) %>% 
  group_by(CHR)%>%
  mutate(BP = as.numeric(BP)) %>% 
  mutate(start = BP-100,   ## can change this to whatever window of interest you want around the site of interest
         end = BP+100) %>% 
  stickle_CHR_reorder2() %>% 
  rename(POS = BP) %>% 
  mutate(status = 'Neutral')


Body_WGP_out = read_csv('BODY_WGP_clean_RDA_CAND_corr.csv') %>% 
  separate(col = loc, 
           into = c('CHR', 
                    'BP'), 
           sep = '-', 
           remove = F) %>% 
  filter(Association == 'BODY_WGP_clean') %>% 
  select(loc, 
         axis, 
         CHR, 
         BP, 
         scores)%>% 
  arrange(CHR,
          BP) %>% 
  group_by(CHR)%>%
  mutate(BP = as.numeric(BP)) %>% 
  mutate(start = BP-100,   ## can change this to whatever window of interest you want around the site of interest
         end = BP+100) %>% 
  stickle_CHR_reorder2() %>% 
  rename(POS = BP) %>% 
  mutate(status = 'Outlier')


BODY_WGP_COMBO = bind_rows(Body_WGP_out, 
                           Body_WGP_nonout)%>% 
  dist_cal()

Body_WGP_axisdf = axis_df(BODY_WGP_COMBO)

outs = BODY_WGP_COMBO %>% 
  filter(status == 'Outlier')
non_outs = BODY_WGP_COMBO %>% 
  filter(status == 'Neutral')

non_outs %>% 
  group_by(CHR) %>% 
  distinct()


BODY_WGP_manhattan = ggplot(non_outs, 
       aes(x = POS, 
           y = scores))+
  # plot the non outliers in grey
  geom_point(aes(color = as.factor(CHR)), 
             alpha = 0.8, 
             size = 1.3)+
  ## alternate colors per chromosome
  scale_color_manual(values = rep(c("grey", "dimgrey"), 24))+
  ## plot the outliers on top of everything
  ## currently digging this hot pink colour
  geom_point(data = outs,
             col = '#fb8500',
             alpha=0.8, 
             size=1.3)+
  scale_x_continuous(label = Body_WGP_axisdf$CHR, 
                     breaks = Body_WGP_axisdf$center)+
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(-0.05,0.05))+
  facet_grid(~CHR, 
             scales = 'free')+
  # geom_hline(yintercept = 0.00043, 
  #            linetype = 2, 
  #            col = 'Black')+
  # ylim(0,1.0)+
  # scale_y_reverse(expand = c(0, 0))+
  # remove space between plot area and x axis
  labs(x = 'Cumulative base pair', 
       y = 'RDA score', 
       title = 'BODY WGP')+
  theme(legend.position="none",
        # panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        # axis.text.x = element_text(size = 9, 
        #                            angle = 90), 
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'))


ggsave('BODY_WGP_Outliers.svg', 
       plot = BODY_WGP_manhattan, 
       dpi = 'retina',
       units = 'cm',
       height = 15, 
       width = 30)  
ggsave('BODY_WGP_Outliers.tiff', 
       plot = BODY_WGP_manhattan, 
       dpi = 'retina',
       units = 'cm',
       height = 15, 
       width = 30)  



# Body TGP manhattan plot -------------------------------------------------

Body_TGP_nonout = read_csv('BODY_TGP_clean_RDA_PCaxes_nonoutliers_methylation.csv')%>% 
  separate(col = loc, 
           into = c('CHR', 
                    'BP'), 
           sep = '-', 
           remove = F) %>% 
  arrange(CHR,
          BP) %>% 
  group_by(CHR)%>%
  mutate(BP = as.numeric(BP)) %>% 
  mutate(start = BP-100,   ## can change this to whatever window of interest you want around the site of interest
         end = BP+100) %>% 
  stickle_CHR_reorder2() %>% 
  rename(POS = BP) %>% 
  mutate(status = 'Neutral')


Body_TGP_out = read_csv('BODY_TGP_clean_RDA_CAND_corr.csv') %>% 
  separate(col = loc, 
           into = c('CHR', 
                    'BP'), 
           sep = '-', 
           remove = F) %>% 
  filter(Association == 'BODY_TGP_clean') %>% 
  select(loc, 
         axis, 
         CHR, 
         BP, 
         scores)%>% 
  arrange(CHR,
          BP) %>% 
  group_by(CHR)%>%
  mutate(BP = as.numeric(BP)) %>% 
  mutate(start = BP-100,   ## can change this to whatever window of interest you want around the site of interest
         end = BP+100) %>% 
  stickle_CHR_reorder2() %>% 
  rename(POS = BP) %>% 
  mutate(status = 'Outlier')


BODY_TGP_COMBO = bind_rows(Body_TGP_out, 
                           Body_TGP_nonout)%>% 
  dist_cal()

Body_TGP_axisdf = axis_df(BODY_TGP_COMBO)

tgp_outs = BODY_TGP_COMBO %>% 
  filter(status == 'Outlier')
tgp_non_outs = BODY_TGP_COMBO %>% 
  filter(status == 'Neutral')

tgp_non_outs %>% 
  group_by(CHR) %>% 
  distinct()


BODY_TGP_manhattan = ggplot(tgp_non_outs, 
                            aes(x = POS, 
                                y = scores))+
  # plot the non outliers in grey
  geom_point(aes(color = as.factor(CHR)), 
             alpha = 0.8, 
             size = 1.3)+
  ## alternate colors per chromosome
  scale_color_manual(values = rep(c("grey", "dimgrey"), 24))+
  ## plot the outliers on top of everything
  ## currently digging this hot pink colour
  geom_point(data = tgp_outs,
             col = '#2a9d8f',
             alpha=0.8, 
             size=1.3)+
  scale_x_continuous(label = Body_TGP_axisdf$CHR, 
                     breaks = Body_TGP_axisdf$center)+
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(-0.05,0.05))+
  facet_grid(~CHR, 
             scales = 'free')+
  # geom_hline(yintercept = 0.00043, 
  #            linetype = 2, 
  #            col = 'Black')+
  # ylim(0,1.0)+
  # scale_y_reverse(expand = c(0, 0))+
  # remove space between plot area and x axis
  labs(x = 'Cumulative base pair', 
       y = 'RDA score', 
       title = 'BODY TGP')+
  theme(legend.position="none",
        # panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        # axis.text.x = element_text(size = 9, 
        #                            angle = 90), 
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'))


ggsave('BODY_TGP_Outliers.svg', 
       plot = BODY_TGP_manhattan, 
       dpi = 'retina',
       units = 'cm',
       height = 15, 
       width = 30)  
ggsave('BODY_TGP_Outliers.tiff', 
       plot = BODY_TGP_manhattan, 
       dpi = 'retina',
       units = 'cm',
       height = 15, 
       width = 30)  


# BODY TGP and WGP overlap ------------------------------------------------

BODY_TGP_WGP_overlap = inner_join(Body_WGP_out, 
           Body_TGP_out, 
           by = c('loc', 
                  'axis', 
                  'CHR', 
                  'POS', 
                  'start',
                  'end',
                  'status')) %>% 
  mutate(outlier_type = 'Overlap') %>% 
  select(-scores.y) %>% 
  rename(scores = scores.x)

## quantify unique TGP outliers
BODY_TGP_out_unique = Body_TGP_out %>% 
  anti_join(., 
            BODY_TGP_WGP_overlap, 
            by = 'loc') %>% 
  mutate(outlier_type = 'TGP unique')

## quantify unique WGP outliers
BODY_WGP_out_unique = Body_WGP_out %>% 
  anti_join(., 
            BODY_TGP_WGP_overlap, 
            by = 'loc') %>% 
  mutate(outlier_type = 'WGP unique')


All_outliers = bind_rows(BODY_TGP_WGP_overlap, 
                         BODY_TGP_out_unique, 
                         BODY_WGP_out_unique)


Neutral = read_csv('BODY_TGP_clean_RDA_PCaxes_nonoutliers_methylation.csv')%>% 
  separate(col = loc, 
           into = c('CHR', 
                    'BP'), 
           sep = '-', 
           remove = F) %>% 
  arrange(CHR,
          BP) %>% 
  group_by(CHR)%>%
  mutate(BP = as.numeric(BP)) %>% 
  mutate(start = BP-100,   ## can change this to whatever window of interest you want around the site of interest
         end = BP+100) %>% 
  anti_join(.,
            All_outliers, 
            by = 'loc') %>% 
  stickle_CHR_reorder2() %>% 
  rename(POS = BP) %>% 
  mutate(status = 'Neutral')


Big_data_set = bind_rows(All_outliers, 
                         Neutral)%>% 
  dist_cal()

axisdf = axis_df(Big_data_set)

outliers = Big_data_set %>% 
  filter(status == 'Outlier')
neutral_meth = Big_data_set %>% 
  filter(status == 'Neutral')

overlap_outs = outliers %>% 
  filter(outlier_type == 'Overlap')

tgp_outs = outliers %>% 
  filter(outlier_type == 'TGP unique')

wgp_outs = outliers %>% 
  filter(outlier_type == 'WGP unique')


# manhat_cols = c('#4895ef', 
#                 '#b5179e', 
#                 '#f72585')

 # big_man_ting = ggplot(neutral_meth, 
 #       aes(x = POS, 
 #           y = scores))+
 #  # plot the non outliers in grey
 #  geom_point(aes(color = as.factor(CHR)), 
 #             alpha = 0.8, 
 #             size = 1.3)+
 #  ## alternate colors per chromosome
 #  scale_color_manual(values = rep(c("grey", "dimgrey"), 24))+
 #  ## plot the outliers on top of everything
 #  ## currently digging this hot pink colour
 #  geom_point(data = overlap_outs,
 #             col = '#4895ef',
 #             alpha=0.8, 
 #             size=1.3)+
 #  geom_point(data = wgp_outs, 
 #             col = '#b5179e', 
 #             alpha = 0.8, 
 #             size = 1.3)+
 #  geom_point(data = tgp_outs, 
 #             col = '#f72585', 
 #             alpha = 0.8, 
 #             size = 1.3)+
 #  # scale_color_manual(values = manhat_cols)+
 #  scale_x_continuous(label = axisdf$CHR, 
 #                     breaks = axisdf$center)+
 #  scale_y_continuous(expand = c(0, 0), 
 #                     limits = c(-0.05,0.05))+
 #  facet_grid(~CHR, 
 #             scales = 'free')+
 #  # geom_hline(yintercept = 0.00043, 
 #  #            linetype = 2, 
 #  #            col = 'Black')+
 #  # ylim(0,1.0)+
 #  # scale_y_reverse(expand = c(0, 0))+
 #  # remove space between plot area and x axis
 #  labs(x = 'Cumulative base pair', 
 #       y = 'RDA score')+
 #  theme(legend.position="none",
 #        # panel.border = element_blank(),
 #        panel.grid.major.x = element_blank(),
 #        panel.grid.minor.x = element_blank(),
 #        axis.text.x = element_blank(), 
 #        axis.ticks.x = element_blank(),
 #        # axis.text.x = element_text(size = 9, 
 #        #                            angle = 90), 
 #        axis.title = element_text(size = 14),
 #        axis.title.x = element_blank(),
 #        axis.text.y = element_text(size = 12), 
 #        strip.background = element_rect(fill = 'white'), 
 #        strip.text = element_text(face = 'bold'))
 # 
 # 
 # ggsave('BODY_Big_manhattan_plot.svg', 
 #        plot = big_man_ting, 
 #        dpi = 'retina',
 #        units = 'cm',
 #        height = 15, 
 #        width = 30)  
 # ggsave('Body_Big_manhattan_plot.tiff', 
 #        plot = big_man_ting, 
 #        dpi = 'retina',
 #        units = 'cm',
 #        height = 15, 
 #        width = 30)  
 # 




big_man_overlap = ggplot(neutral_meth,
      aes(x = POS,
          y = scores))+
 # plot the non outliers in grey
 geom_point(aes(color = as.factor(CHR)),
            alpha = 0.8,
            size = 1.3)+
 ## alternate colors per chromosome
 scale_color_manual(values = rep(c("grey", "dimgrey"), 24))+
 ## plot the outliers on top of everything
 ## currently digging this hot pink colour
 geom_point(data = overlap_outs,
            col = '#4895ef',
            alpha=0.8,
            size=1.3)+
 # geom_point(data = wgp_outs,
 #            col = '#b5179e',
 #            alpha = 0.8,
 #            size = 1.3)+
 # geom_point(data = tgp_outs,
 #            col = '#f72585',
 #            alpha = 0.8,
 #            size = 1.3)+
 # scale_color_manual(values = manhat_cols)+
 scale_x_continuous(label = axisdf$CHR,
                    breaks = axisdf$center)+
 scale_y_continuous(expand = c(0, 0),
                    limits = c(-0.05,0.05))+
 facet_grid(~CHR,
            scales = 'free')+
 # geom_hline(yintercept = 0.00043,
 #            linetype = 2,
 #            col = 'Black')+
 # ylim(0,1.0)+
 # scale_y_reverse(expand = c(0, 0))+
 # remove space between plot area and x axis
 labs(x = 'Cumulative base pair',
      y = 'RDA score')+
 theme(legend.position="none",
       # panel.border = element_blank(),
       panel.grid.major.x = element_blank(),
       panel.grid.minor.x = element_blank(),
       axis.text.x = element_blank(),
       axis.ticks.x = element_blank(),
       # axis.text.x = element_text(size = 9,
       #                            angle = 90),
       axis.title = element_text(size = 14),
       axis.title.x = element_blank(),
       axis.text.y = element_text(size = 12),
       strip.background = element_rect(fill = 'white'),
       strip.text = element_text(face = 'bold'))

big_man_wgp = ggplot(neutral_meth,
      aes(x = POS,
          y = scores))+
 # plot the non outliers in grey
 geom_point(aes(color = as.factor(CHR)),
            alpha = 0.8,
            size = 1.3)+
 ## alternate colors per chromosome
 scale_color_manual(values = rep(c("grey", "dimgrey"), 24))+
 ## plot the outliers on top of everything
 ## currently digging this hot pink colour
 # geom_point(data = overlap_outs,
 #            col = '#4895ef',
 #            alpha=0.8,
 #            size=1.3)+
 geom_point(data = wgp_outs,
            col = '#b5179e',
            alpha = 0.8,
            size = 1.3)+
 # geom_point(data = tgp_outs,
 #            col = '#f72585',
 #            alpha = 0.8,
 #            size = 1.3)+
 # scale_color_manual(values = manhat_cols)+
 scale_x_continuous(label = axisdf$CHR,
                    breaks = axisdf$center)+
 scale_y_continuous(expand = c(0, 0),
                    limits = c(-0.05,0.05))+
 facet_grid(~CHR,
            scales = 'free')+
 # geom_hline(yintercept = 0.00043,
 #            linetype = 2,
 #            col = 'Black')+
 # ylim(0,1.0)+
 # scale_y_reverse(expand = c(0, 0))+
 # remove space between plot area and x axis
 labs(x = 'Cumulative base pair',
      y = 'RDA score')+
 theme(legend.position="none",
       # panel.border = element_blank(),
       panel.grid.major.x = element_blank(),
       panel.grid.minor.x = element_blank(),
       axis.text.x = element_blank(),
       axis.ticks.x = element_blank(),
       # axis.text.x = element_text(size = 9,
       #                            angle = 90),
       axis.title = element_text(size = 14),
       axis.title.x = element_blank(),
       axis.text.y = element_text(size = 12),
       strip.background = element_rect(fill = 'white'),
       strip.text = element_text(face = 'bold'))

big_man_tgp = ggplot(neutral_meth,
      aes(x = POS,
          y = scores))+
 # plot the non outliers in grey
 geom_point(aes(color = as.factor(CHR)),
            alpha = 0.8,
            size = 1.3)+
 ## alternate colors per chromosome
 scale_color_manual(values = rep(c("grey", "dimgrey"), 24))+
 ## plot the outliers on top of everything
 ## currently digging this hot pink colour
 # geom_point(data = overlap_outs,
 #            col = '#4895ef',
 #            alpha=0.8,
 #            size=1.3)+
 # geom_point(data = wgp_outs,
 #            col = '#b5179e',
 #            alpha = 0.8,
 #            size = 1.3)+
 geom_point(data = tgp_outs,
            col = '#f72585',
            alpha = 0.8,
            size = 1.3)+
 # scale_color_manual(values = manhat_cols)+
 scale_x_continuous(label = axisdf$CHR,
                    breaks = axisdf$center)+
 scale_y_continuous(expand = c(0, 0),
                    limits = c(-0.05,0.05))+
 facet_grid(~CHR,
            scales = 'free')+
 # geom_hline(yintercept = 0.00043,
 #            linetype = 2,
 #            col = 'Black')+
 # ylim(0,1.0)+
 # scale_y_reverse(expand = c(0, 0))+
 # remove space between plot area and x axis
 labs(x = 'Cumulative base pair',
      y = 'RDA score')+
 theme(legend.position="none",
       # panel.border = element_blank(),
       panel.grid.major.x = element_blank(),
       panel.grid.minor.x = element_blank(),
       axis.text.x = element_blank(),
       axis.ticks.x = element_blank(),
       # axis.text.x = element_text(size = 9,
       #                            angle = 90),
       axis.title = element_text(size = 14),
       axis.title.x = element_blank(),
       axis.text.y = element_text(size = 12),
       strip.background = element_rect(fill = 'white'),
       strip.text = element_text(face = 'bold'))

big_man_ting = big_man_overlap/big_man_wgp/big_man_tgp


# ggsave('BODY_Big_manhattan_plot.svg', 
#        plot = big_man_ting, 
#        dpi = 'retina',
#        units = 'cm',
#        height = 15, 
#        width = 30)  
ggsave('Body_Big_manhattan_plot.tiff',
       plot = big_man_ting,
       dpi = 'retina',
       units = 'cm',
       height = 30,
       width = 30)


# Filter based on position ------------------------------------------------
overlap_outs = outliers %>%
  filter(outlier_type == 'Overlap')
wgp_outs = outliers %>%
  filter(outlier_type == 'WGP unique')
tgp_outs = outliers %>%
  filter(outlier_type == 'TGP unique')

# unique vs overlap GO terms ----------------------------------------------

## overlap outlier gene names
overlap_outs = outliers %>% 
  filter(outlier_type == 'Overlap')

setDT(overlap_outs)
setDT(annotation_data)

setkey(overlap_outs, 
       CHR, 
       start, 
       end)

## aligns the sites of interest with the annotated genome
overlap_genes = foverlaps(annotation_data,
                                 overlap_outs,
                                 # by.x = start,
                                 # by.y = end,
                                 type="any")

overlap_genes_tib = as_tibble(overlap_genes) %>% 
  na.omit() %>% 
  filter(CHR != 'chrUn') %>% 
  arrange(CHR, 
          BP)

overlap_genes_tib$gene_name %>%
  as_tibble() %>% 
  distinct() %>% 
  write_tsv('BODY_TGP_WGP_overlap_outlier_Genes_100bp_window.txt', 
            col_names = F)


## WGP outliers unique to WGP
wgp_outs = outliers %>% 
  filter(outlier_type == 'WGP unique')

setDT(wgp_outs)
setDT(annotation_data)

setkey(wgp_outs, 
       CHR, 
       start, 
       end)

## aligns the sites of interest with the annotated genome
WGP_genes = foverlaps(annotation_data,
                          wgp_outs,
                          # by.x = start,
                          # by.y = end,
                          type="any")

WGP_genes_tib = as_tibble(WGP_genes) %>% 
  na.omit() %>% 
  filter(CHR != 'chrUn') %>% 
  arrange(CHR, 
          BP)

WGP_genes_tib$gene_name %>%
  as_tibble() %>% 
  distinct() %>% 
  write_tsv('BODY_WGP_unique_outlier_Genes_100bp_window.txt', 
            col_names = F)

## tgp unique genes
tgp_outs = outliers %>% 
  filter(outlier_type == 'TGP unique')

setDT(tgp_outs)
setDT(annotation_data)

setkey(tgp_outs, 
       CHR, 
       start, 
       end)

## aligns the sites of interest with the annotated genome
TGP_genes = foverlaps(annotation_data,
                      tgp_outs,
                      # by.x = start,
                      # by.y = end,
                      type="any")

TGP_genes_tib = as_tibble(TGP_genes) %>% 
  na.omit() %>% 
  filter(CHR != 'chrUn') %>% 
  arrange(CHR, 
          BP)

TGP_genes_tib$gene_name %>%
  as_tibble() %>% 
  distinct() %>% 
  write_tsv('BODY_TGP_unique_outlier_Genes_100bp_window.txt', 
            col_names = F)




