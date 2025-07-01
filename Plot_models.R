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
  filter(feature == 'gene') %>% 
  select(position,
         chromosome,
         feature, 
         start, 
         end)

ensemlbe_annotation_data = gene_annotation %>% 
  filter(feature == 'gene') %>% 
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
  select(-trash) %>% 
  separate(ensemble_id, 
           into = c('ensemble_name', 
                    'trash'), 
           sep = '.CDS') %>% 
  select(-trash) 

annotation_data = bind_cols(gene_metadata, 
                            ensemble_annotation_genes) %>% 
  rename(CHR = chromosome, 
         BP = position)


# Body WGP genome data -------------------------------------------------------


Body_WGP_out = read_csv('BODY_WGP_clean_RDA_CAND_corr.csv')%>% 
  separate(col = loc, 
           into = c('CHR', 
                    'BP'), 
           sep = '-') %>% 
  arrange(CHR,
          BP) %>% 
  group_by(CHR)%>%
  mutate(BP = as.numeric(BP)) %>% 
  mutate(start = BP-100,   ## can change this to whatever window of interest you want around the site of interest
         end = BP+100) %>% 
  filter(Association == 'BODY_WGP_clean') %>% 
  select(CHR, 
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
           sep = '-') %>% 
  arrange(CHR,
          BP) %>% 
  group_by(CHR)%>%
  mutate(BP = as.numeric(BP)) %>% 
  mutate(start = BP-100,   ## can change this to whatever window of interest you want around the site of interest
         end = BP+100) %>% 
  filter(Association == 'BODY_TGP_clean') %>% 
  select(CHR, 
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
           sep = '-') %>% 
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
           sep = '-') %>% 
  filter(Association == 'BODY_WGP_clean') %>% 
  select(axis, 
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
           sep = '-') %>% 
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
           sep = '-') %>% 
  filter(Association == 'BODY_TGP_clean') %>% 
  select(axis, 
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

outs = BODY_TGP_COMBO %>% 
  filter(status == 'Outlier')
non_outs = BODY_TGP_COMBO %>% 
  filter(status == 'Neutral')

non_outs %>% 
  group_by(CHR) %>% 
  distinct()


BODY_TGP_manhattan = ggplot(non_outs, 
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


