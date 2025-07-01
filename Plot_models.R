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


Body_WGP_out = read_csv('BODY_WGP_clean_RDA_outliers_AXIS1_RAW_PCaxes_methylation.csv')%>% 
  separate(col = loc, 
           into = c('CHR', 
                    'BP'), 
           sep = '-') %>% 
  arrange(CHR,
          BP) %>% 
  group_by(CHR)%>%
  mutate(BP = as.numeric(BP)) %>% 
  mutate(start = BP-100,   ## can change this to whatever window of interest you want around the site of interest
         end = BP+100)


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
  distinct() 


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
  dist_cal()


Body_WGP_out = read_csv('BODY_WGP_clean_RDA_outliers_AXIS1_RAW_PCaxes_methylation.csv')%>% 
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
  dist_cal()

Body_WGP_out_axisdf = axis_df(Body_WGP_out)


Fst_manhattan(outs = Body_WGP_out, 
              axisdf = Body_WGP_out_axisdf, 
              xval = BPcum, 
              yval = scores, 
              chr = Body_WGP_out$CHR, 
              out_col = '#439a86')


