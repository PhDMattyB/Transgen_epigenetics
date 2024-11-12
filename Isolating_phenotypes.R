##############################
## Meth paper - transgenerational phenotypes
##
## Matt Brachmann (PhDMattyB)
##
## 12.11.2024
##
##############################


setwd('~/Parsons_Postdoc/Methylation_data/')

library(tidyverse)
library(geomorph)


# Body shape data ---------------------------------------------------------

bs_data = readland.tps('F2_All_aligned_withsliders.tps', 
                       specID = 'imageID', 
                       readcurves = T)

sliders = define.sliders(bs_data, 
               10)

gpa = gpagen(bs_data, 
             print.progress = T, 
             curves = sliders)
