##############################
## Meth paper - transgenerational phenotypes
##
## Matt Brachmann (PhDMattyB)
##
## 12.11.2024
##
##############################


setwd('~/Parsons_Postdoc/Methylation_data/')
setwd('~/Methylation_data/')

library(tidyverse)
library(geomorph)


# metadata ----------------------------------------------------------------

meta_data = read_csv('formattedDataEU.csv') %>% 
  select(fish, 
             ID, 
             Full_ID, 
             F1, 
             F2, 
             poppair, 
             ecotype...12, 
         csize_real)


# Body shape data ---------------------------------------------------------

bs_data = readland.tps('F2_All_aligned_withsliders.tps', 
                       specID = 'imageID', 
                       readcurves = T)

# bs_data = read_csv("epigenetic_landmark_data.csv") %>% 
#   select(-na)
# 
# bs_data = as.matrix(bs_data)
# 
# arrayspecs(A = bs_data, 
#            p = 37, 
#            k = 2)


sliders = define.sliders(c(28:37,1))

gpa = gpagen(bs_data, 
             print.progress = T, 
             curves = sliders)



# test for allometry ------------------------------------------------------

allo_mod = procD.lm(gpa$coords ~ log(meta_data$csize_real), 
                    iter = 999, 
                    RRPP = T)

summary(allo_mod)

## Allometry is not significant. 


# mean shape data ---------------------------------------------------------

mean_shape = mshape(gpa$coords)
matrix_mean_shape = as.matrix(mean_shape)
mean_shape_array = array(matrix_mean_shape, 
                         dim = c(37, 2, 1))



# Isolate F1 effects ------------------------------------------------------

F1_temp_mod = procD.lm(gpa$coords ~ meta_data$F1, 
                       iter = 999, 
                       RRPP = T)

F1_fitted_12deg = F1_temp_mod$GM$fitted[,,1]
F1_fitted_mat_12deg = as.matrix(F1_fitted_12deg)
F1_12deg_array = array(F1_fitted_mat_12deg, 
                       dim = c(37, 2, 1))

F1_fitted_18deg = F1_temp_mod$GM$fitted[,,64]
F1_fitteed_18deg_mat = as.matrix(F1_fitted_18deg)
F1_18deg_array = array(F1_fitteed_18deg_mat, 
                       dim = c(37, 2, 1))

F1_12deg_range = c(1:104, 206:370, 411:511, 612:700, 789:897, 999:1097, 1198:1297, 1398:1474)

F1_18deg_range = c(105:205, 371:410, 512:611, 701:788, 898:998, 1098:1197, 1298:1297, 1475:1575)




