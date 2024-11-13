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

F1_array = array(0, dim = c(37, 2, 1575))

for(i in F1_12deg_range){
  F1_array[,,i] = gpa$coords[,,i] - F1_12deg_array[,,1]
}

for(i in F1_18deg_range){
  F1_array[,,i] = gpa$coords[,,i] - F1_18deg_array[,,1]
}


F1_array_consensus = array(0, dim = c(37, 2, 1575))

for(i in 1:1575){
  F1_array_consensus[,,i] = F1_array[,,i] + mean_shape_array[,,1]
}

# writeland.tps(F1_array_consensus, 
#               file = 'F1_effect_landmarks_all_individuals.tps',
#               scale = NULL, 
#               specID = T)



# Isolate F2 effects ------------------------------------------------------


F2_temp_mod = procD.lm(gpa$coords ~ meta_data$F2, 
                       iter = 999, 
                       RRPP = T)


F2_fitted_12deg = F2_temp_mod$GM$fitted[,,1]
F2_fitted_mat_12deg = as.matrix(F2_fitted_12deg)
F2_12deg_array = array(F2_fitted_mat_12deg, 
                       dim = c(37, 2, 1))

F2_fitted_18deg = F2_temp_mod$GM$fitted[,,35]
F2_fitteed_18deg_mat = as.matrix(F2_fitted_18deg)
F2_18deg_array = array(F2_fitteed_18deg_mat, 
                       dim = c(37, 2, 1))

F2_12deg_range = c(1:50, 105:155, 206:259, 310:359, 412:461, 512:563, 612:650, 701:746, 798:847, 898:947, 998:1047,
                   1098:1147, 1198:1249, 1298:1347, 1398:1447, 1475:1525)

F2_18deg_range = c(51:104, 156:205, 260:309, 360:411, 462:511, 564:611, 651:700, 747:797, 848:897, 948:997,
                   1048:1097, 1148:1197, 1250:1297, 1348:1397, 1448:1474, 1525:1575)

F2_array = array(0, dim = c(37, 2, 1575))

for(i in F2_12deg_range){
  F2_array[,,i] = gpa$coords[,,i] - F2_12deg_array[,,1]
}

for(i in F2_18deg_range){
  F2_array[,,i] = gpa$coords[,,i] - F2_18deg_array[,,1]
}


F2_array_consensus = array(0, dim = c(37, 2, 1575))

for(i in 1:1575){
  F2_array_consensus[,,i] = F2_array[,,i] + mean_shape_array[,,1]
}

# writeland.tps(F2_array_consensus,
#               file = 'F2_effect_landmarks_all_individuals.tps',
#               scale = NULL,
#               specID = T)


# isolate ecotype effects -------------------------------------------------

## isolate effects due to warm cold divergence axis regardless of population
ecotype_mod1 = procD.lm(gpa$coords ~ meta_data$ecotype...12, 
                       iter = 999, 
                       RRPP = T)

ecotype_lm = ecotype_mod1$GM$fitted

ecotype_consensus = array(0, dim = c(37, 2, 1575))

for(i in 1:1575){
  ecotype_consensus[,,i] = ecotype_lm[,,i] + mean_shape_array[,,1]
}

writeland.tps(ecotype_consensus,
              file = 'ecotype_effect_landmarks_all_individuals.tps',
              scale = NULL,
              specID = T)

ecotype_mod2 = procD.lm(gpa$coords ~ meta_data$poppair*meta_data$ecotype...12, 
                        iter = 999, 
                        RRPP = T)

ecotype2_fitted = ecotype_mod2$GM$fitted

ecotype2_consensus = array(0, dim = c(37, 2, 1575))

for(i in 1:1575){
  ecotype2_consensus[,,i] = ecotype2_fitted[,,i] + mean_shape_array[,,1]
}

writeland.tps(ecotype2_consensus,
              file = 'ecotype_effect_per_population_landmarks_all_individuals.tps',
              scale = NULL,
              specID = T)




# Phenotype PCA -----------------------------------------------------------



## pca of the raw landmark data
raw = readland.tps('F2_All_aligned_withsliders.tps', 
                   specID = 'imageID', 
                   readcurves = T)

sliders = define.sliders(c(28:37,1))

raw_gpa = gpagen(raw, 
             print.progress = T, 
             curves = sliders)

raw_pca = gm.prcomp(A = raw_gpa$coords)

summary(raw_pca)

raw_pca_vals = raw_pca$x %>% 
  as_tibble() %>% 
  select(1:5)

raw_pca_data = bind_cols(meta_data, 
      raw_pca_vals) 

raw_pca_data$F1 = as.character(raw_pca_data$F1)
raw_pca_data$F2 = as.character(raw_pca_data$F2)

ggplot(data = raw_pca_data)+
  geom_point(aes(x = Comp1, 
                 y = Comp2, 
                 col = F1, 
                 shape = F2))


## pca of the f1 effects
F1_effects = readland.tps('F1_effect_landmarks_all_individuals.tps', 
                          specID = 'imageID', 
                          readcurves = T)

f1_gpa = gpagen(F1_effects, 
                curves = sliders)
f1_pca = gm.prcomp(f1_gpa$coords)

summary(f1_pca)

f1_pca_vals = f1_pca$x %>% 
  as_tibble() %>% 
  select(1:5)

f1_pca_data = bind_cols(meta_data, 
                        f1_pca_vals)

f1_pca_data$F1 = as.character(f1_pca_data$F1)
f1_pca_data$F2 = as.character(f1_pca_data$F2)

ggplot(data = f1_pca_data)+
  geom_point(aes(x = Comp1, 
                 y = Comp2, 
                 col = F1, 
                 shape = F2))


## pca of the f2 effects
f2_effects = readland.tps('F2_effect_landmarks_all_individuals.tps', 
                          specID = 'imageID', 
                          readcurves = T)

f2_gpa = gpagen(f2_effects, 
                curves = sliders)
f2_pca = gm.prcomp(f2_gpa$coords)

summary(f2_pca)

f2_pca_vals = f2_pca$x %>% 
  as_tibble() %>% 
  select(1:5)

f2_pca_data = bind_cols(meta_data, 
                        f2_pca_vals)

f2_pca_data$F1 = as.character(f2_pca_data$F1)
f2_pca_data$F2 = as.character(f2_pca_data$F2)

ggplot(data = f2_pca_data)+
  geom_point(aes(x = Comp1, 
                 y = Comp2, 
                 col = F1, 
                 shape = F2))


## pca of the cold vs warm ecotype effects
eco1_effects = readland.tps('ecotype_effect_per_population_landmarks_all_individuals.tps', 
                          specID = 'imageID', 
                          readcurves = T)

eco1_gpa = gpagen(eco1_effects, 
                curves = sliders)
eco1_pca = gm.prcomp(eco1_gpa$coords)

summary(eco1_pca)

eco1_pca_vals = eco1_pca$x %>% 
  as_tibble() %>% 
  select(1:5)

eco1_pca_data = bind_cols(meta_data, 
                        eco1_pca_vals)

eco1_pca_data$F1 = as.character(eco1_pca_data$F1)
eco1_pca_data$F2 = as.character(eco1_pca_data$F2)

ggplot(data = eco1_pca_data)+
  geom_point(aes(x = Comp1, 
                 y = Comp2, 
                 col = F1, 
                 shape = F2))

ggplot(data = eco1_pca_data)+
  geom_jitter(aes(x = Comp1, 
                  y = Comp2, 
                  col = F1, 
                  shape = F2))
