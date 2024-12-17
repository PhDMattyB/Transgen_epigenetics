##############################
## Meth paper - transgenerational phenotypes
##
## Matt Brachmann (PhDMattyB)
##
## 12.11.2024
##
##############################

## Office computer
setwd('~/Parsons_Postdoc/Methylation_data/')

## personal laptop
setwd('~/Methylation_data/')


library(tidyverse)
library(geomorph)
library(PCDimension)


# metadata ----------------------------------------------------------------

meta_data = read_csv('formattedDataEU.csv') %>% 
  select(fish, 
             ID, 
             Full_ID, 
             F1, 
             F2, 
             poppair, 
             ecotype...13, 
         csize_real) %>% 
  rename(ecotype = ecotype...13)
# %>% 
#   unite(Ecotype, 
#         c(poppair, 
#           ecotype...12))


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

allo_mod2 = procD.lm(gpa$coords ~ log(meta_data$csize_real)*meta_data$ecotype*meta_data$poppair, 
                    iter = 999, 
                    RRPP = T)

summary(allo_mod2)
## Theres an ecotype*lake allometric scaling effect
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

writeland.tps(F1_array_consensus,
              file = 'F1_effect_landmarks_all_individuals.tps',
              scale = NULL,
              specID = T)


# F1 by ecotype effects ---------------------------------------------------
F1_temp_mod_eco = procD.lm(gpa$coords ~ meta_data$F1*meta_data$ecotype, 
                       iter = 999, 
                       RRPP = T)

F1_fitted_12deg_cold = F1_temp_mod_eco$GM$fitted[,,1]
F1_fitted_12deg_cold_mat = as.matrix(F1_fitted_12deg_cold)
F1_12deg_cold_array = array(F1_fitted_12deg_cold_mat, 
                            dim = c(37, 2, 1))

F1_fitted_12deg_warm = F1_temp_mod_eco$GM$fitted[,,206]
F1_fitted_12deg_warm_mat = as.matrix(F1_fitted_12deg_warm)
F1_12deg_warm_array = array(F1_fitted_12deg_warm_mat, 
                            dim = c(37, 2, 1))

F1_fitted_18deg_cold = F1_temp_mod_eco$GM$fitted[,,105]
F1_fitted_18deg_cold_mat = as.matrix(F1_fitted_18deg_cold)
F1_18deg_cold_array = array(F1_fitted_18deg_cold_mat, 
                            dim = c(37, 2, 1))

F1_fitted_18deg_warm = F1_temp_mod_eco$GM$fitted[,,310]
F1_fitted_18deg_warm_mat = as.matrix(F1_fitted_18deg_warm)
F1_18deg_warm_array = array(F1_fitted_18deg_warm_mat, 
                            dim = c(37, 2, 1))

F1_12deg_cold_range = c(1:104, 412:511, 798:897, 1198:1297)
F1_12deg_warm_range = c(206:309, 612:700, 998:1097, 1398:1474)
F1_18deg_cold_range = c(105:205, 512:611, 898:997, 1298:1397)
F1_18deg_warm_range = c(310:411, 701:797,1098:1197, 1475:1575)

F1_temp_eco_array = array(0,
                          dim = c(37, 2, 1575))

for(i in F1_12deg_cold_range){
  F1_temp_eco_array[,,i] = gpa$coords[,,i] - F1_12deg_cold_array[,,1]
}

for(i in F1_12deg_warm_range){
  F1_temp_eco_array[,,i] = gpa$coords[,,i] - F1_12deg_warm_array[,,1]
}

for(i in F1_18deg_cold_range){
  F1_temp_eco_array[,,i] = gpa$coords[,,i] - F1_18deg_cold_array[,,1]
}

for(i in F1_18deg_warm_range){
  F1_temp_eco_array[,,i] = gpa$coords[,,i] - F1_18deg_warm_array[,,1]
}

F1_array_consensus = array(0, dim = c(37, 2, 1575))

for(i in 1:1575){
  F1_array_consensus[,,i] = F1_temp_eco_array[,,i] + mean_shape_array[,,1]
}

writeland.tps(F1_array_consensus,
              file = 'TGP_ecotype_variation_landmarks.tps',
              scale = NULL,
              specID = T)




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

writeland.tps(F2_array_consensus,
              file = 'F2_effect_landmarks_all_individuals.tps',
              scale = NULL,
              specID = T)


# F2 effects by ecotype ---------------------------------------------------

F2_temp_mod_eco = procD.lm(gpa$coords ~ meta_data$F2*meta_data$ecotype, 
                           iter = 999, 
                           RRPP = T)

F2_fitted_12deg_cold = F2_temp_mod_eco$GM$fitted[,,1]
F2_fitted_12deg_cold_mat = as.matrix(F2_fitted_12deg_cold)
F2_12deg_cold_array = array(F2_fitted_12deg_cold_mat, 
                            dim = c(37, 2, 1))

F2_fitted_12deg_warm = F2_temp_mod_eco$GM$fitted[,,206]
F2_fitted_12deg_warm_mat = as.matrix(F2_fitted_12deg_warm)
F2_12deg_warm_array = array(F2_fitted_12deg_warm_mat, 
                            dim = c(37, 2, 1))

F2_fitted_18deg_cold = F2_temp_mod_eco$GM$fitted[,,105]
F2_fitted_18deg_cold_mat = as.matrix(F2_fitted_18deg_cold)
F2_18deg_cold_array = array(F2_fitted_18deg_cold_mat, 
                            dim = c(37, 2, 1))

F2_fitted_18deg_warm = F2_temp_mod_eco$GM$fitted[,,310]
F2_fitted_18deg_warm_mat = as.matrix(F2_fitted_18deg_warm)
F2_18deg_warm_array = array(F2_fitted_18deg_warm_mat, 
                            dim = c(37, 2, 1))

F2_12deg_cold_range = c(1:50, 105:155, 412:461, 512:561, 798:847, 898:947, 1198:1247, 1298:1347)
F2_12deg_warm_range = c(206:259, 310:359, 612:650, 701:746, 998:1047, 1098:1147, 1398:1447, 1475:1524)
F2_18deg_cold_range = c(51:104, 156:205, 462:511, 562:611, 848:897, 948:997, 1248:1297, 1348:1397)
F2_18deg_warm_range = c(260:309, 360:411, 651:700, 747:797, 1048:1097, 1148:1197, 1448:1474, 1525:1575)

F2_temp_eco_array = array(0,
                          dim = c(37, 2, 1575))

for(i in F2_12deg_cold_range){
  F2_temp_eco_array[,,i] = gpa$coords[,,i] - F2_12deg_cold_array[,,1]
}

for(i in F2_12deg_warm_range){
  F2_temp_eco_array[,,i] = gpa$coords[,,i] - F2_12deg_warm_array[,,1]
}

for(i in F2_18deg_cold_range){
  F2_temp_eco_array[,,i] = gpa$coords[,,i] - F2_18deg_cold_array[,,1]
}

for(i in F2_18deg_warm_range){
  F2_temp_eco_array[,,i] = gpa$coords[,,i] - F2_18deg_warm_array[,,1]
}

Eco_array_consensus = array(0, dim = c(37, 2, 1575))

for(i in 1:1575){
  Eco_array_consensus[,,i] = F2_temp_eco_array[,,i] + mean_shape_array[,,1]
}

writeland.tps(Eco_array_consensus,
              file = 'Ecotype_variation_landmarks.tps',
              scale = NULL,
              specID = T)




# isolate ecotype effects -------------------------------------------------

## isolate effects due to warm cold divergence axis regardless of population
ecotype_mod1 = procD.lm(gpa$coords ~ meta_data$ecotype, 
                       iter = 999, 
                       RRPP = T)

eco1_fitted_cold = ecotype_mod1$GM$fitted[,,1]
eco1_fitted_mat_cold = as.matrix(eco1_fitted_cold)
eco1_cold_array = array(eco1_fitted_mat_cold, 
                       dim = c(37, 2, 1))

eco1_fitted_warm = ecotype_mod1$GM$fitted[,,292]
eco1_fitteed_warm_mat = as.matrix(eco1_fitted_warm)
eco1_warm_array = array(eco1_fitteed_warm_mat, 
                       dim = c(37, 2, 1))

eco1_cold_range = c(1:204, 412:611, 798:997, 1198:1397)

eco1_warm_range = c(205:411, 612:797, 998:1197, 1397:1575)

eco1_array = array(0, dim = c(37, 2, 1575))

for(i in eco1_cold_range){
  eco1_array[,,i] = gpa$coords[,,i] - eco1_cold_array[,,1]
}

for(i in eco1_warm_range){
  eco1_array[,,i] = gpa$coords[,,i] - eco1_warm_array[,,1]
}


eco1_array_consensus = array(0, dim = c(37, 2, 1575))

for(i in 1:1575){
  eco1_array_consensus[,,i] = eco1_array[,,i] + mean_shape_array[,,1]
}

writeland.tps(eco1_array_consensus,
              file = 'Ecotype_effect_landmarks_all_individuals.tps',
              scale = NULL,
              specID = T)

# ecotype_mod2 = procD.lm(gpa$coords ~ meta_data$poppair*meta_data$ecotype...12, 
#                         iter = 999, 
#                         RRPP = T)
# 
# ecotype2_fitted = ecotype_mod2$GM$fitted
# 
# ecotype2_consensus = array(0, dim = c(37, 2, 1575))
# 
# for(i in 1:1575){
#   ecotype2_consensus[,,i] = ecotype2_fitted[,,i] + mean_shape_array[,,1]
# }
# 
# writeland.tps(ecotype2_consensus,
#               file = 'ecotype_effect_per_population_landmarks_all_individuals.tps',
#               scale = NULL,
#               specID = T)
# 
# 


# Phenotype PCA -----------------------------------------------------------

shape_data = readmulti.tps(c('F2_All_aligned_withsliders.tps',
                             'F1_effect_landmarks_all_individuals.tps',
                             'F2_effect_landmarks_all_individuals.tps',
                             'Ecotype_effect_landmarks_all_individuals.tps'),
                           specID = 'imageID')

sliders = define.sliders(c(28:37,1))

shape_gpa = gpagen(shape_data,
                   print.progress = T,
                   curves = sliders)

id = read_csv('shape_data_id.csv')

coord_sub = coords.subset(shape_gpa$coords,
                          id$shape_data)

raw_F2_coords = coord_sub$raw
TGP_coords = coord_sub$TGP
WGP_coods = coord_sub$WGP
Eco_coords = coord_sub$Eco

## pca of the raw landmark data
# raw = readland.tps('F2_All_aligned_withsliders.tps',
#                    specID = 'imageID',
#                    readcurves = T)
# 
# 
# raw_gpa = gpagen(raw,
#              print.progress = T,
#              curves = sliders)

raw_pca = gm.prcomp(A = raw_F2_coords)
# raw_pca = gm.prcomp(A = raw_gpa$coords)

summary(raw_pca)

raw_pca_dim = bsDimension(raw_pca$x)
raw_pca_scree = screeplot(raw_pca)

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

# raw_pca_data %>%
#   write_csv('Common_GPA_Unfilered_PCA_data.csv')



## pca of the f1 effects
# F1_effects = readland.tps('F1_effect_landmarks_all_individuals.tps', 
#                           specID = 'imageID', 
#                           readcurves = T)
# 
# f1_gpa = gpagen(F1_effects, 
#                 curves = sliders)
f1_pca = gm.prcomp(TGP_coords)
summary(f1_pca)


TGP_dim = bsDimension(f1_pca$x)
TGP_pca_scree = screeplot(f1_pca)

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

# f1_pca_data %>% 
#   write_csv('F1_effect_PCA_data.csv')

## PCA of the TGP by ecotype effects
TGP_eco = readland.tps('TGP_ecotype_variation_landmarks.tps', 
                          specID = 'imageID', 
                          readcurves = T)

TGP_eco_gpa = gpagen(TGP_eco, 
                curves = sliders)
TGP_eco_pca = gm.prcomp(TGP_eco_gpa$coords)

summary(TGP_eco_pca)

bsDimension(TGP_eco_pca$x)

TGP_eco_pca_vals = TGP_eco_pca$x %>% 
  as_tibble() %>% 
  select(1:5)

TGP_eco_pca_data = bind_cols(meta_data, 
                        TGP_eco_pca_vals)

TGP_eco_pca_data$F1 = as.character(TGP_eco_pca_data$F1)
TGP_eco_pca_data$F2 = as.character(TGP_eco_pca_data$F2)

ggplot(data = TGP_eco_pca_data)+
  geom_point(aes(x = Comp1, 
                 y = Comp2, 
                 col = F1, 
                 shape = ecotype))

TGP_eco_pca_data %>%
  write_csv('TGP_ecotype_variation_PCA_data.csv')


## pca of the f2 effects
f2_effects = readland.tps('F2_effect_landmarks_all_individuals.tps', 
                          specID = 'imageID', 
                          readcurves = T)

f2_gpa = gpagen(f2_effects, 
                curves = sliders)
f2_pca = gm.prcomp(f2_gpa$coords)

summary(f2_pca)
bsDimension(f2_pca$x)

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

# f2_pca_data %>% 
#   write_csv('F2_effect_pca_data.csv')

## PCA of the WGP by ecotype effects
WGP_eco = readland.tps('WGP_ecotype_variation_landmarks.tps', 
                       specID = 'imageID', 
                       readcurves = T)

WGP_eco_gpa = gpagen(WGP_eco, 
                     curves = sliders)
WGP_eco_pca = gm.prcomp(WGP_eco_gpa$coords)

summary(WGP_eco_pca)

WGP_eco_pca_vals = WGP_eco_pca$x %>% 
  as_tibble() %>% 
  select(1:5)

WGP_eco_pca_data = bind_cols(meta_data, 
                             WGP_eco_pca_vals)

WGP_eco_pca_data$F1 = as.character(WGP_eco_pca_data$F1)
WGP_eco_pca_data$F2 = as.character(WGP_eco_pca_data$F2)

ggplot(data = WGP_eco_pca_data)+
  geom_point(aes(x = Comp1, 
                 y = Comp2, 
                 col = F2, 
                 shape = ecotype))


WGP_eco_pca_data %>%
  write_csv('WGP_ecotype_variation_PCA_data.csv')


## pca of the cold vs warm ecotype effects
eco1_effects = readland.tps('Ecotype_effect_landmarks_all_individuals.tps',
                            specID = 'imageID',
                            readcurves = T)


# eco1_effects = readland.tps('ecotype_effect_per_population_landmarks_all_individuals.tps', 
#                           specID = 'imageID', 
#                           readcurves = T)

eco1_gpa = gpagen(eco1_effects, 
                curves = sliders)
eco1_pca = gm.prcomp(eco1_gpa$coords)

summary(eco1_pca)

bsDimension(eco1_pca$x)

eco1_pca_vals = eco1_pca$x %>% 
  as_tibble() %>% 
  select(1:5)

eco1_pca_data = bind_cols(meta_data, 
                        eco1_pca_vals)

# eco1_pca_data$F1 = as.character(eco1_pca_data$F1)
# eco1_pca_data$F2 = as.character(eco1_pca_data$F2)

ggplot(data = eco1_pca_data)+
  geom_point(aes(x = Comp1, 
                 y = Comp2, 
                 col = ecotype))

# eco1_pca_data %>% 
#   write_csv('Ecotype_effect_pca_data.csv')


# interlandmark distances -------------------------------------------------

## raw body shape data


bs_data = readland.tps('F2_All_aligned_withsliders.tps', 
                       specID = 'imageID', 
                       readcurves = T)
raw_pca = read_csv('Unfilered_PCA_data.csv')

sliders = define.sliders(c(28:37,1))
gpa = gpagen(bs_data, 
             print.progress = T, 
             curves = sliders)

raw_lmks = data.frame(body_width = c(12, 21), 
                      body_length = c(1, 16), 
                      row.names = c('start', 
                                    'end'))

raw_lmk_dist = interlmkdist(gpa$coords, 
             raw_lmks)

raw_lmk_dist %>% 
  as_tibble() %>% 
  mutate(fineness_ratio = body_length/body_width) %>% 
  bind_cols(raw_pca, 
            .) %>% 
  write_csv('Raw_RDA_phenotypes.csv')

## F1 effects

F1_data = readland.tps('F1_effect_landmarks_all_individuals.tps', 
                       specID = 'imageID', 
                       readcurves = T)
F1_pca = read_csv('F1_effect_PCA_data.csv')

sliders = define.sliders(c(28:37,1))
f1_gpa = gpagen(F1_data, 
             print.progress = T, 
             curves = sliders)

lmks = data.frame(body_width = c(12, 21), 
                      body_length = c(1, 16), 
                      row.names = c('start', 
                                    'end'))

F1_lmk_dist = interlmkdist(f1_gpa$coords, 
                            lmks)

F1_lmk_dist %>% 
  as_tibble() %>% 
  mutate(fineness_ratio = body_length/body_width) %>% 
  bind_cols(F1_pca, 
            .) %>% 
  write_csv('F1_RDA_phenotypes.csv')


## F2 effects
F2_data = readland.tps('F2_effect_landmarks_all_individuals.tps', 
                       specID = 'imageID', 
                       readcurves = T)
F2_pca = read_csv('F2_effect_PCA_data.csv')

sliders = define.sliders(c(28:37,1))
F2_gpa = gpagen(F2_data, 
                print.progress = T, 
                curves = sliders)

lmks = data.frame(body_width = c(12, 21), 
                  body_length = c(1, 16), 
                  row.names = c('start', 
                                'end'))

F2_lmk_dist = interlmkdist(F2_gpa$coords, 
                           lmks)

F2_lmk_dist %>% 
  as_tibble() %>% 
  mutate(fineness_ratio = body_length/body_width) %>% 
  bind_cols(F2_pca, 
            .) %>% 
  write_csv('F2_RDA_phenotypes.csv')

## ecotype effects

Ecotype_data = readland.tps('Ecotype_effect_landmarks_all_individuals.tps', 
                       specID = 'imageID', 
                       readcurves = T)
Ecotype_pca = read_csv('Ecotype_effect_PCA_data.csv')

sliders = define.sliders(c(28:37,1))
Ecotype_gpa = gpagen(Ecotype_data, 
                print.progress = T, 
                curves = sliders)

lmks = data.frame(body_width = c(12, 21), 
                  body_length = c(1, 16), 
                  row.names = c('start', 
                                'end'))

Ecotype_lmk_dist = interlmkdist(Ecotype_gpa$coords, 
                           lmks)

Ecotype_lmk_dist %>% 
  as_tibble() %>% 
  mutate(fineness_ratio = body_length/body_width) %>% 
  bind_cols(Ecotype_pca, 
            .) %>% 
  write_csv('Ecotype_RDA_phenotypes.csv')
