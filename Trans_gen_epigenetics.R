##############################
## Transgenerational epigenetics
##
## Matt Brachmann (PhDMattyB)
##
## 07.10.2024
##
##############################


setwd('C:/Users/phdma/OneDrive/OneDrive - University of Glasgow/Documents/Parsons_Postdoc/Stickleback_Genomic/Transgen_epigenetics/')

library(tidyverse)

methy = read_rds(file = "Methylated.cov.noZero.rds")

unmethy = read_rds(file = 'UnMethylated.cov.noZero.rds')
