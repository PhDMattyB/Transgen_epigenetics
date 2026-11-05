### Methylation topGO 

library(tidyverse)
library(topGO)

setwd("C:/Users/phdma/OneDrive/OneDrive - University of Glasgow/Documents/Parsons_Postdoc/Methylation_data/GLM_results")

topDiffGenes <- function(allScore) {
  return(allScore == '1') # adjust based on your p-value cutoff
}


stick_GO_df = read_tsv('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/GCF_964276395.1-RS_2025_12_gene_ontology.gaf')


# GLMM F1*Eco topGO -------------------------------------------------------

F1_eco_genes = read_tsv("GLMM_res_F1_Ecotype_inter_gene_names.txt") %>% 
  dplyr::rename(Gene = gene_name)

F1_eco_genelist = F1_eco_genes$Gene

geneID2GO = readMappings(file = "~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/Stickleback_Gene_GO_Map.txt") 
all_genes_universe = names(geneID2GO)

all_genelist = rep(0, length(all_genes_universe))
names(all_genelist) = all_genes_universe
all_genelist[names(all_genelist) %in% F1_eco_genelist] <- 1
F1_eco_genelist <- as.factor(all_genelist)

F1_eco_GOdata_BP <- new("topGOdata",
                        ontology = 'BP',              # Choose "BP", "MF", or "CC"
                        allGenes = F1_eco_genelist,          # Named vector of genes
                        geneSelectionFun = topDiffGenes,
                        nodeSize = 5,                # Minimum number of annotated genes per GO term (e.g., >= 10)
                        annot = annFUN.gene2GO,       # Use custom gene-to-GO mapping
                        gene2GO = geneID2GO)

F1_eco_GOdata_MF <- new("topGOdata",
                        ontology = 'MF',              # Choose "BP", "MF", or "CC"
                        allGenes = F1_eco_genelist,          # Named vector of genes
                        geneSelectionFun = topDiffGenes,
                        nodeSize = 5,                # Minimum number of annotated genes per GO term (e.g., >= 10)
                        annot = annFUN.gene2GO,       # Use custom gene-to-GO mapping
                        gene2GO = geneID2GO)

F1_eco_GOdata_CC <- new("topGOdata",
                        ontology = 'CC',              # Choose "BP", "MF", or "CC"
                        allGenes = F1_eco_genelist,          # Named vector of genes
                        geneSelectionFun = topDiffGenes,
                        nodeSize = 5,                # Minimum number of annotated genes per GO term (e.g., >= 10)
                        annot = annFUN.gene2GO,       # Use custom gene-to-GO mapping
                        gene2GO = geneID2GO)

F1_eco_resultFisher_BP <- runTest(
  F1_eco_GOdata_BP,
  algorithm = "weight01",
  statistic = "fisher")

F1_eco_resultFisher_MF <- runTest(
  F1_eco_GOdata_MF,
  algorithm = "weight01",
  statistic = "fisher")

F1_eco_resultFisher_CC <- runTest(
  F1_eco_GOdata_CC,
  algorithm = "weight01",
  statistic = "fisher")

F1_eco_GO_BP <- GenTable(
  F1_eco_GOdata_BP,
  weight01Fisher = F1_eco_resultFisher_BP,
  orderBy = "weight01Fisher",
  ranksOf = "weight01Fisher",
  topNodes = length(score(F1_eco_resultFisher_BP)))

# Molecular Function
F1_eco_GO_MF <- GenTable(
  F1_eco_GOdata_MF,
  weight01Fisher = F1_eco_resultFisher_MF,
  orderBy = "weight01Fisher",
  ranksOf = "weight01Fisher",
  topNodes = length(score(F1_eco_resultFisher_MF)))

# Cellular Component
F1_eco_GO_CC <- GenTable(
  F1_eco_GOdata_CC,
  weight01Fisher = F1_eco_resultFisher_CC,
  orderBy = "weight01Fisher",
  ranksOf = "weight01Fisher",
  topNodes = length(score(F1_eco_resultFisher_CC)))


F1_eco_GO_BP = F1_eco_GO_BP %>% 
  as_tibble() %>% 
  # mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
  #         FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(weight01Fisher <= 0.05)

F1_eco_GO_MF = F1_eco_GO_MF %>% 
  as_tibble() %>% 
  # mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
  #         FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(weight01Fisher <= 0.05) 

F1_eco_GO_CC = F1_eco_GO_CC %>% 
  as_tibble() %>% 
  # mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
  #         FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(weight01Fisher <= 0.05)

# F1_eco_BP_GO_genes = genesInTerm(
#   F1_eco_GOdata_BP,
#   F1_eco_GO_BP$GO.ID
# )
# 
# F1_eco_BP_GO_genes_sig <- lapply(
#   F1_eco_BP_GO_genes,
#   function(genes) dplyr::intersect(genes, F1_eco_genes$Gene))
# 
# 
# F1_eco_MF_GO_genes = genesInTerm(
#   F1_eco_GOdata_MF,
#   F1_eco_GO_MF$GO.ID
# )
# 
# F1_eco_MF_GO_genes_sig <- lapply(
#   F1_eco_MF_GO_genes,
#   function(genes) dplyr::intersect(genes, F1_eco_genes$Gene))
# 
# 
# F1_eco_CC_GO_genes = genesInTerm(
#   F1_eco_GOdata_CC,
#   F1_eco_GO_CC$GO.ID
# )
# 
# F1_eco_CC_GO_genes_sig <- lapply(
#   F1_eco_CC_GO_genes,
#   function(genes) dplyr::intersect(genes, F1_eco_genes$Gene))

