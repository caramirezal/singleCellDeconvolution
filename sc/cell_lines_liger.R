## liger test
library(liger)
library(dplyr)

## loading Tirosh et, 2016 data. GSE70630  
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102130
rna <- readRDS('data/cell_lines_rna.RDS')

## loading Filbin et, 2018. GSE70630
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70630
atac <- readRDS('data/cell_lines_activity_matrix.RDS')
colnames(atac) <- paste0(colnames(atac), '.atac')

## preprocessing data
liger <- createLiger(list(rna=rna, atac=atac))
liger <- normalize(liger)
liger <- selectGenes(liger, combine = 'union')
liger <- scaleNotCenter(liger)
 

## running liger
liger <- optimizeALS(liger, k = 4) 
liger <- quantileAlignSNF(liger) #SNF clustering and quantile alignment

## saving results
saveRDS(liger, 'data/cell_lines_liger.RDS')
