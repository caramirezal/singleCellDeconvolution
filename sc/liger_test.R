## liger test
library(liger)
library(dplyr)

## loading Tirosh et, 2016 data. GSE70630  
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102130
idh_oligo <- read.table('data/GSE70630_OG_processed_data_v2.txt')
dim(idh_oligo)

## loading Filbin et, 2018. GSE70630
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70630
h3k27m_glioma <- read.table('data/GSE102130_K27Mproject.RSEM.vh20170621.txt',
                            header = TRUE)
row.names(h3k27m_glioma) <- h3k27m_glioma$Gene
h3k27m_glioma <- select(h3k27m_glioma, -Gene)
dim(h3k27m_glioma) 

## preprocessing data
liger <- createLiger(list(idh=idh_oligo, h3k27=h3k27m_glioma))
liger <- normalize(liger)
liger <- selectGenes(liger, var.thresh = 0.1)
liger <- scaleNotCenter(liger)

## running liger
liger <- optimizeALS(liger, k = 20) 
liger <- quantileAlignSNF(liger) #SNF clustering and quantile alignment

## saving results
saveRDS(liger, 'data/liger_results.RDS')
