## Transform SNARE data to NMF Andres Q pipeline
library(dplyr)
library()

##############################################################################################
## Cells Mixture cDNA
mixcDNA <- read.table('../data/chenS2019sameCellsData/GSE126074_CellLineMixture_SNAREseq_cDNA_counts.tsv',
                      sep = '\t',
                      stringsAsFactors = FALSE) 
rnaseq <- mixcDNA[ , 1:50] %>% 
        as.matrix(mixcDNA) 
saveRDS(rnaseq, file = '../data/NMF/input/rnaseq.RDS')

## annotations
clines.ann <- read.table('../data/cell_lines_metadata.tsv',
                         sep = '\t',
                         stringsAsFactors = FALSE,
                         nrows = 50, 
                         header = TRUE)

assigned_cols <- plyr::mapvalues(clines.ann$Ident, 
                                 from = c('BJ', 'K562', 'GM12878', 'H1'),
                                 to = c('forestgreen', 'indianred', 'salmon', 'steelblue'))
clines.ann <- mutate(clines.ann, 
                     color = assigned_cols,
                     sampleID = Barcode,
                     Celltype = as.factor(setNames(nm = (Ident))),
                     rna.sampleID = Barcode,
                     atac.sampleID = Barcode,
                     original.atacID = Barcode)
str(clines.ann)
saveRDS(clines.ann, file = '../data/NMF/input/rna_annotation.RDS')

##############################################################################################
## Cells Mixture chromatin
mixchr <- read.table('../data/chenS2019sameCellsData/GSE126074_CellLineMixture_SNAREseq_chromatin_counts.tsv',
                     sep = '\t',
                     stringsAsFactors = FALSE) 
atacseq <- mixchr[ , 1:50] %>%
        as.matrix(mixchr)
saveRDS(atacseq, file = '../data/NMF/input/atacseq.RDS')

#############################################################################################
## rna and atac integration
multiview <- list(rnaseq, atacseq)

saveRDS(multiview, '../data/NMF/input/multiview.RDS')





