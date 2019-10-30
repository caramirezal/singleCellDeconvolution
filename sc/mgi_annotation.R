## Annotation of mgi symbols
library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
riks <- readLines('../data/chenS2019sameCellsData/GSE126074_AdBrainCortex_SNAREseq_cDNA.genes.tsv')
riksAnn <- getLDS(attributes = c("mgi_symbol", 'ensembl_gene_id'),
                  filters = c("mgi_symbol"),
                  values = riks,
                  mart = mouse,
                  attributesL = 'hgnc_symbol', 
                  martL = human)

write.table(riksAnn, 
            file = '../data/mgi_ann.tsv', 
            sep = '\t',
            row.names = FALSE)