## Data integration of mRNA quantifications along with chromatin accesibility on the same cells

## Dependencies
library(Seurat)

#############################################################################################

## Loading data into R

## Loading Adult brain cortex mRNA quantified transcripts
cDNADir <- '../data/chenS2019sameCellsData/GSE126074_AdBrainCortex_SNAREseq_cDNA/'
adBrcDNA <- Read10X(data.dir = cDNADir, 
                    gene.column = 1      ## necessary parameter, otherwise runs error
                    )
adBrcDNA.Seu <- CreateSeuratObject(counts = adBrcDNA,
                                   project = 'sameCellsData', 
                                   min.cells = 10,
                                   min.features = 500)

## Loading Adult brain cortex chromatin accesibility reads
chrDir <- '../data/chenS2019sameCellsData/GSE126074_AdBrainCortex_SNAREseq_chromatin/'
adBrchr <- Read10X(data.dir = chrDir, 
                    gene.column = 1      ## necessary parameter, otherwise runs error
                    )
adBrchr.Seu <- CreateSeuratObject(counts = adBrcDNA,
                                   project = 'sameCellsData', 
                                   min.cells = 10,
                                   min.features = 500)


###########################################################################################

adBrcDNA.Seu[['percent.mt']] <- PercentageFeatureSet(adBrcDNA.Seu, pattern = '^MT-')
VlnPlot(adBrcDNA.Seu, 
        features = c('nFeature_RNA', 
                     'nCount_RNA'),
        ncol = 1)





