## Data integration of mRNA quantifications along with chromatin accesibility on the same cells

## Dependencies
library(Seurat)

## Loading data into R
dataDir <- '../data/chenS2019sameCellsData/GSE126074_AdBrainCortex_SNAREseq_cDNA/'
adBrcDNA <- Read10X(data.dir = dataDir, 
                    gene.column = 1      ## necessary parameter, otherwise runs error
                    )
adBrcDNA.Seu <- CreateSeuratObject(counts = adBrcDNA,
                                   project = 'sameCellsData', 
                                   min.cells = 10,
                                   min.features = 500)

VlnPlot(adBrcDNA.Seu, 
        features = c('nFeature_RNA', 
                     'nCount_RNA',
                     "percent.mt"),
        ncol = 1)




