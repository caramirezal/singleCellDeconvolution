## Data integration of mRNA quantifications along with chromatin accesibility on the same cells

## Dependencies
library(Seurat)

#############################################################################################

## Loading data into R

#########################################################################################
## Adult brain cells

## Loading Adult brain cortex mRNA quantified transcripts
cDNADir <- '../data/chenS2019sameCellsData/GSE126074_AdBrainCortex_SNAREseq_cDNA/'
list.files(cDNADir)
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
## Embryonic cells 
## GSE126074_P0_BrainCortex
cDNADir <- '../data/chenS2019sameCellsData/GSE126074_P0_BrainCortex_SNAREseq_cDNA/'
list.files(cDNADir)
p0cDNA <- Read10X(data.dir = cDNADir, 
                            gene.column = 1      ## necessary parameter, otherwise runs error
        )
p0cDNA.Seu <- CreateSeuratObject(counts = p0cDNA,
                                 project = 'sameCellsData', 
                                 min.cells = 10,
                                 min.features = 500)


##############################################################################################
## Cells Mixture cDNA
mixcDNA <- read.table('../data/chenS2019sameCellsData/GSE126074_CellLineMixture_SNAREseq_cDNA_counts.tsv',
                      sep = '\t',
                      stringsAsFactors = FALSE) 
mixcDNA.seu <- CreateSeuratObject(counts = mixcDNA, 
                                  assay = 'mixCDNA')

##############################################################################################

## Cells Mixture chromatin
mixchr <- read.table('../data/chenS2019sameCellsData/GSE126074_CellLineMixture_SNAREseq_chromatin_counts.tsv',
                      sep = '\t',
                      stringsAsFactors = FALSE) 
mixchr.Seu <- CreateSeuratObject(counts = mixchr,
                                 assay = 'mixChr')

######################################################################################################################
#####################################################################################################################

## Preprocessing
mixchr.Seu <- NormalizeData(mixchr.Seu)
mixchr.Seu <- FindVariableFeatures(mixchr.Seu, 
                                   selection.method = 'vst',
                                   nfeatures = 30)
mixchr.anchors <- FindIntegrationAnchors(object.list = list(mixchr.Seu),
                                         dims = 1:3)

