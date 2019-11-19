## Scripts containing the pipeline for analyzing single cell RNA seq data

## Input: A directory containing 10x cell Ranger files
## Load the data into a seurat object, performs normalization
## scaling, pca, umap and calculates feature variability
## Output: A seurat object 
seurat_pipe_10x <- function(
        data_dir='.',
        gene.column = 1,
        projectName = 'sc_data_analysis',
        min.cells = 1,
        min.features = 1,
        n_features = 1000,
        annotations = FALSE,
        pca_dims = 15) {
        cat('Loading data\n')
        counts_10x <- Read10X(data.dir = data_dir, 
                              ## necessary parameter, otherwise runs error
                              gene.column = gene.column)
        seurat_pipe <- CreateSeuratObject(counts = counts_10x,
                                          project = projectName, 
                                          min.cells = min.cells,
                                          min.features = min.features)
        
        cat('Normalizing and finding variable features\n')
        seurat_pipe <- NormalizeData(seurat_pipe) %>%
                FindVariableFeatures(selection.method = 'vst',
                                     nfeatures = n_features)
        
        cat('Adding annotations\n')
        if ( class(annotations) %in% c('factor', 'character') &
             length(annotations) == length(colnames(seurat_pipe))
        ) {
                seurat_pipe$'annotations' <- annotations
        }
        
        cat('Scaling and projection\n')
        seurat_pipe <- ScaleData(seurat_pipe, 
                                 verbose = FALSE) %>% 
                RunPCA(npcs = pca_dims, 
                       verbose = FALSE) %>%
                RunUMAP(reduction = "pca", 
                        dims = 1:15)
        
        return(seurat_pipe)
}
## Example
## setting pathway
#dir_10x <- '../data/chenS2019sameCellsData/GSE126074_AdBrainCortex_SNAREseq_cDNA/'
## Loading annotations
#metadata  <- readRDS("AdBrainCortex_SNAREseq_metadata.rds")
#adBr_mRNA <- seurat_pipe_10x(
#        dir_10x,
#        projectName = 'mRNA',
#        annotations = metadata$Ident
#)

## Input: A directory with files in 10x cell Ranger format
## Loads data, creates an activity matrix using peaks counts and
## sets up a seurate object
## Output: A seurat object
seurat_pipe_atac <- function(
        data_dir = '.',
        n_features = 3000,
        annotation_file = '.',
        seq_levels = c(1:19, 'X','Y')
){
        cat('Loading data \n')
        seurat_pipe_data <- Read10X(data.dir = data_dir,
                                    gene.column = 1)
        seurat_pipe <- CreateSeuratObject(
                counts = seurat_pipe_data,
                project = 'ad_chr',
                min.features = 1,
                min.cells = 1
        )
        
        cat('Feature selection of ', n_features,' most variables features\n')
        seurat_pipe <- FindVariableFeatures(seurat_pipe) %>% 
                subset(features = 1:n_features)
        
        cat('Calculating activity matrix from peaks counts\n')
        activity.matrix <- CreateGeneActivityMatrix(
                peak.matrix = GetAssayData(seurat_pipe), 
                annotation.file = annotation_file, 
                seq.levels = seq_levels, 
                upstream = 2000, 
                verbose = TRUE)
        
        cat('Setting the Seurat object containing scATAC counts\n')
        atac_seurat <- CreateSeuratObject(
                counts = GetAssayData(seurat_pipe, slot = 'counts'),
                assay = 'ATAC',
                project = 'atac'
        )
        
        cat('Adding estimated gene activity matrix\n')
        atac_seurat[['activity']] <- CreateAssayObject(counts = activity.matrix)
        
        return(atac_seurat)
}
## Loading Adult brain cortex chromatin accesibility reads
#chrDir <- '../data/chenS2019sameCellsData/GSE126074_AdBrainCortex_SNAREseq_chromatin/'
#ann_file_path <- "../data/Mus_musculus.GRCm38.98.gtf.gz"
#ad_br_chr <- seurat_pipe_atac(data_dir = chrDir,
#                              annotation_file = ann_file_path)

## Input: A seurat object
## Performs normalization, scaling, pca and umap 
## Output: A processed seurat object
st_workflow <- function(
        seurat_object,
        n_features = 3000,
        n_pca_dims = 15
){
        cat('Normalizing and finding variable features\n')
        seurat.p <- NormalizeData(seurat_object) %>%
                FindVariableFeatures(selection.method = 'vst',
                                     nfeatures = n_features)
        cat('Scaling and projection\n')
        seurat.p <- ScaleData(seurat.p, 
                              verbose = FALSE) %>% 
                RunPCA(npcs = n_pca_dims, 
                       verbose = FALSE) %>%
                RunUMAP(reduction = "pca", 
                        dims = 1:n_pca_dims) 
        return(seurat.p)
}
## Example:
#ad_br_chr <- st_workflow(ad_br_chr)