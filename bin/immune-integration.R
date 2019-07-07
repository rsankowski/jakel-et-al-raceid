#from url https://satijalab.org/seurat/v3.0/immune_alignment.html

library(tidyverse)
library(Seurat)
library(Matrix)
library(RaceID)
library(cowplot)

date <- Sys.Date()


#load human data

              load('/home/roman/Documents/Single cell analysis/20180911-Healthy-microglia-for-resub/RaceID-Cl-14-fused/20181128-Healthy-human-micr-cl14-fused-sc-RaceID3.Robj')
              use_all_genes = T
              if (use_all_genes) {              
                #micr_data_cel <- Matrix(as.matrix(sc@ndata), sparse = T)
                micr_data_cel <- Matrix(as.matrix(sc@expdata), sparse = T)[,colnames(sc@ndata)]
                
              } else {
                micr_data_cel <- Matrix(as.matrix(sc@fdata), sparse = T)
              }
              micr_data_cel <- micr_data_cel[,!colnames(micr_data_cel) %in% names(sc@cpart)[sc@cpart %in% c(4,10:15) ]]
              rownames(micr_data_cel) <- gsub('_.*', '', rownames(micr_data_cel))
              
              #add batch info
              data_t <- data.frame(ID = colnames(sc@ndata), cell_ID = colnames(sc@ndata))
              
              #Add Cluster numbers
              data_t$Region_gm_wm <- ifelse(grepl('WM', data_t$ID), 'WM', 
                                            ifelse(grepl('GM', data_t$ID), 'GM', 'Both')) 
              
              data_t$ID <- gsub('_.*', '', data_t$ID)
              data_t$ID <- gsub('-.*', '', data_t$ID)
              
              data_t$ID <- gsub('GM|WM|all|micr|pos|17Pl1|17Pl2', '', data_t$ID)
              
              batch_info <- read.csv('/home/roman/Documents/Single cell analysis/20170921 Healthy microglia/20171210-batch-information.csv', stringsAsFactors = F, header = T)[,-1]
              data_t <- left_join(data_t, batch_info)
              data2 <- data_t[,c('cell_ID', 'Batch')]
              colnames(data2)[2] <- 'batch'
              table(data_t$ID, data_t$Batch)
              
              data_cel <- data.frame(cluster=as.character(sc@cpart), cell_ID=names(sc@cpart), dataset=rep('celseq', nrow(data_t))) %>% left_join(data2)

              
              micr_data_human <- micr_data_cel
              
              #metainfo
              data_human <- data.frame(cluster=as.character(sc@cpart), cell_ID=names(sc@cpart), dataset=rep('human', ncol(sc@ndata)))
              data_human <- data_human[!data_human$cluster %in% c(4,10:15), ]
              
              
              #ms dataset
              #load metadata
              url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE118257&format=file&file=GSE118257%5FMSCtr%5FsnRNA%5FFinalAnnotationTable%2Etxt%2Egz'
              tmp <- tempfile()
              ##
              download.file(url,tmp)
              metadata <- read.delim(gzfile(tmp), sep = '\t', stringsAsFactors = F)
              metadata$Detected <- gsub("\\:", "\\.", metadata$Detected)
              
              #load data
              if (!file.exists("data/jackel-et-al-micr.Robj")) {
              url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE118257&format=file&file=GSE118257%5FMSCtr%5FsnRNA%5FExpressionMatrix%5FR%2Etxt%2Egz'
              tmp <- tempfile()
              ##
              download.file(url,tmp)
              prdata <- read.delim(gzfile(tmp), sep = '\t', stringsAsFactors = F)
              #extract microglia
              micr_data_ms <- prdata[, colnames(prdata) %in% metadata$Detected[metadata$Celltypes == "Microglia_Macrophages"]]
              rownames(micr_data_ms) <- gsub("\\.", "\\-", rownames(micr_data_ms))
              
              } else {
                load("data/jackel-et-al-micr.Robj")
              }
              
              #regev data
              load("/home/roman/Documents/Single cell analysis/20171213-Healthy-vs-tumor-microglia/RaceID/20180416-harvard+stanford-tumor-microglia.Robj")
              
              tirosh <- paste(c("^(MGH36", "MGH53", "MGH54", "MGH60", "MGH93", "MGH97)"), collapse = "|")
              oligo_data <- x[,grepl(tirosh, colnames(x))]
              data_tirosh <- data.frame(cluster=NA, cell_ID=colnames(oligo_data), dataset=rep('tirosh-et-al', ncol(oligo_data)))
              
              #obtain common row names
              genes <- unique(rownames(micr_data_ms))
              #genes <- gsub("_.*", "", rownames(sc@fdata))
              genes2 <- genes[genes %in% rownames(micr_data_human)]
              genes2 <- genes2[genes2 %in% rownames(oligo_data)]
              ms.data <- micr_data_ms
              human.data <- micr_data_human[genes2,]
              human.data <- human.data #[,sample(1:ncol(human.data), 2000)]
              oligo.data <- oligo_data[genes2,]
              
# Set up control object
human <- CreateSeuratObject(counts = human.data, project = "IMMUNE_human", min.cells = 5)
human$stim <- "sankowski-et-al"
human <- subset(human, subset = nFeature_RNA > 450)
human <- NormalizeData(human, verbose = FALSE)
human <- FindVariableFeatures(human, selection.method = "dispersion", nfeatures = 2000)

# Set up stimulated object
ms <- CreateSeuratObject(counts = ms.data, project = "IMMUNE_ms", min.cells = 5)
ms$stim <- "jackel-et-al"
ms <- subset(ms, subset = nFeature_RNA > 450)
ms <- NormalizeData(ms, verbose = FALSE)
ms <- FindVariableFeatures(ms, selection.method = "dispersion", nfeatures = 2000)

# Set up stimulated object
oligo <- CreateSeuratObject(counts = oligo.data, project = "IMMUNE_oligo", min.cells = 5)
oligo$stim <- "tirosh-et-al"
oligo <- subset(oligo, subset = nFeature_RNA > 450)
oligo <- NormalizeData(oligo, verbose = FALSE)
oligo <- FindVariableFeatures(oligo, selection.method = "dispersion", nfeatures = 2000) #vst


#integrate datasets
immune.anchors <- FindIntegrationAnchors(object.list = list(human, ms, oligo), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

#integrated analysis
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

#plot
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

save(immune.combined, file = "data/immune-dataset-integration-control-ms-oligodendroglioma.Robj")






#integration with less clusters

#integrate datasets
immune.anchors <- FindIntegrationAnchors(object.list = list(human, ms, oligo), dims = 1:30)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:30)

#integrated analysis
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.9)

#plot
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

   