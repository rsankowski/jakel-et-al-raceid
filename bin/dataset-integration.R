#devtools::install_github(repo = "satijalab/seurat", ref = "release/3.0")
library(tidyverse)
library(Seurat)
library(Matrix)
library(RaceID)

date <- Sys.Date()

#load data
                #control data
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
                
                micr_data_human <- micr_data_cel
                
                #metainfo
                data_human <- data.frame(cluster=as.character(sc@cpart), cell_ID=names(sc@cpart), dataset=rep('sankowski-et-al', ncol(sc@ndata)))
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
                
                #metadata
                data_ms <- data.frame(cluster=rep("x6", ncol(micr_data_ms)), cell_ID=colnames(micr_data_ms), dataset=rep('jackel-et-al', ncol(micr_data_ms)))
                
                
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
                ms.data <- ms.data[genes2,]

#merge metadata
              metadata <- bind_rows(data_human, data_ms, data_tirosh)
              rownames(metadata) <- metadata$cell_ID

              #extract common rows from Olah et al dataset
              micr.data <- cbind(human.data, ms.data, oligo.data)

#save datasets
              save(metadata, file = 'data/metadata.Robj')
              save(micr.data, file = 'data/merged-counts.Robj')
              
micr <- CreateSeuratObject(counts = micr.data, meta.data = metadata)
micr.list <- SplitObject(object = micr, split.by = "dataset")
for (i in 1:length(x = micr.list)) {
  micr.list[[i]] <- NormalizeData(object = micr.list[[i]], verbose = FALSE)
  micr.list[[i]] <- FindVariableFeatures(object = micr.list[[i]], selection.method = "vst", 
                                             nfeatures = 20000, verbose = FALSE)
}

reference.list <- micr.list[c("sankowski-et-al", "jackel-et-al", "tirosh-et-al")]
micr.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30, anchor.features = 20000)
micr.integrated <- IntegrateData(anchorset = micr.anchors, dims = 1:30)

#save data
save(micr.integrated, file = "data/dataset-integration-ctrl-ms-oligo.Robj")
prdata <- micr.integrated@assays$integrated@data
save(prdata, file = "integrated-counts-ctrl-ms-oligo.Robj")

#plot
library(ggplot2)
library(cowplot)
# switch to integrated assay. The variable features of this assay are automatically set during
# IntegrateData
DefaultAssay(object = micr.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
micr.integrated <- ScaleData(object = micr.integrated, verbose = FALSE)
micr.integrated <- RunPCA(object = micr.integrated, npcs = 20, verbose = FALSE)
micr.integrated <- RunUMAP(object = micr.integrated, reduction = "pca", dims = 1:20)
p1 <- DimPlot(object = micr.integrated, reduction = "umap", group.by = "dataset")
p2 <- DimPlot(object = micr.integrated, reduction = "umap", group.by = "cluster", label = TRUE, 
              repel = TRUE) #+ NoLegend()
plot_grid(p1, p2)
ggsave(paste0('plots/', date, '-integration-mouse-human-dataset-idh-mut.png'))
ggsave(paste0('plots/', date, '-integration-mouse-human-dataset-idh-mut.pdf'))

svg(paste0('plots/', date, '-integration-mouse-human-dataset-idh-mut.svg'))
plot_grid(p1, p2)
dev.off()

#export cell embeddings
write_csv(data.frame(ID=rownames(micr.integrated@reductions$umap@cell.embeddings),micr.integrated@reductions$umap@cell.embeddings), paste0('data/', date, 'idh-mut-umap-cell-embeddings.csv'))



#
micr.query <- micr.list[["human"]]
micr.anchors <- FindTransferAnchors(reference = micr.integrated, query = micr.query, 
                                        dims = 1:30)
predictions <- TransferData(anchorset = micr.anchors, refdata = micr.integrated$celltype, 
                            dims = 1:30)
micr.query <- AddMetaData(object = micr.query, metadata = predictions)
micr.query$prediction.match <- micr.query$predicted.id == micr.query$celltype
table(micr.query$prediction.match)
table(micr.query$predicted.id)

