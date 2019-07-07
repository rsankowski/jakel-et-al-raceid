library(tidyverse)
library(viridis)
library(readxl)
library(RaceID)
library(RColorBrewer)
library(pheatmap)
source("~/Documents/Single cell analysis/RaceID3_StemID2/RaceID3_StemID2_class.R")


date = Sys.Date()

source("/home/roman/Documents/Single cell analysis/Advanced-plots/20181025-sankowski-et-al-functions.R")

load("data/immune-dataset-integration-control-ms-oligodendroglioma.Robj")

#load healthy dataset
              load('/home/roman/Documents/Single cell analysis/20180911-Healthy-microglia-for-resub/RaceID-Cl-14-fused/20181128-Healthy-human-micr-cl14-fused-sc-RaceID3.Robj')
              
              cell_numbers <-as.numeric()
              for (i in 1:length(unique(sc@cpart)))
              {
                cell_numbers[i] <- length(sc@cpart[sc@cpart==i])
              }
              names(cell_numbers) <- c(1:length(unique(sc@cpart)))
              retain_cl <- as.numeric(names(cell_numbers[cell_numbers > dim(sc@ndata)[2]/100]))
              retain_cl <- retain_cl[retain_cl != 4]
              order_clusters <- clustheatmap(sc, final = T)
              order_clusters <- order_clusters[order_clusters %in% retain_cl]
              
              #build a table with cell Ids, cluster and conditions
              df <- data.frame('Cluster' = sc@cpart, sc@tsne) 
              df$ID <- rownames(df)
              df$cell_ID <- rownames(df)
              df$Condition <- gsub('_.*', '', rownames(df))
              df$Patient_ID <- df$Condition
              df <- df[df$Cluster %in% retain_cl,] 
              order_clusters <- order_clusters[order_clusters %in% retain_cl]
              df$Cluster <- factor(df$Cluster, levels = order_clusters)
              #define region
              df$Region <- ifelse(grepl('WM',df$ID), 'WM', 
                                  ifelse(grepl('GM',df$ID), 'GM', 'Mixed'))
              
              df$Region <- factor(df$Region, levels = c('GM', 'WM', 'Mixed'))
              
              batch_info <- read.csv('/home/roman/Documents/Single cell analysis/20170921 Healthy microglia/20171210-batch-information.csv', stringsAsFactors = F, header = T)[,-1]
              colnames(batch_info)[c(1,4)] <- c('Condition', 'Diagnosis')
              df <- left_join(df, batch_info)
              table(df$Condition)
              
              df$Condition <- reorder(df$Condition, df$Age)
              df$anon_ID <- df$Condition
              levels(df$anon_ID) <- paste0('Pat', 1:15)
              df2 <- df
              df2$Patient <- rep("ctrl_micr", nrow(df2))
              colnames(df2)[1] <- "Cluster_orig"
              
#load metadata for MS dataset
              url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE118257&format=file&file=GSE118257%5FMSCtr%5FsnRNA%5FFinalAnnotationTable%2Etxt%2Egz'
              tmp <- tempfile()
              ##
              download.file(url,tmp)
              metadata <- read.delim(gzfile(tmp), sep = '\t', stringsAsFactors = F)[metadata$Celltypes == "Microglia_Macrophages",]
              metadata$Detected <- gsub("\\:", "\\.", metadata$Detected)
              colnames(metadata)[c(1,6,4)] <- c("cell_ID", "Cluster_orig", "Patient")
              metadata$Cluster_orig <- as.factor(paste0('X', metadata$Cluster_orig))

#prepare a dataframe of the merged dataset
              df <- data.frame('Cluster' = as.numeric(immune.combined$seurat_clusters), immune.combined@reductions$umap@cell.embeddings) 
              df$ID <- rownames(df)
              
              df$Condition <- ifelse(grepl("^MGH", df$ID), "tirosh-et-al", 
                                     ifelse(grepl("10X", df$ID), "jackel-et-al", "sankowski-et-al"))
              
              cell_numbers <- numeric()
              for (i in 1:max(df$Cluster, na.rm = T)) {
                cell_numbers[i] <- length(na.omit(df$Cluster[df$Cluster==i]))
              }
              names(cell_numbers) <- c(1:max(df$Cluster, na.rm = T))
              retain_cl <- as.numeric(names(cell_numbers[cell_numbers > dim(immune.combined@assays$integrated)[2]/100]))
              if (!file.exists('data/ord_clust_combined.Robj')) {
                clust <- list()
                for (i in unique(immune.combined$seurat_clusters)) {
                  clust[i] <- immune.combined@assays$integrated[,colnames(immune.combined@assays$integrated) %in% names(immune.combined$seurat_clusters)[immune.combined$seurat_clusters == i]]
                }
                
                clust2 <- clust %>% 
                  map(rowSums) %>% 
                  bind_rows(.id = 'Cluster')
                ord_clust <- hclust(dist(t(as.matrix(clust2)[,-1])))
                
                ord_clust <- unique(immune.combined$seurat_clusters)[ord_clust$order]
                save(ord_clust, file = 'data/ord_clust_combined.Robj')
              } else {
                load('data/ord_clust_combined.Robj')
              }
              ord_clust <- ord_clust[ord_clust %in% retain_cl]
              df$Cluster <- factor(df$Cluster, levels = ord_clust)
              
              df <- df[df$Cluster %in% retain_cl,]
              colnames(df)[2:4] <- c("V1", "V2","cell_ID")


#add metadata to df
              meta2 <- df2[,c("cell_ID","Cluster_orig", "Patient")] %>% 
                bind_rows(metadata[,c("cell_ID","Cluster_orig", "Patient", "Lesion")])
              
              df <- df %>%
                left_join(meta2) 
              
#Plot tsne map
              tsne <- tsne_plot(df, FILL = df$Condition, fill_colors = c(colors_pat), point_outline = "black", point_size = 3, line_width = 0.25) 
              
              tsne
              
              ggsave(paste0('plots/tsne/', date, '-condition-tsne-plot.pdf'), width = 8.57, height = 5.79)  
              
              svg(paste0('plots/tsne/', date, '-condition-tsne-plot.svg'), width = 8.57, height = 5.79)
              tsne
              dev.off()
              
              #marimekko plot - patients
              mosaicGG2(df, "Cluster", "Condition", c(colors_pat, colors_many), rect_col = 'black', line_width = 0.1) +
                scale_fill_brewer(palette = "Set3")
              ggsave(paste0('plots/others/', date,'-compartment-marimekko-cluster-plot.pdf'))
              
              #marimekko stat plot
              mosaicGG(df, "Cluster", "Condition", rect_col = 'black', line_width = 0.1)
              ggsave(paste0('plots/others/', date,'-conditions-marimekko-cluster-stat-plot.pdf'))
              

              #marimekko plot - cluster results
              mosaicGG2(df, "Cluster", "Cluster_orig", c(colors_pat, colors_many), rect_col = 'black', line_width = 0.1) +
                scale_fill_brewer(palette = "Set3")
              ggsave(paste0('plots/others/', date,'-compartment-marimekko-cluster-plot.pdf'))
              
              #marimekko stat plot
              mosaicGG(df, "Cluster", "Cluster_orig", rect_col = 'black', line_width = 0.1)
              ggsave(paste0('plots/others/', date,'-conditions-marimekko-cluster-stat-plot.pdf'))
              
              #original clustering
              tsne <- tsne_plot(df, FILL = df$Cluster_orig, fill_colors = c(colors_pat), point_outline = "black", point_size = 3, line_width = 0.25) 
              
              tsne
              