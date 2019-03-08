library(tidyverse)
library(viridis)
library(readxl)
library(RaceID)
library(RColorBrewer)
library(pheatmap)

date = Sys.Date()
load('data/sc.Robj')
source("/home/roman/Documents/Single cell analysis/Advanced-plots/20181025-sankowski-et-al-functions.R")

#load metadata
url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE118257&format=file&file=GSE118257%5FMSCtr%5FsnRNA%5FFinalAnnotationTable%2Etxt%2Egz'
tmp <- tempfile()
##
download.file(url,tmp)
metadata <- read.delim(gzfile(tmp), sep = '\t', stringsAsFactors = F)
colnames(metadata)[1] <- 'ID'

#build data frame
df <- data.frame('Cluster' = sc@cpart, sc@tsne) 
df$ID <- gsub('\\.', ':', rownames(df))

df <- df %>% 
  inner_join(metadata)

df$Condition <- factor(df$Condition, levels = c('Ctrl', 'MS'))

df$Celltypes <- factor(df$Celltypes, levels = c('OPCs', 'COPs', 'ImOlGs', 'Oligo1', 'Oligo2', 'Oligo3', 'Oligo4', 'Oligo5', 'Oligo6', 'Astrocytes', 'Astrocytes2', 'Neuron1','Neuron2','Neuron3', 'Neuron4', 'Neuron5', 'Endothelial_cells1', 'Endothelial_cells2', 'Macrophages', 'Microglia_Macrophages', 'Immune_cells', 'Vasc_smooth_muscle','Pericytes'))

df$Lesion <- factor(df$Lesion, levels = c('Ctrl', 'NAWM', 'A', 'CA', 'CI', 'RM'))

cell_numbers <- numeric()
for (i in 1:max(sc@cpart, na.rm = T)) {
  cell_numbers[i] <- length(na.omit(sc@cpart[sc@cpart==i]))
}
names(cell_numbers) <- c(1:max(sc@cpart, na.rm = T))
retain_cl <- as.numeric(names(cell_numbers[cell_numbers > dim(sc@ndata)[2]/100]))
load('ord_clust.Robj')
ord_clust <- ord_clust[ord_clust %in% retain_cl]
df$Cluster <- factor(df$Cluster, levels = ord_clust)

df <- df[df$Cluster %in% retain_cl,]

#export means per celltype gene expression table
colnames(sc@ndata) <- gsub('\\.', ':', colnames(sc@ndata))
ndata_summary <- data.frame(row.names = rownames(sc@ndata))
for (i in unique(metadata$Celltypes)) {
  step1 <- data.frame(rowMeans(as.matrix(sc@ndata)[,colnames(sc@ndata) %in% metadata$ID[metadata$Celltypes==i]]))
  colnames(step1) <- i
  ndata_summary <- cbind(ndata_summary, step1)
}
write.csv(ndata_summary, paste0('data/',date,'normalized-gene-counts-per-celltype.csv'))


#condition tsne map
tsne <- tsne_plot(df, FILL = df$Condition, fill_colors = c(colors_pat), point_outline = "black", point_size = 2.5, line_width = 0.05) +
  scale_fill_brewer(palette = "Set1")

tsne

ggsave(paste0('plots/tsne/', date, '-Conditions-tsne-plot.pdf'))  

svg(paste0('plots/tsne/', date, '-Conditions-tsne-plot.svg'), width = 8.57, height = 5.79)
tsne
dev.off()

#Plot tsne map - Celltypes
tsne <- tsne_plot(df[!df$Celltypes %in% c('Macrophages', 'Vasc_smooth_muscle'),], FILL = df$Celltypes[!data$Celltypes %in% c('Macrophages', 'Vasc_smooth_muscle')], fill_colors = c(colors_many, colors_pat), point_outline = "black", point_size = 2.5, line_width = 0.05) +
  theme(legend.position = 'None')

tsne

ggsave(paste0('plots/tsne/', date, '-stages-tsne-plot.pdf'))  

svg(paste0('plots/tsne/', date, '-stages-tsne-plot.svg'), width = 8.57, height = 5.79)
tsne
dev.off()

#tsne with legend
tsne <- tsne_plot(df[!df$Celltypes %in% c('Macrophages', 'Vasc_smooth_muscle'),], FILL = df$Celltypes[!data$Celltypes %in% c('Macrophages', 'Vasc_smooth_muscle')], fill_colors = c(colors_many, colors_pat), point_outline = "black", point_size = 2.5, line_width = 0.05)

tsne

ggsave(paste0('plots/tsne/', date, '-stages-tsne-plot-with-legend.pdf'))  

svg(paste0('plots/tsne/', date, '-stages-tsne-plot-with-legend.svg'), width = 8.57, height = 5.79)
tsne
dev.off()

#clusters
tsne <- tsne_plot(df, FILL = df$Cluster, fill_colors = c(colors_pat, colors_many), point_outline = "black", point_size = 2.5, line_width = 0.05) 

tsne


ggsave(paste0('plots/tsne/', date, '-clusters-tsne-plot.pdf'))  

svg(paste0('plots/tsne/', date, '-clusters-tsne-plot.svg'), width = 8.57, height = 5.79)
tsne
dev.off()


#marimekko plot - Conditions
mosaicGG2(df, "Cluster", "Condition", c(colors_pat, colors_many), rect_col = 'black', line_width = 0.1) +
  scale_fill_brewer(palette = "Set3")
ggsave(paste0('plots/tsne/', date,'-conditions-marimekko-cluster-plot.pdf'))

#marimekko stat plot
mosaicGG(df, "Cluster", "Condition", rect_col = 'black', line_width = 0.1)
ggsave(paste0('plots/tsne/', date,'-genotypes-marimekko-cluster-stat-plot.pdf'))

#Celltypes
mosaicGG2(df, "Cluster", "Celltypes", c(colors_many, colors_pat), rect_col = 'black', line_width = 0.1) 
ggsave(paste0('plots/tsne/', date,'-stages-marimekko-cluster-plot.pdf'))

#marimekko stat plot
mosaicGG(df, "Cluster", "Celltypes", rect_col = 'black', line_width = 0.1)
ggsave(paste0('plots/tsne/', date,'-stages-marimekko-cluster-stat-plot.pdf'))

#violin plot of gene expressions
igf1r <- data.frame('IGF1R' = as.matrix(sc@ndata)['IGF1R', ]* min(sc@counts))
rownames(igf1r) <- gsub('\\.', ':', rownames(igf1r))
data <- cbind(df, igf1r[df$ID,])
colnames(data)[12] <- 'IGF1R' 


sorted <- data %>% group_by(Lesion, Celltypes) %>%
  summarise(expression = mean(IGF1R)) %>%
  arrange(expression)
sorted2 <- data %>% group_by(Celltypes) %>%
  summarise(expression = mean(IGF1R)) %>%
  arrange(expression)


data$Celltypes <- factor(data$Celltypes, levels = as.character(sorted2$Celltypes))

ggplot(data[!data$Celltypes %in% c('Macrophages', 'Vasc_smooth_muscle'),], aes(Celltypes, IGF1R, fill=Celltypes)) +
  stat_summary(fun.y = mean, geom = "bar", color = 'black', lwd=0.25) + 
  stat_summary(fun.data = mean_se, geom = "errorbar",width = 0) +
  coord_flip() +
  theme_minimal() +
  labs(y='IGF1R expression') +
  #facet_wrap(~Lesion, scales = 'free_x')#+
  #scale_fill_manual(values = rev(colors_many))

ggsave(paste0('plots/others/',date, '-celltype-dependent-igf1r-expression.pdf'))


ggplot(data[!data$Celltypes %in% c('Macrophages', 'Vasc_smooth_muscle'),], aes(Celltypes, IGF1R, fill=Celltypes)) +
  stat_summary(fun.y = mean, geom = "bar", color = 'black', lwd=0.25) + 
  stat_summary(fun.data = mean_se, geom = "errorbar",width = 0) +
  coord_flip() +
  theme_minimal() +
  labs(y='IGF1R expression') +
  facet_wrap(~Lesion, scales = 'free_x')#+
#scale_fill_manual(values = rev(colors_many))

ggsave(paste0('plots/others/',date, '-celltype-and-lesion-dependent-igf1r-expression.pdf'))


#pdgfra
PDGFRA <- data.frame('PDGFRA' = as.matrix(sc@ndata)['PDGFRA', ]* min(sc@counts))
rownames(PDGFRA) <- gsub('\\.', ':', rownames(PDGFRA))
data <- cbind(df, PDGFRA[df$ID,])
colnames(data)[11] <- 'PDGFRA' 

data$Celltypes <- factor(data$Celltypes, levels = as.character(sorted$Celltypes))

ggplot(data[!data$Celltypes %in% c('Macrophages', 'Vasc_smooth_muscle'),], aes(Celltypes, PDGFRA, fill=Celltypes)) +
  #geom_violin(scale='width') +
  #geom_boxplot()+
  stat_summary(fun.y = mean, geom = "bar", color = 'black', lwd=0.25) + 
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar",width = 0) +
  coord_flip() +
  theme_minimal() +
  labs(y='PDGFRA expression') #+
#scale_fill_manual(values = rev(colors_many))

ggsave(paste0('plots/others/',date, '-celltype-dependent-PDGFRA-expression.pdf'))

#highlight lesions
fplt <- ggplot(df, aes(V1, V2, fill = Lesion)) +
  geom_point(pch=21, size=6, stroke=0.1) +
  scale_fill_brewer(palette = 'Dark2') +
  facet_wrap(~Lesion, nrow = 1) +
  theme_void() +
  theme(legend.position = 'None',
        text=element_text(size=20))

svg('plots/tsne/lesion-types.svg', width = 30, height = 3.3)
fplt
dev.off()

#marimekko stat plot
mosaicGG2(df[!df$Celltypes %in% c('Macrophages', 'Vasc_smooth_muscle'),], "Celltypes", "Lesion", c(colors_many, colors_pat), rect_col = 'black', line_width = 0.1) +
  scale_fill_brewer(palette = 'Dark2') 
ggsave(paste0('plots/others/', date,'-stages-marimekko-cluster-plot.pdf'))

#normalized cell counts
cells_table <- df %>% group_by(Lesion, Celltypes) %>%
  summarize(Freq = n()) %>%
  dplyr::mutate(Ratio = Freq/sum(Freq))

df$sum_cells_bin <- length(cells_table$Celltypes)

cells_table_wide <- spread(cells_table, Celltypes, Freq, fill = 0)
colSums(cells_table_wide[,3:25])

#normalize the counts from url:https://rpubs.com/bbolker/sweep_divide
cells_table_wide2 <- sweep(cells_table_wide[,3:25],MARGIN=2,FUN="/",STATS=colSums(cells_table_wide[,3:25]))
cells_table_wide2 <- (sweep(cells_table_wide2, MARGIN = 2, FUN ="*", STATS=min(colSums(cells_table_wide[,3:25])))) 

cells_table_wide[,3:25] <- cells_table_wide2
colSums(cells_table_wide[,3:25])

cells_table_scale <- gather(cells_table_wide, Celltypes, normalized_counts, OPCs:Pericytes)
#cells_table_scale$Lesion <- factor(cells_table_scale$Lesion, levels = c('Ctrl', 'NAWM', 'A', 'CA', 'CI', 'RM'))
cells_table <- left_join(cells_table, cells_table_scale)


cells_table <- cells_table %>% 
  group_by(Lesion) %>%
  dplyr::mutate(Freq = sum(normalized_counts)) %>%
  group_by(Celltypes) %>%
  dplyr::mutate(normalized_Ratio = (normalized_counts/Freq))

cells_table$Celltypes <- factor(cells_table$Celltypes, levels = c('OPCs', 'COPs', 'ImOlGs', 'Oligo1', 'Oligo2', 'Oligo3', 'Oligo4', 'Oligo5', 'Oligo6', 'Astrocytes', 'Astrocytes2', 'Neuron1','Neuron2','Neuron3', 'Neuron4', 'Neuron5', 'Endothelial_cells1', 'Endothelial_cells2', 'Macrophages', 'Microglia_Macrophages', 'Immune_cells', 'Vasc_smooth_muscle','Pericytes'))

#plot normalized ratios
cluster_age_stack_plot <- ggplot(cells_table[!cells_table$Celltypes %in% c('Macrophages', 'Vasc_smooth_muscle'),], aes(Lesion, normalized_Ratio, fill = Celltypes)) + #, levels = c(4, 12, 7, 6, 8, 5, 3, 1, 9, 10, 2, 15, 19)
  geom_bar(position = 'fill', stat = 'identity', color = 'black', lwd = 0.25) +
  #geom_line(size = 2) +
  theme_minimal() +
  #scale_x_continuous(limits = c(0.5,19.5), breaks = c(1:19)) +
  theme(panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        text=element_text(size=17)) +
  scale_fill_manual(values = c(colors_many, colors_pat)) + 
  labs(title = 'Cell contribution to Lesions', y='#cells in Lesion / total #cells', x='Lesion') +
  scale_y_reverse() 

cluster_age_stack_plot

ggsave(paste0('plots/others/', date,'-normalized-stack-plot-celltypes-per-lesion.pdf'))

#plot cell signatures
signature_genes <- data.frame("monocytes"=c('CCR2', 'CLEC12A', 'PLAC8', 'FCN1', 'S100A9'),
                              "microglia"= c('P2RY12', 'CX3CR1', 'CSF1R', 'TMEM119', 'SLC2A5'),
                              "tcell"=c('TRAC', 'TRBC2', 'CD52', 'IL32', NA),
                              "myeloid"=c('ITGAM',  'MS4A6A', 'TYROBP', 'CD14', NA),
                              "oligodendrocyte"=c('MBP',  'MOG', 'MAG', 'PLP1', NA),
                              "bcells"=c('CD79A', 'IGHG4', 'IGLL5', NA, NA))

source('~/Documents/Single cell analysis/Advanced-plots/20190102_plot_expmap.R')

for (i in colnames(signature_genes)) {
  tryCatch({
    svg(paste0('plots/tsne/', i, '-gene_signature.svg'), width = 8.57, height = 5.79)
    pl <- plot_expmap(gene=c(na.omit(as.character(signature_genes[[i]]))), point_size = 2.5)
    print(pl)
    dev.off()
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  #on.exit(dev.off())
}     


#plot single cell gene expression
load_data <- function(path) { 
  files <- dir(path, pattern = '\\.csv', full.names = TRUE)
  tables <- lapply(files, read.csv)
  do.call(rbind, tables)
}

up_genes <- load_data(file.path('data/Cluster specific genes/Up'))


source('~/Documents/Single cell analysis/Advanced-plots/20190102_plot_expmap.R')

plot_genes <- unique(as.character(up_genes$GENEID))
#plot_genes <- 'PDGFRA'
  
for (i in plot_genes) {
  tryCatch({
    svg(paste0('plots/tsne/',date,'-', i, '.svg'), width = 8.57, height = 5.79)
    pl <- plot_expmap(gene=c(i), point_size = 2.5)
    print(pl)
    dev.off()
    
    pdf(paste0('plots/tsne/',date,'-', i, '.pdf'), width = 8.57, height = 5.79)
    print(pl)
    dev.off()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  #on.exit(dev.off())
}     

#heatmaps

load_data <- function(path) { 
  files <- dir(path, pattern = '\\.csv', full.names = TRUE)
  tables <- lapply(files, read.csv)
  do.call(rbind, tables)
}

up_genes <- load_data(file.path('data/Cluster specific genes/Up'))

#up_genes <- up_genes[!grepl("(^Gm|^RP[0-9]|Rik|Hist)", up_genes$GENEID),]

gene_names <- up_genes[up_genes$padj<0.05 & up_genes$log2FoldChange >1,] %>%
  group_by(Cluster) %>%
  dplyr::arrange(Cluster, log2FoldChange) %>%
  top_n(n = 30, wt=log2FoldChange) 
genes <- unique(as.character(gene_names$GENEID))

exp <- as.data.frame(t(as.matrix(sc@ndata[genes,])+0.1))
nam <- as.data.frame(sc@cpart[sc@cpart %in% retain_cl])
nam$newid <- paste(rownames(nam),nam$`sc@cpart[sc@cpart %in% retain_cl]`,sep = "_cl")
exp.new <- merge(exp, nam, by = "row.names")
rownames(exp.new) <- exp.new$newid 
exp.new <- exp.new[,c(3:ncol(exp.new)-1)] # here you have to increase 11 to one more number per gene you will add
exp.new$`sc@cpart[sc@cpart %in% retain_cl]` <- factor(exp.new$`sc@cpart[sc@cpart %in% retain_cl]`, levels = ord_clust)
exp.new <- exp.new[order(exp.new$`sc@cpart[sc@cpart %in% retain_cl]`),]
log10exp <- log10(exp.new[,-ncol(exp.new)])
colnames(log10exp) <- gsub('_.*', '', colnames(log10exp))

#find position of column breaks for heatmap (basically where the new cluster starts)
breaks <- c(table(exp.new$`sc@cpart[sc@cpart %in% retain_cl]`))
breaks <- as.numeric(breaks)
a <- c()
for (i in 1:length(breaks)) {
  if (i==1) {
    a <- breaks[1]
  } else {a <- c(a, a[i-1]+breaks[i])}
}

#define annotation colors from url: https://stackoverflow.com/questions/33292067/pheatmap-annotation-colors-and-border
annotation_col = data.frame(ID = factor(exp.new$`sc@cpart[sc@cpart %in% retain_cl]`))
rownames(annotation_col)<-rownames(exp.new)
cols <- c(colors_pat, colors_many)[1:length(retain_cl)] #colorRampPalette(c(brewer.pal(length(breaks) ,name = 'Set3'), brewer.pal(length(breaks) - 10,name = 'Set2')))
names(cols) <- unique(annotation_col$ID)
annotation_colors <- list(ID = cols)


pheat <- pheatmap(t(log10exp),
                  show_colnames = F, 
                  #color = inferno(1000), 
                  #color = colorRampPalette(c("#FC8D59", "#FFFFBF", "#91CF60"))(length(palette.breaks) - 1),
                  #color = viridis(1000),
                  color = colorRampPalette(c('blue', 'white', 'red'))(1000),
                  cluster_cols = F, 
                  annotation_legend = T, 
                  cluster_rows = T, 
                  fontsize_row = 10, 
                  scale = 'column',
                  gaps_col=a[-length(a)],
                  annotation_col = annotation_col, 
                  annotation_colors = annotation_colors[1])


pdf('plots/heatmaps/single-cell-top-30-heatmap.pdf',height = 12, width = 12)
pheat
dev.off()


mat <- as.matrix(sc@ndata)
nd <- mat[genes,]
nd[is.na(nd)] <- 0
nd  <- t(nd[complete.cases(nd),])
clust_n <- sc@cpart[sc@cpart %in% retain_cl]
#clust_n <- clust_n[clust_n != 7 & clust_n != 8 & clust_n != 9]
mnd <- as.data.frame(cbind(clust_n,nd[rownames(nd) %in% names(clust_n),]))
mean_mnd <- aggregate(mnd[, 2:dim(mnd)[2]], list(mnd$clust_n), mean)
row.names(mean_mnd) <- paste("C",mean_mnd$Group.1,sep = "")
mean_mnd$Group.1 <- NULL
mean_mnd <- as.matrix(mean_mnd) + 0.1
gene <- as.data.frame(t(log10(mean_mnd)))
gene <- gene[complete.cases(gene),]
#row.names(gene) <- id2name(rownames(gene))
row.names(gene) <- gsub('__.*', '', row.names(gene))
pheatmap(gene, cluster_cols = T, cluster_rows=T,fontsize_row = 8, border_color = F, show_rownames = T,show_colnames = T, scale = "row")

pdf(paste0('plots/heatmaps/top-30-mean-heatmap.pdf'),height = 12, width = 6)
pheatmap(gene, cluster_cols = T, cluster_rows=T,fontsize_row = 8, border_color = F, show_rownames = T,show_colnames = T, scale = "row")
dev.off()

svg(paste0('plots/heatmaps/top-30-mean-heatmap.svg'),height = 12, width = 6)
pheatmap(gene, cluster_cols = T, cluster_rows=T,fontsize_row = 8, border_color = F, show_rownames = T,show_colnames = T, scale = "row")
dev.off()

