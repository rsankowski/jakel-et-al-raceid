#Differential genes
library(tidyverse)
library(RaceID)

date = Sys.Date()
load('data/sc.Robj')

subDir <- "data/Cluster specific genes"
dir.create('data/Cluster specific genes')
dir.create('data/Cluster specific genes/Up')
dir.create('data/Cluster specific genes/Down')

for (i in unique(sc@cpart)) {
  
  tryCatch({cl <- names(sc@cpart[sc@cpart %in% c(i)])
  rest <- names(sc@cpart[sc@cpart %in% sc@cpart[sc@cpart != i]])
  diffexp <- diffexpnb(getfdata(sc,n=c(rest,cl)), A=rest, B=cl )
  diffexpgenes <- diffexp[["res"]]
  diffexpgenes <- subset(diffexpgenes, diffexpgenes$pval < 0.05)
  diffexpgenes <- subset(diffexpgenes, abs(diffexpgenes$log2FoldChange) > 0)
  diffexpgenes$GENEID <- rownames(diffexpgenes)
  diffexpgenes <- diffexpgenes %>% dplyr::arrange(padj)
  #diffexpgenes$GENEID <- gsub('__.*', '', diffexpgenes$GENEID)
  diffexpgenes$Cluster <- rep(i, nrow(diffexpgenes))
  diffgene_up <- diffexpgenes[diffexpgenes$log2FoldChange > 0, ]
  diffgene_down <- diffexpgenes[diffexpgenes$log2FoldChange < 0, ]
  write.csv(diffgene_up, file = paste0(subDir, '/Up/', as.character(date), "-diffgenes_cl", as.character(i), "_up_rest.csv"))
  write.csv(diffgene_down, file = paste0(subDir, '/Down/', as.character(date), "-diffgenes_cl", as.character(i), "_down_rest.csv"))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  tryCatch({
    png(paste0(subDir, '/Up/', as.character(date),'-MA-plot-Cl', as.character(i) ,'_rest.png'), res = 300, width =7, height = 7, units = 'in')
    plotdiffgenesnb(diffexp,show_names = T, pthr=.01, lthr = .0001, mthr = .0001, padj=T)
    dev.off()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  for (j in unique(sc@cpart)) {
    tryCatch({
      cl <- names(sc@cpart[sc@cpart %in% i])
      rest <- names(sc@cpart[sc@cpart %in% j])
      diffexp <- diffexpnb(getfdata(sc,n=c(rest,cl)), A=rest, B=cl )
      diffexpgenes <- diffexp[["res"]]
      diffexpgenes <- subset(diffexpgenes, diffexpgenes$pval < 0.05)
      diffexpgenes <- subset(diffexpgenes, abs(diffexpgenes$log2FoldChange) > 0)
      diffexpgenes$GENEID <- rownames(diffexpgenes)
      diffexpgenes <- diffexpgenes %>% dplyr::arrange(padj)
      #diffexpgenes$GENEID <- gsub('__.*', '', diffexpgenes$GENEID)
      diffexpgenes$Cluster <- rep(i, nrow(diffexpgenes))
      diffgene_up <- diffexpgenes[diffexpgenes$log2FoldChange > 0, ]
      diffgene_down <- diffexpgenes[diffexpgenes$log2FoldChange < 0, ]
      write.csv(diffgene_up, file = paste0(subDir, '/Up/', as.character(date), "-diffgenes_cl", as.character(i), "_up_vs_cl", as.character(j), ".csv"))
      write.csv(diffgene_down, file = paste0(subDir, '/Down/', as.character(date), "-diffgenes_cl", as.character(i), "_down_vs_cl", as.character(j), ".csv"))
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    tryCatch({
      png(paste0(subDir, '/Up/', as.character(date), '-MA-plot-Cl', as.character(i) , "_vs_cl", as.character(j), ".png"), res = 300, width =7, height = 7, units = 'in')
      plotdiffgenesnb(diffexp,show_names = T, pthr=.01, lthr = .0001, mthr = .0001, padj=T)
      dev.off()
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    
  }
  
}

#export diffgenes

load_data <- function(path) { 
  files <- dir(path, pattern = '\\.csv', full.names = TRUE)
  tables <- lapply(files, read.csv)
  do.call(rbind, tables)
}

up_genes <- load_data(file.path('data/Cluster specific genes/Up'))
up_genes <- up_genes[!grepl("(^Gm|^Rp|Rik)", up_genes$GENEID),]
up_genes %>% 
  filter(padj<0.05, log2FoldChange > 1) %>% 
  dplyr::arrange(Cluster) %>%
  dplyr::distinct(Cluster, GENEID) %>%
  write.csv('data/up-genes.csv')


up_genes %>% 
  filter(padj<0.05, log2FoldChange > 1) %>% #, up_genes$log2FoldChange > .5
  dplyr::arrange(Cluster) %>%
  dplyr::distinct(Cluster, GENEID) %>%
  write.csv('data/unique-up-genes.csv')


