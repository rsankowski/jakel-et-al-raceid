#RaceID4 
library(tidyverse)
library(viridis)
library(RaceID)

date = Sys.Date()
#download counts file from url: https://stackoverflow.com/questions/28986150/downloading-and-extracting-gz-data-file-using-r
url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE118257&format=file&file=GSE118257%5FMSCtr%5FsnRNA%5FExpressionMatrix%5FR%2Etxt%2Egz'
tmp <- tempfile()
##
download.file(url,tmp)
prdata <- read.delim(gzfile(tmp), sep = '\t', stringsAsFactors = F)

sc <- SCseq(prdata)
# filtering of expression data
a <- apply(prdata, 2, sum)
sc <- filterdata(sc, mintotal=quantile(a, 0.25)) # exlcude the lower quartile of the cells
sc <- CCcorrect(sc, 
                dimR = T, 
                nComp = 20,
                CGenes = c('JUN',
                           'FOS',
                           'ZFP36',
                           'HSPA1A|HSPA1B',
                           'DUSP1',
                           'EGR1',
                           'MALAT1'))

sc <- compdist(sc,metric="pearson")
sc <- clustexp(sc) 

plotsaturation(sc,disp=FALSE)
plotsaturation(sc,disp=TRUE)
plotjaccard(sc)

sc <- clustexp(sc,cln=22,sat=FALSE) 
sc <- findoutliers(sc)
plotbackground(sc)
plotsensitivity(sc)
plotoutlierprobs(sc)
ord_clust <- clustheatmap(sc)
save(ord_clust, file = 'ord_clust.Robj')

pdf(paste0('plots/heatmaps/clustheatmap.pdf'))
clustheatmap(sc, final = T)
dev.off()

sc <- comptsne(sc)
sc <- compfr(sc,knn=10)

plotmap(sc)
plotmap(sc,fr=TRUE)
dev.off()

name2id <- function(x,id) {
  ##  id[sub("\\_\\_chr\\w+","",id) %in% x]
  n <- c()
  for ( j in x ){ n <- append(n,id[grep(paste(j,"(\\_\\_chr\\w+|$|\\|)",sep=""),id)])
  }
  n
}

plotexpmap(sc,name2id("Mrc1", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("Lyve1", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("Cd163", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("Tmem119", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("Cx3cr1", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("Hexb", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("Ptprc", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("Cd3e", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("Itgam", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("Cd8a", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("Cd4", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("H2-Aa", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("Zbtb46", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("Ly6c2", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("Cd177", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("Igkc", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("Wfdc17", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("Cd79a", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("Cst3", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("Nkg7", rownames(sc@ndata)),logsc=F,fr=F)

plotexpmap(sc,name2id("Mrc1", rownames(sc@ndata)),logsc=F,fr=T)
plotexpmap(sc,name2id("Lyve1", rownames(sc@ndata)),logsc=F,fr=T)
plotexpmap(sc,name2id("Cd163", rownames(sc@ndata)),logsc=F,fr=T)
plotexpmap(sc,name2id("Tmem119", rownames(sc@ndata)),logsc=F,fr=T)
plotexpmap(sc,name2id("Cx3cr1", rownames(sc@ndata)),logsc=F,fr=T)
plotexpmap(sc,name2id("Ptprc", rownames(sc@ndata)),logsc=F,fr=T)
plotexpmap(sc,name2id("Cd3e", rownames(sc@ndata)),logsc=F,fr=T)
plotexpmap(sc,name2id("Itgam", rownames(sc@ndata)),logsc=F,fr=T)
plotexpmap(sc,name2id("Cd8a", rownames(sc@ndata)),logsc=F,fr=T)
plotexpmap(sc,name2id("H2-Aa", rownames(sc@ndata)),logsc=F,fr=T)
plotexpmap(sc,name2id("Zbtb46", rownames(sc@ndata)),logsc=F,fr=T)
plotexpmap(sc,name2id("Ly6c2", rownames(sc@ndata)),logsc=F,fr=T)
plotexpmap(sc,name2id("Cd177", rownames(sc@ndata)),logsc=F,fr=T)
plotexpmap(sc,name2id("Igkc", rownames(sc@ndata)),logsc=F,fr=T)
plotexpmap(sc,name2id("Wfdc17", rownames(sc@ndata)),logsc=F,fr=T)
plotexpmap(sc,name2id("Plp1", rownames(sc@ndata)),logsc=F,fr=T)
plotexpmap(sc,name2id("Mog", rownames(sc@ndata)),logsc=F,fr=T)
#plot marker genes

dg <- clustdiffgenes(sc,4,pvalue=.01)
head(dg,25)
types <- sub("(\\_\\d+)$","", colnames(sc@ndata))
genes <- head(rownames(dg)[dg$fc>1],10)
plotmarkergenes(sc,genes,samples=types)

#Save sc file
save(sc, file = 'data/sc.Robj')

