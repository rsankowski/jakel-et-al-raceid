#help from url: https://stackoverflow.com/questions/6713973/how-do-i-calculate-correlation-between-corresponding-columns-of-two-matrices-and
library(tidyverse)
library(viridis)
library(corrplot)

date <- Sys.Date()

#load(data)
        jakel <- read.csv('data/2019-02-23normalized-gene-counts-per-celltype.csv', stringsAsFactors = F)
        sank <- read.csv('data/2019-02-23-healthy-micr-normalized-gene-counts.csv', stringsAsFactors = F)
        
#adjust gene names and remove duplicates in sank
        sank$X <- gsub('_.*', '', sank$X)
        sank <- sank[!duplicated(sank$X),]
        sank[,2:20] <- as.data.frame(as.matrix(sank[,2:20])-0.1)
        
#adjust gene names in jakel
        jakel$X <- gsub('\\.', '-', jakel$X)
        
#define common genes        
        a <- sank$X %in% jakel$X
        genes <- sank$X[sank$X %in% jakel$X]
#remove batch effect genes
        genes <- genes[!grepl('^(HSPA1A|MTRNR2L8|MTRNR2L12|HSP90AA1|MALAT1|ZFP36L1|ZFP36|FOS|MALAT1|HSPB*|DUSP1|HSPH1|HSPA*|JUN|HSP90B1|RPS16|DNAJB1|H3F3B|HERPUD1|NEAT1|IVNS1ABP|HIST1H2BG|RP*|XIST|KLF*)', genes)]
        
#trim the gene lists        
        sank <- sank[sank$X %in% genes,]
        jakel <- jakel[jakel$X %in% genes,]
#rename rownames
        rownames(sank) <- sank$X
        rownames(jakel) <- jakel$X
        
#adjust order of rows        
        jakel <- jakel[ order(match(jakel$X,sank$X)), ] # from https://stackoverflow.com/questions/27362718/reordering-rows-in-a-dataframe-according-to-the-order-of-rows-in-another-datafra
        
        
        sank <- sank[,-1]
        jakel <- jakel[,-1]
        

#tutorial from url: https://hemberg-lab.github.io/scRNA.seq.course/comparingcombining-scrnaseq-datasets.html

#define find_vargenes function
       find_vargenes <- function(counts, cv2_cutoff = .3, ...) {
         require(statmod)
         require(DESeq2)  
         lib.size <- estimateSizeFactorsForMatrix(counts)
         ed <-t(t(counts)/lib.size)
         means <- rowMeans(ed)
         vars <- apply(ed,1,var)
         cv2 <- vars/means^2
         
         #fit regression line
         minMeanForFit <- unname( quantile(means[which( cv2 > cv2_cutoff)], 0.95))
         useForFit <- means >= minMeanForFit
         fit <- glmgam.fit( cbind( a0=1, altitude= 1/means[useForFit]), cv2[useForFit])
         a0 <- unname(fit$coefficients["a0"])
         a1 <- unname(fit$coefficients["altitude"])
         
         #plot regression
         xg <- exp(seq(min(log(means[means>0])),max(log(means)), length.out = 1000 ))
         vfit <- a1/xg + a0
         
         #rank genesby the significance of deviation from the fit
         afit <- a1/means + a0
         varFitRatio <- vars/(afit*means^2)
         varorder <- order(varFitRatio, decreasing = T)
         smoothScatter(log(means), log(cv2))
         lines(log(xg), log(vfit), col = "black", lwd=3)
         df<- ncol(ed) -1
         #add conf interval
         lines(log(xg), log(vfit * qchisq(0.975, df)/df), lty=2, col="black")
         lines(log(xg), log(vfit * qchisq(0.075, df)/df), lty=2, col="black")
         
         oed <- ed[varorder,]
         return(oed)
       }
       
       vargenes_sank <- find_vargenes(sank)
       vargenes_jakel <- find_vargenes(jakel)
       
       var_genes <- rownames(vargenes_sank)[1:1000]
       
       #correlation matrices
       jakel2 <- jakel[rownames(jakel) %in% var_genes,] %>% t %>% scale %>% t
       sank2 <- sank[rownames(sank) %in% var_genes,] %>% t %>% scale %>% t
       
       jakel2 <- as.matrix(na.omit(jakel2))
       sank2 <- as.matrix(sank2[rownames(jakel2),])
       
       both_datasets <- cor(cbind(sank2, jakel2), method = "pearson")
       dist(both_datasets)
       
       #Hmisc method url: http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
       library("Hmisc")
       res2 <- rcorr(as.matrix(cbind(sank2, jakel2)), type = "pearson")
       res2
       
       res2$r
       res2$P
       
       #flatten correlation matrix
       flattenCorrMatrix <- function(cormat, pmat) {
         ut <- upper.tri(cormat)
         data.frame(
           row = rownames(cormat)[row(cormat)[ut]],
           column = rownames(cormat)[col(cormat)[ut]],
           cor  =(cormat)[ut],
           p = pmat[ut]
         )
       }
       
       flattenCorrMatrix(res2$r, res2$P)
      
       # Insignificant correlation are crossed
       corrplot(res2$r, type="upper", order="hclust", 
                p.mat = res2$P, sig.level = 0.0001, insig = "blank")
       # Insignificant correlations are leaved blank
       corrplot(res2$r, type="upper", order="hclust", 
                p.mat = res2$P, sig.level = 0.0001, insig = "blank")
       corrplot(res2$r, method= "square", type="upper", order="hclust", 
                p.mat = res2$P, sig.level = 0.0001, insig = "blank")
       
       
       
       
       
       #remove duplicate values from correlation matrix
       pear_flat <- flattenCorrMatrix(res2$r, res2$P)
       
       pear_flat <- pear_flat[grepl("Cluster", pear_flat$row) & !grepl("Cluster", pear_flat$column),]
       
       #ggplot
       clustering <- hclust(dist(both_datasets, method = "maximum"), method = "ward.D2")
       label_order <- clustering$labels[clustering$order]
       
       pear_flat$row <- factor(pear_flat$row, levels = label_order)
       pear_flat$column <- factor(pear_flat$column, levels = label_order)
       pear_flat$padj <- p.adjust(pear_flat$p, method = 'BH')
       
       pear_flat2 <- pear_flat %>% 
         filter(padj <0.01, cor > 0.1)
       pear_flat2$padj[pear_flat2$padj==0] <- 1e-13
       #pear_flat2$p[pear_flat2$padj>0.05] <- NA
       
       #exclude clusters abobe X9
       
       
       #pear_flat2 <- na.omit(pear_flat2[pear_flat2$cor > 0,])
       ggplot(pear_flat2, aes(row, column, fill=-log10(padj))) +
         geom_tile() +
         coord_flip() +
         theme_minimal() +
         theme(axis.text.x = element_text(angle = 90, hjust = 1))+
         scale_fill_viridis(option = 5) +
         labs(y= 'jakel-dataset', x='sankowski-dataset')
ggsave(paste0("plots/others/",date, "-cluster-correlation-jakel-sankowski.pdf"))      

#line plot
pear_flat2 <- pear_flat2[1:10,]
pear_flat2$row <- gsub('luster', '', pear_flat2$row)
trigrams2 <- data.frame('trigram' = paste(pear_flat2$row, pear_flat2$column))

                                      #pear3 <- pear_flat2 %>% gather('dataset', 'cluster' ,row:column)
                                      #pear3$dataset <- ifelse(pear3$dataset == 'row', 'Prinz', 'jakel')
                                      #pear3$dataset <- factor(pear3$dataset, levels = c('Prinz', 'jakel'))
                                      #pear3$direction <- rep(1:nrow(pear_flat2))

#code from # from: https://www.hvitfeldt.me/2018/01/visualizing-trigrams-with-the-tidyverse/
                                      n_word <- nrow(pear_flat2)
                                      n_top <- nrow(pear_flat2)
                                      n_gramming <- 2
                                      
                                      str_nth_word <- function(x, n, sep = " ") {
                                        str_split(x, pattern = " ") %>%
                                          map_chr(~ .x[n])
                                      }
                                      
                                      nodes <- map_df(seq_len(n_gramming),
                                                      ~ trigrams2 %>%
                                                        dplyr::mutate(word = str_nth_word(trigram, .x)) %>%
                                                        dplyr::count(word, sort = TRUE) %>%
                                                        dplyr::slice(seq_len(n_word)) %>% 
                                                        dplyr::mutate(y = seq(from = n_word + 1, to = 0, 
                                                                       length.out = n() + 2)[seq_len(n()) + 1],
                                                               x = .x))


                                      sigmoid <- function(x_from, x_to, y_from, y_to, scale = 5, n = 100) {
                                        require(tidyverse)
                                        library(tidytext)
                                        library(purrrlyr)
                                        
                                        x <- seq(-scale, scale, length = n)
                                        y <- exp(x) / (exp(x) + 1)
                                        tibble(x = (x + scale) / (scale * 2) * (x_to - x_from) + x_from,
                                               y = y * (y_to - y_from) + y_from)
                                      }
                                      
                                      egde_lines <- function(trigram, from_word, to_word, scale = 5, n = 38, 
                                                             x_space = 0) {
                                        require(tidyverse)
                                        library(tidytext)
                                        library(purrrlyr)
                                        
                                        from_word <- from_word %>%
                                          select(-n) %>%
                                          set_names(c("from", "y_from", "x_from"))
                                        
                                        to_word <- to_word %>%
                                          select(-n) %>%
                                          set_names(c("to", "y_to", "x_to"))
                                        
                                        links <- crossing(from = from_word$from, 
                                                          to = to_word$to) %>%
                                          mutate(word_pair = paste(from, to),
                                                 number = map_dbl(word_pair, 
                                                                  ~ sum(str_detect(trigram$trigram, .x)))) %>%
                                          left_join(from_word, by = "from") %>%
                                          left_join(to_word, by = "to")
                                        
                                        links %>%
                                          by_row(~ sigmoid(x_from = .x$x_from + 0.2 + x_space,
                                                           x_to = .x$x_to - 0.05, 
                                                           y_from = .x$y_from, y_to = .x$y_to, 
                                                           scale = scale, n = n) %>%
                                                   mutate(word_pair = .x$word_pair,
                                                          number = .x$number,
                                                          from = .x$from)) %>%
                                          pull(.out) %>%
                                          bind_rows()
                                      }
                                      
                                      
                                      #lines
                                      egde_lines(trigram = trigrams2, 
                                                 from_word = filter(nodes, x == 1), 
                                                 to_word = filter(nodes, x == 2)) %>%
                                        filter(number > 0) %>%
                                        ggplot(aes(x, y, group = word_pair, alpha = number, color = from)) +
                                        geom_line()
                                      
                                      
                                      # egdes between first and second column
                                      edges <- egde_lines(trigram = trigrams2, 
                                                          from_word = filter(nodes, x == 1), 
                                                          to_word = filter(nodes, x == 2), 
                                                          n = 50) %>%
                                        filter(number > 0) %>%
                                        mutate(id = word_pair)
                  
                                      edges$significance <- cut(-log10(edges$padj), 10)
                                      
                                      colnames(trigrams2) <- "word_pair"
                                      edges <- pear_flat2 %>% 
                                        bind_cols(trigrams2) %>% 
                                        full_join(edges)
                                      
                                      #visualization
                                      p <- nodes %>% 
                                        ggplot(aes(x, y, label = word, size = n)) +
                                        geom_text(hjust = 0) + #, color = "#DDDDDD"
                                        theme_void() +
                                        geom_line(data = edges,
                                                  aes(x, y, group = id, color = -log10(padj)), #,, lwd = sqrt(as.numeric(significance)), alpha = sqrt(number)
                                                  lwd = 2,
                                                  inherit.aes = FALSE) +
                                        theme(text = element_text(size = 15)) + #color = "#EEEEEE", plot.background = element_rect(fill = "#666666", colour = 'black'),
                                        guides(alpha = "none", size = "none") +
                                        xlim(c(0.9, 2.2)) +
                                        scale_color_viridis(option = 5) +
                                        labs(title = "Correlation between Prinz and jakel Datasets") + 
                                        scale_size(range = c(3, 8))
                                      p
                                      
                                      ggsave(paste0(date, '-correlation-between-prinz-and-jakel-datasets.pdf'))
                                      sank <- as.data(sank)                                                                            
                                      svg(paste0(date, '-correlation-between-prinz-and-jakel-datasets.svg'), width = 6.23, height = 3.54)
                                      p
                                      dev.off()





























