##### params #####
load.data  <- T
do.init <- T
mircol     <- "MIRTAR.TS.0.9"
date <- Sys.Date()


csamp <- c("samples")

############################################ functions start ############################################

if ( load.data ){
  if (do.init){
    data   <- list()
    data.add <- list()
    cdiff <- list()
    sco   <- list()
    pca   <- list()
    entr  <- list()
    ltree <- list()
    gcl   <- list()
    lgres <- list()
    somd   <- list()
    marker <- list()
    net    <- list()
    regentr <- list()
    hyperv  <- list()
    traj    <- list()
  }
  
  gene2iso <- read.csv("/home/roman/Documents/Single cell analysis/wgEncodeGencodeBasicVM9_clean_genes2groups.tsv",sep="\t",header=FALSE)
  ercc     <- read.csv("/home/roman/Documents/Single cell analysis/ERCC_Controls_Analysis_length_ext.txt",sep="\t",header=TRUE)
  for ( n in grep("nb.cell",names(ercc)) ){ercc[,n] <- 5/6 * ercc[,n]} # original CEL-seq
  
  n1 <- 1:97
  n2 <- c(1,98:193)
  n4 <- 1:193
  
  for ( s in csamp ){  
    if (s %in% c("samples")){
      data.add[[s]] <- list()
      for ( m in c("t","b","c") ){
        z <- list()
        #if ( s == "ALL" )  { nl <- c(sDN1M1_1="DN1M1_1",sDN1M1_2="DN1M1_2",sDN1M1_3="DN1M1_3",sDN1M1_4="DN1M1_4",sDN1M2_5="DN1M2_5",sDN1M2_6="DN1M2_6",sDN1M2_7="DN1M2_7",sDN1M2_8="DN1M2_8",sDN2M2_1="DN2M2_1",sDN2M2_2="DN2M2_2",sDN2M2_3="DN2M2_3",sDN2M2_4="DN2M2_4",sDN2_9="DN2_9",sDN2_10="DN2_10",sDN2_11="DN2_11",sDN2_12="DN2_12",sDN3_1="DN3_1",sDN3_2="DN3_2",sDN3_3="DN3_3",sDN3_4="DN3_4",sDN3A_1="DN3A_1",sDN3A_2="DN3A_2",sDN3A_3="DN3A_3",sDN3A_4="DN3A_4",sDN34INT_5="DN34INT_5",sDN34INT_6="DN34INT_6",sDN34INT_7="DN34INT_7",sDN34INT_8="DN34INT_8",sDN4_9="DN4_9",sDN4_10="DN4_10",sDN4_11="DN4_11",sDN4_12="DN4_12",sDP_5="DPold_5",sDP_6="DPold_6",sDP_7="DPold_7",sDP_8="DPold_8",sDP_9="DPnew_9",sDP_10="DPnew_10",sDP_11="DPnew_11",sDP_12="DPnew_12",s4HI8INT_9="4HI8INT_9",s4HI8INT_10="4HI8INT_10",s4HI8INT_11="4HI8INT_11",s4HI8INT_12="4HI8INT_12",sCD4_1="CD4_1",sCD4_2="CD4_2",sCD4_3="CD4_3",sCD4_4="CD4_4",sCD4P2_5="CD4P2_5",sCD4P2_6="CD4P2_6",sCD8_5="CD8_5",sCD8_8="CD8_8",sCD8_9="CD8_9",sCD8_10="CD8_10",sCD8P2_7="CD8P2_7",sCD8P2_8="CD8P2_8",sGD1_1="GD1_1",sGD1_2="GD1_2",sGD1_3="GD1_3",sGD1_4="GD1_4",sGD2_5="GD2_5",sGD2_6="GD2_6",sGD2_7="GD2_7",sGD2_8="GD2_8"); fl <- list(sDN1M1_1=n4,sDN1M1_2=n4,sDN1M1_3=n4,sDN1M1_4=n4,sDN1M2_5=n4,sDN1M2_6=n4,sDN1M2_7=n4,sDN1M2_8=n4,sDN2M2_1=n4,sDN2M2_2=n4,sDN2M2_3=n4,sDN2M2_4=n4,sDN2_9=n4,sDN2_10=n4,sDN2_11=n4,sDN2_12=n4,sDN3_1=n4,sDN3_2=n4,sDN3_3=n4,sDN3_4=n4,sDN3A_1=n4,sDN3A_2=n4,sDN3A_3=n4,sDN3A_4=n4,sDN34INT_5=n4,sDN34INT_6=n4,sDN34INT_7=n4,sDN34INT_8=n4,sDN4_9=n4,sDN4_10=n4,sDN4_11=n4,sDN4_12=n4,sDP_5=n4,sDP_6=n4,sDP_7=n4,sDP_8=n4,sDP_9=n4,sDP_10=n4,sDP_11=n4,sDP_12=n4,s4HI8INT_9=n4,s4HI8INT_10=n4,s4HI8INT_11=n4,s4HI8INT_12=n4,sCD4_1=n4,sCD4_2=n4,sCD4_3=n4,sCD4_4=n4,sCD4P2_5=n4,sCD4P2_6=n4,sCD8_5=n4,sCD8_8=n4,sCD8_9=n4,sCD8_10=n4,sCD8P2_7=n4,sCD8P2_8=n4,sGD1_1=n4,sGD1_2=n4,sGD1_3=n4,sGD1_4=n4,sGD2_5=n4,sGD2_6=n4,sGD2_7=n4,sGD2_8=n4); di <- "" }
        if ( s == "samples" )  { nl <- c(P1_ESN_neg_1="P1_ESN_neg", 
                                         P1_ESN_pos_2="P1_ESN_pos", 
                                         P2_ES_LPS_neg_3="P2_ES_LPS_neg", 
                                         P2_ES_LPS_pos_4="P2_ES_LPS_pos", 
                                         P3_ES_Abx_neg_5="P3_ES_Abx_neg", 
                                         P3_ES_Abx_pos_6="P3_ES_Abx_pos", 
                                         P4_LSN_neg_7 = 'P4_LSN_neg', 
                                         P4_LSN_pos_8 = 'P4_LSN_pos', 
                                         P5_LS_Abx_neg_3 = 'P6_LS_Abx_neg', 
                                         P5_LS_Abx_pos_4 = 'P6_LS_Abx_pos', 
                                         P5_LS_LPS_neg_1 = 'P5_LS_LPS_neg', 
                                         P5_LS_LPS_pos_2 = 'P5_LS_LPS_pos'
        );
        fl <- list()
        for ( i in names(nl) ) fl[[i]] <- n4; di <- "" }
        if ( s == "CO_CB_Micr" )  { nl <- c(Plate1_CB_3="CB_3", 
                                            Plate1_CO_4="CO_4" 
        );
        fl <- list()
        for ( i in names(nl) ) fl[[i]] <- n4;  di <- "" }
        
        
        for ( sl in names(nl) ){
          if ( length(di) > 0 ){
            x <- read.csv(paste("/home/roman/Documents/Single cell analysis/fad-alberto/data/counts",di,"/",sl,".cout",m,".csv",sep=""),sep="\t",header=TRUE)
          }else{
            x <- read.csv(paste("/home/roman/Documents/Single cell analysis/fad-alberto/data/counts",sl,".cout",m,".csv",sep=""),sep="\t",header=TRUE)
          }
          x <- x[,fl[[sl]]]
          x <- merge(data.frame(GENEID=c(as.vector(gene2iso[,1]),as.vector(ercc[,1])),GROUP=c(as.vector(gene2iso[,2]),as.vector(ercc[,1]))),x,by="GENEID",all=TRUE)[,-1]
          names(x)[1] <- "GENEID"
          x[is.na(x[,2]),-1] <- 0
          x <- x[order(x$GENEID),]
          z[[sl]]  <- x
          names(z[[sl]])  <- c("GENEID",paste(nl[[sl]],sub("X","",names(z[[sl]])[-1]),sep="_"))
        }
        for ( i in 1:length(z) ) y <- if ( i == 1 ) z[[i]] else merge(y,z[[i]],by="GENEID")
        row.names(y) <- y$GENEID
        y <- y[,-1]
        y <- y[,apply(y,2,sum)>0]
        if ( m == "t" ){
          data[[s]] <- y
        }else{
          data.add[[s]][[m]] <- y
        }
      }
    }
  }
}

#Create prdata file
prdata <- data[[s]][grep("^(ERCC|mt-)",row.names(data[[s]]),invert=TRUE),]
cs <- apply(prdata,2,sum)
prdata <- prdata[,cs > 500]
cs <- cs[cs > 500]
f <- t(prdata["Kcnq1ot1",])/cs < .02

sigcor <- function(x,y,cthr=.4){
  if ( min(var(x),var(y)) == 0 ) return(NA)
  fit <- lm(x ~ y)
  pv <- as.data.frame(summary(fit)[4])[2,4]
  y <- as.data.frame(summary(fit)[4])[2,1]
  if ( is.na(pv) | is.na(y) ) return( NA )
  z <- sign(y)*sqrt(summary(fit)$r.square)
  if ( is.na(z) ) return(NA)
  if ( pv < .01 & abs(z) >= cthr ) return(z) else return(NA)
}


h <- rep(TRUE,nrow(prdata))
for ( g in c("Kcnq1ot1","Gm10715","Gm42418","Gm10800","Mid1")){
  z <- apply(prdata,1,function(x,y) sigcor(x,y),y=t(prdata[g,]))
  h <- h & ( is.na(z) | z < .65 )
}
prdata <- prdata[h,]

# save prdata file
save(prdata, file =paste0('data/prdata.Robj'))


