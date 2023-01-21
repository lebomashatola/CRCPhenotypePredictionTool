



############################################################################################################################



cbind.fill<-function (...){
  
  nm <- list(...) 
  nm<-lapply(nm, as.matrix)
  
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
  
}


############################################################################################################################



unzip_files <- function(paths, type){
  
  
  setwd( paths )
  control <- list.dirs ( path = ".", full.names = F, recursive = F )
  
  out <- NULL
  df_meta <- data.frame()
  
  count_c <- c()
  count <- 1
  
  for ( file_1 in control ){
    
    setwd( paste(paths, file_1, sep='' ))
    
    files <- list.files(path = ".", recursive = FALSE,
                        pattern = "\\.gz$", 
                        full.names = FALSE)
    
    for ( file in files ){
      
      data <- read.table(gzfile(file[[1]]))
      df_meta <- cbind.fill( df_meta, data[,2])
      
      count_c <- append(count_c, paste(type, count, sep=""))
      count = count + 1
    }
    
    rownames(df_meta) <- data$V1
    all_genes <- as.character(rownames(df_meta))
    all_genes <- sub("\\.\\d+", "", all_genes)
    
    rownames(df_meta) <- all_genes
    
  }
  
  colnames(df_meta) <- count_c
  df_meta <- as.data.frame(df_meta)
  
  return(list(df_meta, count_c))
  
}
    

sta1 <- unzip_files('~/documents/phd/samples/count_data/crc_stage1/', 'S')
sta2 <- unzip_files('~/documents/phd/samples/count_data/crc_stage2/', 'R')
sta3 <- unzip_files('~/documents/phd/samples/count_data/crc_stage3/', 'R')
sta4 <- unzip_files('~/documents/phd/samples/count_data/crc_stage4/', 'R')
rec <- unzip_files('~/documents/phd/samples/count_data/recurrent_crc/', 'R')



prep_data <- function(ana1, num1, ana2, num2){
  
  
  all <- cbind(ana2[[1]], ana1[[1]])
  coldata_stage <- data.frame(Sample=ana2[[2]], Group=rep(c('R'), each=num2))
  coldata_norm <- data.frame(Sample=ana1[[2]], Group=rep(c('S'), each=num1))
  
  coldata <- rbind(coldata_stage, coldata_norm)
  rownames(coldata) <- coldata$Sample
  
  return(list(all, coldata))
  
  
}

stage2 <- prep_data(sta1, 113, sta2, 235)
stage3 <- prep_data(sta1, 113, sta3, 210)
stage4 <- prep_data(sta1, 113, sta4, 95)
rec <- prep_data(sta1, 113, rec, 4)


############################################################################################################################





de <- function(stage){
  
  library(DESeq2)
  library(dplyr)
  library(devtools)
  library(BiocParallel)
  
  register(MulticoreParam())
  
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(stage[[1]]),
                                  colData = stage[[2]],
                                  design = ~ Group)
  
  dds$Group <- relevel(dds$Group, ref = "S")
    
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds <- DESeq(dds, parallel=TRUE)

  res <- results(dds, contrast=c("Group", "R", "S"), alpha=0.05, parallel=TRUE)
  resOrdered <- res[order(res$padj), ]
  
  all_genes <- as.character(rownames(resOrdered))
  all_genes <- sub("\\.\\d+", "", all_genes)
  rownames(resOrdered) <- all_genes
    
  up <- resOrdered[which(resOrdered$log2FoldChange > 1), ]
  dwn <- resOrdered[which(resOrdered$log2FoldChange < 1), ]
  
  up <- up[which(up$padj < 0.05), ]
  dwn <- dwn[which(dwn$padj < 0.05), ]
  
  all <- stage[[1]] 
  up_genes <- all[rownames(up), ]
  up_genes <- varianceStabilizingTransformation(as.matrix(up_genes), blind = TRUE, fitType = "local")
  
  return(up_genes)

}

cnts_2 <- de(stage2)
cnts_3 <- de(stage3)
cnts_4 <- de(stage4)
cnts_rec <- de(rec)




############################################################################################################################



get_files <- function(files_in, label){
  
  
  library(clusterProfiler)
  library(dplyr)
  library(org.Hs.eg.db)
  library(enrichplot)
  
  
  up_genes <- files_in[[1]]
  dwn_genes <- files_in[[2]]
  
  eg_up <- bitr(rownames(up_genes), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  eg_dwn <- bitr(rownames(dwn_genes), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  write.csv(eg_up, paste('~/downloads/up', label, '.csv', sep=''))
  write.csv(eg_dwn, paste('~/downloads/dwn', label, '.csv', sep=''))
  
}



get_files(cnts_2, 'stage2')
get_files(cnts_3, 'stage3')
get_files(cnts_4, 'stage4')
get_files(cnts_rec, 'recurrent')


############################################################################################################################





functional <- function(){
  
  
  library(clusterProfiler)
  library(dplyr)
  library(org.Hs.eg.db)
  library(enrichplot)
  
  
  dat_dwn <- read.csv('~/documents/phd/samples/count_data/dwn_stage.csv', sep=';')
  dat_up <- read.csv('~/documents/phd/samples/count_data/up_stage.csv', sep=';')
  dat <- read.csv('~/documents/phd/samples/count_data/rec.csv', sep=';')
  
  
  ck_up <- compareCluster(geneCluster = dat, 
                          fun = "enrichKEGG", pvalueCutoff = 0.05, pAdjustMethod = "BH")
  
  
  cbp_up <- compareCluster(geneCluster = dat, fun = "enrichGO", 
                           pvalueCutoff = 0.05, OrgDb = org.Hs.eg.db, 
                           ont = "BP", pAdjustMethod = "BH")
  
  
  cmf_up <- compareCluster(geneCluster = dat, fun = "enrichGO", 
                           pvalueCutoff = 0.05, OrgDb = org.Hs.eg.db, 
                           ont = "MF", pAdjustMethod = "BH")
  
  
  ccc_up <- compareCluster(geneCluster = dat, fun = "enrichGO", 
                           pvalueCutoff = 0.05, OrgDb = org.Hs.eg.db, 
                           ont = "CC", pAdjustMethod = "BH")
  
  
  
  
  return(list(ck_up, cbp_up, cmf_up, ccc_up))
  
}

enrich <- functional()


dotplot(enrich[[4]])




############################################################################################################################





pwr <- function(dataset){
  
  library(WGCNA)
  
  enableWGCNAThreads()
  allowWGCNAThreads()
  options(stringsAsFactors = FALSE)
  
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  sft = pickSoftThreshold(dataset, powerVector = powers, verbose = 5)
  
  par(mfrow = c(1,2))
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",
       ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
       labels=powers,cex=0.9, col="red")
  abline(h=0.85,col="red")
  
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  
  text(sft$fitIndices[, 1], sft$fitIndices[, 5],
       labels=powers, cex=0.9, col="red")
  
  
}


pwr(t(cnts_2))         
pwr(t(cnts_3))    
pwr(t(cnts_4))    
pwr(t(cnts_rec))    


############################################################################################################################




TOM <- function(dataset){
  
  library(WGCNA)
  
  enableWGCNAThreads()
  allowWGCNAThreads()
  options(stringsAsFactors = FALSE)
  
  adj <- abs(cor(dataset, use="p"))^2
  dissTOM <- TOMdist(adj)
  
  hierTOMa <- hclust(as.dist(dissTOM), method="complete")
  
  Gene_Modules <- labels2colors(cutreeDynamic(hierTOMa, method="tree", cutHeight=0.999))
  Gene_Clusters <- labels2colors(cutreeDynamic(hierTOMa, distM= dissTOM , cutHeight = 0.9999,
                                               deepSplit=3, pamRespectsDendro = FALSE))
  
  
  par(mfrow = c(1,1))
  plotDendroAndColors(hierTOMa,
                      colors = data.frame(Gene_Clusters),
                      dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                      cex.axis = 1.2)
  
  
  return(list(adj, dissTOM, Gene_Modules))
  
}



cnt_TOM_2 <- TOM(t(cnts_2))
cnt_TOM_3 <- TOM(t(cnts_3))
cnt_TOM_4 <- TOM(t(cnts_4))
cnt_TOM_rec <- TOM(t(cnts_rec))



############################################################################################################################




eigenetic_network <- function(dataset, colorh1){
  
  library(purrr)
  
  datME <- moduleEigengenes(t(dataset),colorh1)$eigengenes
  MET <- orderMEs(cbind(datME))
  datKME <- signedKME(t(dataset), datME, outputColumnName="")
  
  return(datKME)
  
}

ein_2 <- eigenetic_network(cnts_2, cnt_TOM_2[[3]])
ein_3 <- eigenetic_network(cnts_3, cnt_TOM_3[[3]])
ein_4 <- eigenetic_network(cnts_4, cnt_TOM_4[[3]])
ein_rec <- eigenetic_network(cnts_rec, cnt_TOM_rec[[3]])



############################################################################################################################




enrichment <- function(dataset, colorh1, datKME){
  
  library(clusterProfiler)
  library(dplyr)
  library(org.Hs.eg.db)
  library(enrichplot)
  
  
  intModules <- table(colorh1)
  intModules <- as.data.frame(intModules)
  intModules <-intModules$colorh1
  intModules <- as.character(intModules)
  
  
  dat <- data.frame() 
  dat_new <- data.frame()
  colrs <- c()
  
  for (color in intModules){
    
    FilterGenes <- abs(subset(datKME, select=c(color))) > 0.3
    
    genes <- dimnames(data.frame(dataset))[[2]][FilterGenes]
    
    dat <- cbind.fill(dat, genes, fill = NA)
    colrs <- append(color, colrs)
    
  }
  
  dat <- dat[,seq(1,ncol(dat),2)]
  colnames(dat) <- colrs
  
  dat <- as.data.frame(dat)
  
  for (j in 1:ncol(dat)){
    
    gene <- dat[, j]
    eg <- bitr(gene, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    genes <- eg$ENTREZID
    
    dat_new <- cbind.fill(dat_new, genes, fill = NA)
    
  }
  
  dat_new <- dat_new[,seq(1,ncol(dat_new),2)]
  
  colnames(dat_new) <- colrs
  dat_new <- as.data.frame(dat_new)

  dat_new <- dat_new[,!names(dat_new) %in% c("grey")]
  
  ck <- compareCluster(geneCluster = dat_new, fun = "enrichKEGG",pvalueCutoff = 0.05, pAdjustMethod = "BH")
  
  #cbp <- compareCluster(geneCluster = dat_new, fun = "enrichGO", 
                        #pvalueCutoff = 0.05, OrgDb = org.Hs.eg.db, 
                        #ont = "BP", pAdjustMethod = "BH")
  
  #dotplot(cbp, showCategory=5)
  dotplot(ck)
  
  
}


enrichment(t(cnts_rec), cnt_TOM_rec[[3]], ein_rec)

enrichment(t(cnts_2), cnt_TOM_2[[3]], ein_2)
enrichment(t(cnts_3), cnt_TOM_3[[3]], ein_3)
enrichment(t(cnts_4), cnt_TOM_4[[3]], ein_4)


############################################################################################################################













