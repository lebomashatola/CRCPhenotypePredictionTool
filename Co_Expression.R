


cbind.fill <- function (...){

  nm <- list(...)
  nm<-lapply(nm, as.matrix)

  n <- max(sapply(nm, nrow))
  do.call(cbind, lapply(nm, function (x)
    rbind(x, matrix(, n-nrow(x), ncol(x)))))

}


unzip_files <- function(paths, type, tech){

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


prep <- function(vars){

  library(dplyr)

  ana1 <- vars[[1]]

  all <- cbind(ana1[[1]])
  all <- all %>% filter(row_number() <= n()-5)

  coldata_resistant <- data.frame(Sample=ana1[[2]], Group=rep(c('R'), each=vars[[2]]))
  coldata <- rbind(coldata_resistant)

  rownames(coldata) <- coldata$Sample
  return(list(all, coldata))

}

Resis <- unzip_files('/Users/lebomash/documents/phd_2021/scripts2020/ml/biomarkers/data/resistant/', 'R')
resistant <- prep(list(Resis, 766))


library(ggfortify)
pca_res <- prcomp(t(resistant[[1]]))
autoplot(pca_res, colour='blue', label = TRUE, label.size = 5)


preprocess <- function(dataset){

  library(ggfortify)
  pca_res <- prcomp(dataset[[1]])
  pca_res <- unclass(pca_res)
  pca_res <- as.data.frame(pca_res$rotation)

  pca_res <- select(pca_res, PC1, PC2)
  resis <- subset(dataset[[1]], select=rownames(pca_res))

  pca_res <- prcomp(resis)
  print(autoplot(pca_res, colour='blue', label=TRUE))

  return(resis)

}


resistant_2 <- preprocess(resistant)



de <- function(stage){

  library(DESeq2)
  library(devtools)
  library(BiocParallel)
  library(dplyr)
  library(data.table)
  library(ggplot2)


  register(MulticoreParam())

  dds <- DESeqDataSetFromMatrix(countData = as.matrix(stage[[1]]),
                                colData = stage[[2]],
                                design = ~ Group)

  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]

  a <- DESeq2::vst(dds, blind = TRUE, fitType = "local")
  print(plotPCA(a, intgroup=c("Group"), returnData=FALSE))

  pcaData <- DESeq2::plotPCA(a, intgroup=c("Group"), returnData=TRUE)

  pcaData_resistant <- as.data.frame(pcaData[pcaData$Group %like% "R", ])
  #pcaData_resistant <- as.data.frame(subset(pcaData_resistant, PC1 > -50 & PC1 < 10))
  pcaData <- rbind(pcaData_resistant)

  sampls <- counts(dds)
  newData <- sampls[ , rownames(pcaData)]

  coldata <- stage[[2]][colnames(newData),]

  return(list(newData, coldata))

}


resistant_processed <- de(resistant)


de_ <- function(stage){


  register(MulticoreParam())

  dds <- DESeqDataSetFromMatrix(countData = as.matrix(stage[[1]]),
                                colData = stage[[2]],
                                design = ~ Group)

  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]

  a <- DESeq2::vst(dds, blind = TRUE, fitType = "local")
  print(plotPCA(a, intgroup=c("Group"), returnData=FALSE))

  pcaData <- DESeq2::plotPCA(a, intgroup=c("Group"), returnData=TRUE)

  pcaData_resistant <- as.data.frame(pcaData[pcaData$Group %like% "R", ])
  pcaData_sensitive <- as.data.frame(pcaData[pcaData$Group %like% 'A', ])
  pcaData_healthy <- as.data.frame(pcaData[pcaData$Group %like% 'B', ])

  pcaData_healthy <- as.data.frame(subset(pcaData_healthy, PC2 < -20))
  pcaData_sensitive <- as.data.frame(subset(pcaData_sensitive, PC1 > 5 & PC2 < 30))
  pcaData_resistant <- as.data.frame(subset(pcaData_resistant, PC2 > 26))

  pcaData <- rbind(pcaData_resistant, pcaData_sensitive, pcaData_healthy)

  sampls <- counts(dds)
  newData <- sampls[ , rownames(pcaData)]

  coldata <- stage[[2]][colnames(newData),]

  return(list(newData, coldata))

}



resistant_processed_2 <- de_(resistant_processed)



df <- read.csv(file='/Users/lebo/Documents/PhD_2022/ML/Biomarkers/Res.csv', sep='\t')
rownames(df) <- df[,1]
df <- df[,-1]


pwr <- function(dataset){

  library(WGCNA)

  enableWGCNAThreads()
  allowWGCNAThreads()
  options(stringsAsFactors = FALSE)

  dataset <- t(dataset)

  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  sft = pickSoftThreshold(dataset, powerVector = powers, verbose = 5)

  par(mfrow = c(1,2))
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)", cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2,
       ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))

  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
       labels=powers,cex=1, col="red")
  abline(h=0.85,col="red")

  plot(sft$fitIndices[,1], sft$fitIndices[,5], cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2,
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))

  text(sft$fitIndices[, 1], sft$fitIndices[, 5],
       labels=powers, cex=1, col="red")

}

pwr(df)



TOM <- function(dataset, pwr){


  library(WGCNA)

  enableWGCNAThreads()
  allowWGCNAThreads()
  options(stringsAsFactors = FALSE)

  dataset <- t(dataset)

  adj <- abs(cor(dataset, use="p"))^pwr

  dissTOM <- TOMdist(adj)
  hierTOMa <- hclust(as.dist(dissTOM), method="average")

  Gene_Modules <- labels2colors(cutreeDynamic(hierTOMa, method="tree", cutHeight=0.98))
  Gene_Clusters <- labels2colors(cutreeDynamic(hierTOMa, distM= dissTOM , cutHeight = 0.98,
                                               deepSplit=3, pamRespectsDendro = FALSE))


  par(mfrow = c(2,4))
  plotDendroAndColors(hierTOMa,
                      colors = data.frame(Gene_Clusters),
                      dendroLabels = FALSE, marAll = c(2, 19, 9, 13),
                      cex.axis = 1.2)


  return(list(adj, dissTOM, Gene_Modules))

}


cnt3_TOM <- TOM(df, 2)

eigenetic_network <- function(dataset, colorh1, pwr){

  library(purrr)

  ADJ1 <- abs(cor(t(dataset), use="p"))^pwr

  colors <- unique(colorh1)
  Alldegrees1 <- intramodularConnectivity(ADJ1, colorh1)

  datME <- moduleEigengenes(t(dataset),colorh1)$eigengenes
  MET <- orderMEs(cbind(datME))
  datKME <- signedKME(t(dataset), datME, outputColumnName="")

  return(datKME)

}



ein_3 <- eigenetic_network(df, cnt3_TOM[[3]], 2)


enrichment <- function(dataset, colorh1, datKME, threshold){

  library(clusterProfiler)
  library(ReactomePA)
  library(dplyr)
  library(org.Hs.eg.db)
  library(enrichplot)


  intModules <- table(colorh1)
  intModules <- as.data.frame(intModules)
  intModules <-intModules$colorh1
  intModules <- as.character(intModules)

  dat <- data.frame()
  dat_new <- data.frame()
  dat_symbol <- data.frame()

  colrs <- c()

  for (color in intModules){

    color <-  color

    FilterGenes <- abs(subset(datKME, select=c(color))) > threshold
    genes <- dimnames(data.frame(dataset))[[2]][FilterGenes]

    dat <- cbind.fill(dat, genes, fill = NA)
    colrs <- append(color, colrs)

  }

  dat <- dat[,seq(1,ncol(dat),2)]
  colnames(dat) <- colrs
  dat <- as.data.frame(dat)


  dat <- dat[,!names(dat) %in% c("grey")]
  colrs <- colnames(dat)

  for (j in 1:ncol(dat)){

    gene <- dat[, j]
    eg <- bitr(gene, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    eg_symbol <- bitr(gene, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")

    genes <- eg$ENTREZID
    genes_symbol <- eg_symbol$SYMBOL

    dat_new <- cbind.fill(dat_new, genes, fill = NA)
    dat_symbol <- cbind.fill(dat_symbol, genes_symbol, fill=NA)

  }

  dat_new <- dat_new[,seq(1, ncol(dat_new), 2)]
  dat_symbol <- dat_symbol[,seq(1, ncol(dat_symbol), 2)]

  colnames(dat_new) <- colrs
  colnames(dat_symbol) <- colrs

  dat_new <- as.data.frame(dat_new)
  dat_symbol <- as.data.frame(dat_symbol)

  return(list(dat_new, dat_symbol))

}


genes_all <- enrichment(t(de_genes), cnt3_TOM[[3]], ein_3, 0.65)



clusterPr <- function(dat_new){

  dat_new <- dat_new[,!names(dat_new) %in% c("grey")]

  ck <- compareCluster(geneCluster = dat_new, fun = "enrichKEGG",
                       pvalueCutoff = 0.05, pAdjustMethod = "BH")

  cKEGG <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")


  cBp <- compareCluster(geneCluster = dat_new, fun = "enrichGO",
                         pvalueCutoff = 0.05, OrgDb = org.Hs.eg.db,
                         ont = "BP", pAdjustMethod = "BH")

  cGO <- setReadable(cBp, OrgDb = org.Hs.eg.db, keyType="ENTREZID")



  return(list(cKEGG, cGO))

}



genes_all_annotated <- clusterPr(genes_all[[1]])



plots <- function(a){

  dotplot(a[[2]], showCategory=2)
  #cnetplot(a[[2]])

}


plots(genes_all_annotated)
