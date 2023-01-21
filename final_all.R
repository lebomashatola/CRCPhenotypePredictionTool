
#SET-UP


read_files <- function(){
  
  setwd( '~/Documents/rnaseq/datasets' )
  files <- list.dirs( path = ".", full.names = FALSE, recursive = TRUE )
  files <- files[-1]
  
  return(files)
  
}

files <- read_files()



#BUILD REFERENCE GENOME


build_index <- function(){
  
  setwd( '~/Documents/rnaseq' )
  system( 'hisat2-build -p 8 genome.fa genome' )
  
}

build_index()


#QUALITY CONTROL USING TRIMMOMATIC

qc <- function(files){
  
  
  for( i in files ){
    
    if (i == ''){
      
      invisible()
      
    } else {
      
      list <- list.files( path = paste( '~/Documents/rnaseq/datasets/', i, sep='' ), 
                          pattern = '*.fastq.gz$', 
                          all.files = FALSE,full.names = TRUE, recursive = FALSE, 
                          ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE )
      
      setwd( '~/documents/rnaseq/tools/trimmomatic-0.36' )
      system(paste('java -jar trimmomatic-0.36.jar PE -threads 8 -phred33 ', list[1], ' ', list[2], 
                   
                   ' ~/Documents/rnaseq/Trimmed/', i, '/', i, '_forward_paired.fastq.gz ', 
                   '~/Documents/rnaseq/Trimmed/', i, '/', i,'_reverse_paired.fastq.gz ', 
                   
                   'ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 ',
                   'SLIDINGWINDOW:4:15 MINLEN:36', sep=''))
      
    }
  }
}


qc(files)


#ALIGNMENT HISAT2 



align <- function(files){
  

  for( i in files ){
    
    list <- list.files( path = paste( '~/Documents/rnaseq/', i, sep='' ),
                        pattern = '*fastq.gz$', all.files = FALSE, full.names = TRUE, 
                        recursive = FALSE, ignore.case = TRUE, 
                        include.dirs = FALSE, no.. = FALSE )
    
    
    setwd( '~/Documents/rnaseq' )
    system( paste( 'hisat2 -p 8 -x genome -1 ', list[1], ' -2 ', list[2], 
                   ' -S ~/Documents/rnaseq/BAM/', i, '.bam', sep='' ))
    
    
  }
} 

align(files)




#SORT BAM FILES

sorted_bam <- function(files){
  
  
  for( i in files ){
    
    setwd( '~/Documents/rnaseq' )
    system( paste( 'samtools sort ~/Documents/rnaseq/BAM/', i, '.bam -o ~/Documents/rnaseq/BAM_SORTED/', 
                   i, '.sorted.bam -@ 8', sep=''))
      
    
    }
  }
}


sorted_bam(files)


#ASSEMBLY AND QUANTIFICATION USING STRINGTIE

assembly <- function(files){
  
  
  for ( i in files ){
    
    setwd( '~/Documents/rnaseq' )
    system( paste( 'stringtie -p 8 -G ~/Documents/rnaseq/genome.gtf -o ~/Documents/rnaseq/GTF/', i, '.gtf ',
                    '~/Documents/rnaseq/BAM_SORTED/', i, '.sorted.bam', sep='' ))
    
  }
}

assembly(files)


#Merge into non-redundant file 


merge <- function( files ){
  
  list <- list.files( path = '~/Documents/rnaseq/gtf',
                      pattern = '*.gtf', all.files = FALSE, full.names = TRUE, 
                      recursive = FALSE, ignore.case = TRUE, 
                      include.dirs = FALSE, no.. = FALSE )
  
  
  setwd( '~/documents/rnaseq/' )
  fileConn <- file( 'gtf_list.txt' )
  writeLines( list, fileConn )
  close( fileConn )
  
  
  setwd ('~/Documents/rnaseq')
  system( paste( 'stringtie --merge -G ~/Documents/rnaseq/genome.gtf -p 8 -o ',
            'stringtie_merged.gtf gtf_list.txt', sep='' ))
          
}
  

merge(files)




#Finish Assembly and Quantification 


finish <- function ( files ){
  
  
  list <- list.files( path = '~/Documents/rnaseq/gtf',
                      pattern = '*gtf$', all.files = FALSE, full.names = TRUE, 
                      recursive = FALSE, ignore.case = TRUE, 
                      include.dirs = FALSE, no.. = FALSE )
  
  setwd('~/Documents/rnaseq/GTF')
  for ( i in files ){
    
    system( paste( 'mkdir', i, sep = ' ' ))
    
  }
  
  setwd('~/documents/rnaseq')
  
  for ( i in files ){
    
    system( paste( 'stringtie -e -B -p 8 -G stringtie_merged.gtf -o ',
                   '~/Documents/rnaseq/GTF/', i, '/' , i, '.gtf ~/documents/rnaseq/bam_sorted/', i, 
                   '.sorted.bam -A ~/Documents/rnaseq/ABUN/', i, '.tab ',
                   '-G ~/Documents/rnaseq/genome.gtf -p 8', sep=''))

  }
}


finish(files)



#DEG Analysis

ball_gown <- function(files){
  
  
  library(ballgown)
  library(genefilter)
  
  setwd ( '~/documents/rnaseq' )
  bg_chrX <- ballgown(dataDir = "gtf",
                      samplePattern = "ERR")
  
  bg_chrX_filt <- subset(bg_chrX, "rowVars(gexpr(bg_chrX)) >1", genomesubset=TRUE)
  

  
}


ball_gown(files)


#NORMALISATION AND MERGING DATASETS INTO A SINGLE ONE


preprocess <- function(){
  
  library(data.table)
  library(scales)
  
  setDTthreads(8)
  
  setwd( '~/Documents/rnaseq/datasets' )
  files <- list.dirs( path = ".", full.names = FALSE, recursive = TRUE )
  
  cbind.fill<-function(...){
    
    nm <- list(...) 
    nm<-lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
      rbind(x, matrix(, n-nrow(x), ncol(x))))) 
    
  }
  
  df_1 <- data.frame()
  Normalisation <- toupper(as.character(readline(prompt='Enter normalisation method FPKM/TPM: ')))
  
  
  for( i in 2:length(files) ){
    
    
    gtf <- read.table( paste( '~/Documents/rnaseq/ABUN/', files[i], '.tab', sep='' ), 
                      header = F, sep = "\t", fill = TRUE )
    
    
    colnames(gtf) <- gtf[1,]
    gtf <- gtf[-1, ]
    
    gtf <- gtf[,c("Gene Name", Normalisation)]
    gtf <- gtf[!duplicated(gtf[ , c("Gene Name")]),]
    
    colnames(gtf)[2] <- files[i]
    gtf <- gtf[gtf[, 2] != '0.000000', ]
    
    rownames(gtf) <- gtf[,1]
    df_1 <- cbind.fill(df_1, gtf)
    
    
  }
  
  df_1 <- df_1[ , !c(TRUE,FALSE) ]
  df_1 <- as.data.frame(df_1)
  
  
  for (i in 1:ncol(df_1)){
    
    df_1[i] <- as.numeric(unlist(df_1[i]))
    
  }
  
  df_all <- sapply(df_1, rescale, to = c(1, 5))
  rownames(df_all) <- rownames(df_1)
  return(df_all)
  
}


mtrx <- preprocess()



#DIFFERENTIAL GENE EXPRESSION 



sigde <- function(input, fc){
  
  library(limma)
  library(data.table)
  library(dplyr)
  
  design <- model.matrix(~0 + factor(c(rep(1, 1), rep(2, 1), rep(3, 1))))
  colnames(design) <- c("Case", "Control", "Normal")
  
  fit <- lmFit(input, design, method="ls")
  contrast.matrix <- makeContrasts(Case-Control, 
                                   Case-Normal, 
                                   Normal-Control,
                                   levels=design)
  
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE)
  
  results <- decideTests(fit2, method="separate", 
                         adjust.method="BH", 
                         p.value=0.05, lfc=2)
  
  setwd('~/documents/rnaseq/results')
  
  tiff("MA1.tiff", 
       units="in", 
       width=6, 
       height=6, 
       res=150)
  
  par(mfrow=c(1, 1))
  limma::plotMA(fit2, coef=1,status=results[, 1],values=c(1,-1),
                hl.col=c("red","blue"))
  
  dev.off()
  
  tiff("MA2.tiff", 
       units="in", 
       width=6, 
       height=6, 
       res=150)
  
  par(mfrow=c(1, 1))
  limma::plotMA(fit2, coef=2,status=results[, 2],values=c(1,-1),
                hl.col=c("red","blue"))

  
  dfList_up <- c()
  dfList_dwn <- c()
  
  for (k in 1:3){
    
    sig <- topTable(fit2,
                    n=Inf,
                    adjust="BH",
                    coef=k,
                    sort.by="P",
                    p.value=0.05,
                    lfc=fc)
    
    up <- sig[which(sig$logFC > 1), ]
    dwn <- sig[which(sig$logFC < 1), ]
    
    dfList_up <- c(dfList_up, k)
    dfList_dwn <- c(dfList_dwn, k)
    
  }
  
  sig_genes <- list(dfList_up, dfList_dwn)
  genes <- c()
  
  for (r in 1:2){
    
    int_genes <- Reduce(intersect, lapply(sig_genes[[r]], rownames))
    df_data <- input[int_genes, ]
    df_data <- na.omit(df_data)
    genes <- c(genes, df_data)
    
  }
  
  return(genes)
  
}


de_df <- sigde(mtrx)


#RUNNING WGNCA SOFT-THRESHOLD DETERMINATION

pwr <- function(dataset){
  
  library(WGCNA)
  
  enableWGCNAThreads(nThreads = 8)
  allowWGCNAThreads(nThreads = 8)
  options(stringsAsFactors = FALSE)
  
  
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  sft = pickSoftThreshold(dataset, powerVector = powers, verbose = 5)
  
  setwd('~/Documents/Results')
  
  tiff("power1.tiff",
       units="in",
       width=6,
       height=6,
       res=150)
  
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
  
  dev.off()
  
}

pwr(t(mtrx))


#DETERMINING TOM BASED DISSIMALRITY MATRIX AND PLOTTING GENE DENDROGRAM

TOM <- function(dataset, pwr){
  
  require(parallel)
  require(doParallel)
  
  cpucores <- makeCluster(detectCores(), type='PSOCK') 
  registerDoParallel(8) 
  
  library(WGCNA)
  
  dataset <- as.matrix(dataset)
  dataset <- t(dataset)
  
  gsg = goodSamplesGenes(dataset, verbose = 3)
}
  
  
  adj <- adjacency(dataset, power = pwr);
  dissTOM <- TOMdist(adj)
  hierTOM <- hclust(as.dist(dissTOM), method="average")
  HierClust <- labels2colors(cutreeDynamic(hierTOM, method="tree", cutHeight=0.98))
  
  df_in <- fixDataStructure(dataset, verbose = 0, indent = 0)
  
  bnet <- blockwiseConsensusModules(df_in, maxBlockSize = 2000, power = pwr, minModuleSize = 30,
                                    deepSplit = 2,pamRespectsDendro = FALSE,mergeCutHeight = 0.25, 
                                    numericLabels = TRUE,minKMEtoStay = 0,saveTOMs = TRUE, verbose = 5)
  
  moduleLabels <- bnet$colors
  moduleColors <- labels2colors(moduleLabels)
  
  bwLabels <- matchLabels(bnet$colors, moduleLabels, pThreshold = 1e-5)
  bwColors <- labels2colors(bwLabels)
  
  setwd('~/Documents/Results/')
  nBlocks <- length(bnet$dendrograms)
  
  for (block in 1:nBlocks){
    
    tiff(paste("TOMdendro_", block, ".tiff", sep=''), 
         units="in", 
         width=6, 
         height=6, 
         res=150)
    
    plotDendroAndColors(bnet$dendrograms[[block]], 
                        bwColors[bnet$blockGenes[[block]]],
                        "Module colors",
                        main = paste("Gene dendrogram and module colors in block", block),
                        dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, 
                        guideHang = 0.05,setLayout = FALSE)
    
    dev.off()
    
  }
  
  consTree <- bnet$dendrograms
  
  tiff("TOMdendroALL.tiff", 
       units="in", 
       width=6, 
       height=6, 
       res=150)
  
  plotDendroAndColors(consTree,
                      cbind(moduleColors, bwColors),
                      c("Single block", "Blockwise"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Single block consensus gene dendrogram and module colors")
  dev.off()
  
  setwd('~/Documents/Results/')
  
  tiff("TOMdendro.tiff", 
       units="in", 
       width=6, 
       height=6, 
       res=150)
  
  
  par(mfrow = c(1,1))
  
  plotDendroAndColors(hierTOM,
                      colors = data.frame(HierClust),
                      dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                      main = "Gene dendrogram", cex.axis = 1.1)
  
  dev.off()
  
  datME <- moduleEigengenes(dataset, HierClust)$eigengene
  MET <- orderMEs(cbind(datME))
  
  tiff("eigegene.tiff", 
       units="in", 
       width=6, 
       height=6, 
       res=150)
  
  par(mfrow=c(1, 1))
  plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), 
                        marHeatmap = c(3,4,1,2), cex.lab = 0.8, 
                        xLabelsAngle = 90)
  
  dev.off()
  
  datKME <- signedKME(dataset, datME, outputColumnName="")
  
  return(list(datKME, HierClust))
  
}


mcf_tom <- TOM(head(mtrx, 5000), as.integer(readline(prompt="Enter soft-threasholding power: ")))


#Enrichment using CLusterProfiler 

enrichment <- function(dataset, colorh1, datKME, anno, name){
  
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  
  
  intModules <- table(colorh1)
  intModules <- as.data.frame(intModules)
  intModules <-intModules$colorh1
  intModules <- as.character(intModules)
  
  dat <- data.frame() 
  colrs <- c()
  
  cbind.fill<-function(...){
    
    nm <- list(...) 
    nm<-lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
      rbind(x, matrix(, n-nrow(x), ncol(x))))) 
    
  }
  
  for (color in intModules){
    
    FilterGenes <- abs(subset(datKME, select=c(color))) > 0.75
    
    genes <- dimnames(data.frame(dataset))[[2]][FilterGenes]
    entr <- dimnames(data.frame(dataset))[[2]][FilterGenes]
    
    entr <- filter(anno, Symbol %in% entr)
    dat <- cbind.fill(dat, entr$Entrez, fill = NA)
    colrs <- append(color, colrs)
    
  }
  
  dat <- as.data.frame(dat)
  dat <- dat[,seq(1,ncol(dat),2)]
  
  
  colnames(dat) <- colrs
  dat <- dat[, -1]
  dat <- as.data.frame(dat)
  dat <- subset(dat, select = -c(grey))
  
  ck <- compareCluster(geneCluster = dat, fun = "enrichKEGG",
                       pvalueCutoff = 0.05)
  
  dotplot(ck, color = "p.adjust", showCategory = 10, split = NULL,
          font.size = 14, title = "", by = "geneRatio", includeAll = TRUE)
  
  dotplot(ck, color = "p.adjust", showCategory = 10, split = NULL,
          font.size = 14, title = "", by = "geneRatio", includeAll = TRUE)
  
}


ab <- enrichment(t(test), mcf_tom[[2]], mcf_tom[[1]], gpl570)



if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
install.packages(
  pkgs = "bcbioRNASeq",
  repos = c(
    "https://r.acidgenomics.com",
    BiocManager::repositories()
  )
)

