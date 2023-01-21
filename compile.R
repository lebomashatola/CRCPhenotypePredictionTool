

#MACOS CATALINA VERSION 10.15.7 2.5 GHz DUAL-CORE INTEL I5 (3RD GEN)
#INTEL HD GRAPHICS 4000 1536 MB
#16 GB 1600 MHz DDR3 MEMORY


###########################################################################################################################################

#BUILD REFERENCE GENOME USING HISAT-BUILD


build_index <- function(){
  
  setwd('~/Documents/rnaseq')
  system( 'hisat2-build -p 4 genome.fa genome')
  
}

build_index()


###########################################################################################################################################

#RUN FASTQC ON DATASETS PRE-MODIFICATION


fastqc <-function(){
  
  setwd ( '~/Documents/rnaseq/datasets' )
  files <- list.dirs ( path = ".", full.names = FALSE, recursive = TRUE )
  
  files <- files [ c(-1)]
  
  for ( i in files ) {
    
    list <- list.files ( path = paste('~/Documents/rnaseq/datasets/', i, sep='' ), pattern = '*.fastq.gz$', all.files = FALSE,
                         full.names = TRUE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE )
    
    system( paste( 'fastqc -t 2', list[[1]], list[2], sep = ' ' ))
    
  }
}


fastqc()


###########################################################################################################################################

#TRIMMOMATIC ON LOW QUALITY DATASETS


qc <- function(){
  
  setwd ( '~/Documents/rnaseq/datasets' )
  files <- list.dirs ( path = ".", full.names = FALSE, recursive = TRUE )
  
  files <- files [ c(-1)]

  for ( i in files ) {
    
    list <- list.files ( path = paste('~/Documents/rnaseq/datasets/', i, sep='' ), 
                         pattern = '*.fastq.gz$', all.files = FALSE,
                         full.names = TRUE, recursive = FALSE, ignore.case = FALSE, 
                         include.dirs = FALSE, no.. = FALSE )
    
    
    option_1 <- readline( paste ('Trim dataset? (Y/N): ', i, ': ', sep=''))
    
    if (option_1 == 'Y'){
      
      setwd ( '~/documents/rnaseq/tools/trimmomatic-0.39' )
      system ( paste ( 'java -jar trimmomatic-0.39.jar PE -threads 2 -phred33 ', list[1], ' ', list[2], 
                       
                       ' ~/Documents/rnaseq/trimmed/', i, '/', i, '_forward_paired.fastq.gz ', 
                       ' ~/Documents/rnaseq/trimmed/', i, '/', i, '_forward_unpaired.fastq.gz ',
                       
                       '~/Documents/rnaseq/trimmed/', i, '/', i,'_reverse_paired.fastq.gz ', 
                       ' ~/Documents/rnaseq/trimmed/', i, '/', i, '_reverse_unpaired.fastq.gz ',
                       
                       'ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36', sep='' ))
      
    }
  }
}
        
qc()


###########################################################################################################################################

#RUN FASTQC ON DATASETS POST-MODIFICATION


fastqc <-function(){
  
  setwd ( '~/Documents/rnaseq/trimmed' )
  files <- list.dirs ( path = ".", full.names = FALSE, recursive = TRUE )
  
  files <- files [ c(-1)]
  
  for ( i in files ) {
    
    list <- list.files ( path = paste('~/Documents/rnaseq/trimmed/', i, sep='' ), pattern = '\\.fastq.gz$', all.files = FALSE,
                         full.names = TRUE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE )
    
    system( paste( 'fastqc -t 2', list[[1]], list[2], sep = ' ' ))
    
  }
}


fastqc()


###########################################################################################################################################

#PREPARATION FOR ALIGNMENT 


prep_hisat <- function(){
  
  
  setwd ( '~/Documents/rnaseq/trimmed' )
  files <- list.dirs ( path = ".", full.names = FALSE, recursive = TRUE )
  
  files <- files [ c(-1)]

  
  for ( k in files ) {
    
    setwd ( paste('~/Documents/rnaseq/trimmed/', k, sep='' ))
    system( paste('rm -rf ~/Documents/rnaseq/trimmed/', k, '/*_unpaired*', sep='' ))
  
  }
}
  
prep_hisat()


###########################################################################################################################################

#ALIGNMENT USING HISAT2 


align <- function(files){

  get_in <- function(opt){
    
    ll <- list.files( path = paste( '~/Documents/rnaseq/', opt, '/', i, sep='' ), pattern = '*fastq.gz$', all.files = FALSE, 
                        full.names = TRUE, recursive = FALSE, ignore.case = TRUE, include.dirs = FALSE, no.. = FALSE )
    return(ll)
    
  }
  
  setwd( '~/Documents/rnaseq/datasets')
  files <- list.dirs( path = ".", full.names = FALSE, recursive = TRUE )
  files <- files [ c(-1)]

  for( i in files ){
    
    setwd( '~/Documents/rnaseq' )
    option <- readline( paste ('Dataset ', i, ' trimmed? (Y/N): ', i, ': ', sep=''))
    
    if (option == 'Y'){
      
      files.in <- get_in('trimmed')
      system( paste( 'hisat2 -p 2 -x genome -1 ', files.in[1], ' -2 ', files.in[2], ' -S ~/Documents/rnaseq/BAM/', i, '.bam', sep='' ))
      
    } else{
      
      files.in <- get_in('datasets')
      system( paste( 'hisat2 -p 2 -x genome -1 ', files.in[1], ' -2 ', files.in[2], ' -S ~/Documents/rnaseq/BAM/', i, '.bam', sep='' ))
      
    }
  }
} 

align(files)


###########################################################################################################################################

#SORT BAM FILES

sorted_bam <- function(){
  
  setwd( '~/Documents/rnaseq/datasets' )
  files <- list.dirs( path = ".", full.names = FALSE, recursive = TRUE )
  files <- files [ c(-1)]
  
  for( i in files ){
    
    setwd('~/Documents/rnaseq')
    system( paste( 'samtools sort ~/Documents/rnaseq/BAM/', i, '.bam -o ~/Documents/rnaseq/BAM_SORTED/', i, '.sorted.bam -@ 2', sep=''))
    
  }
}


sorted_bam()


###########################################################################################################################################

#ASSEMBLY AND QUANTIFICATION USING STRINGTIE 

assembly <- function(){
  
  
  setwd( '~/Documents/rnaseq/datasets' )
  files <- list.dirs( path = ".", full.names = FALSE, recursive = TRUE )
  files <- files [ c(-1)]
  
  for( i in files ){
    
    setwd('~/Documents/rnaseq')
    system( paste( 'stringtie -e -o ~/Documents/rnaseq/GTF/', i, '.gtf ~/Documents/rnaseq/BAM_SORTED/', i, '.sorted.bam -A ',
                   '~/Documents/rnaseq/ABUN/', i, '.tab -G ~/Documents/rnaseq/genome.gtf -p 2', sep=''))
    
  }
}

assembly()


 ###########################################################################################################################################

#COMPILE FPKM/TPM

preprocess <- function(){
  
  library(data.table)
  
  setDTthreads(4)
  
  setwd( '~/Documents/rnaseq/datasets' )
  files <- list.dirs( path = ".", full.names = FALSE, recursive = TRUE )
  files <- files [ c(-1)]
  
  cbind.fill<-function(...){
    
    nm <- list(...) 
    nm<-lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
      rbind(x, matrix(, n-nrow(x), ncol(x))))) 
    
  }
  
  df_1 <- data.frame()
  Normalisation <- toupper(as.character(readline(prompt='Enter normalisation method FPKM/TPM: ')))
  
  for( i in 1:length(files) ){
    
    
    gtf <- read.table(paste('~/Documents/rnaseq/ABUN/', files[i], '.tab', sep=''), 
                      header = F, sep = "\t", fill = TRUE)
    
    
    colnames(gtf) <- gtf[1,]
    gtf <- gtf[-1, ]
    
    gtf <- gtf[,c("Gene Name", Normalisaton)]
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
  
  sampling <- function(){
    
    a <- as.integer(readline(prompt='Enter number of samples: '))
    df_in <- c()
    
    for (i in 1:a){
      
      df_append <- readline(prompt='Enter sample name: ')
      df_in <- append(df_in, df_append)
    }
    
    return(df_in)
  }
  
  
  sample_names <- toupper(readline(prompt = 'Change Sample Names (Y/N): '))
  
  if (sample_names == 'N'){
    
    return(df_1)
    
  }
  
  if (sample_names == 'Y'){
    
    colnames(df_1) <- sampling()

  }
  
  
  replace_missings <- function(x, replacement) {
    
    is_miss <- is.na(x)
    x[is_miss] <- replacement
    
    message(sum(is_miss), " missings replaced by the value ", replacement)

    
  }
  
  return(replace_missings(mtrx, 0))
  
}

mtrx <- preprocess()





###########################################################################################################################################

#DIFFERENTIAL CO-EXPRESSION ANALYSIS USING DIFFCOEXP 

diff_co <- function(data_in){
  
  library(diffcoexp)
  allowWGCNAThreads(4)
  
  exprs.1 <- data_in[, 1:3] 
  exprs.1 <- head(exprs.1, 5000)
  
  exprs.2 <- data_in[, 4:6]
  exprs.2 <- head(exprs.2, 5000)
  
  res <- diffcoexp(exprs.1, exprs.2, r.method = "pearson", 
                   q.method = "BH", 
                   rth = 0.5, qth = 0.1, r.diffth = 0.5, q.diffth = 0.1,
                   q.dcgth = 0.1)
  
  
}

head(diff_co(mtrx))


###########################################################################################################################################


#IDENTIFY CO-EXPRESSED GENE CLUSTERS 


###########################################################################################################################################


#CLUSTERPROFILER KEGG/GO/DO 