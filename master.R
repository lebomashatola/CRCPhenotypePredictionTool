



###########################################################################################################################################



#RUN FASTQC ON DATASETS PRE-MODIFICATION


fastqc <-function(){
  
  setwd ( '~/Documents/rnaseq/datasets' )
  files <- list.dirs ( path = ".", full.names = FALSE, recursive = FALSE )
  
  for ( i in files ) {
    
    list <- list.files ( path = paste('~/Documents/rnaseq/datasets/', i, sep='' ), 
                         pattern = '*.fastq.gz$', all.files = FALSE,
                         full.names = TRUE, recursive = FALSE, ignore.case = FALSE, 
                         include.dirs = FALSE, no.. = FALSE )
    
    setwd('~/FastQC')
    system( paste( './fastqc -t 12', list[[1]], list[2], sep = ' ' ))
    
  }
}


fastqc()


###########################################################################################################################################



#TRIMMOMATIC ON LOW QUALITY DATASETS


qc <- function(){
  
  setwd ( '~/Documents/rnaseq/datasets' )
  files <- list.dirs ( path = ".", full.names = FALSE, recursive = FALSE )
  
  for ( i in files ) {
    
    setwd('~/Documents/rnaseq/trimmed')
    system(paste('mkdir', i, sep=' '))
    
  }
  
  trim <- function (i) {
    
    list <- list.files ( path = paste('~/Documents/rnaseq/datasets/', i, sep='' ), 
                         pattern = '*.fastq.gz$', all.files = FALSE,
                         full.names = TRUE, recursive = FALSE, ignore.case = FALSE, 
                         include.dirs = FALSE, no.. = FALSE )
    
    setwd ( '~/Trimmomatic-0.36' )
    system ( paste ( 'java -jar trimmomatic-0.36.jar PE -threads 12 -phred33 ', 
                     list[1], ' ', list[2], 
                     
                     ' ~/Documents/rnaseq/trimmed/', i, '/', i, '_forward_paired.fastq.gz ', 
                     ' ~/Documents/rnaseq/trimmed/', i, '/', i, '_forward_unpaired.fastq.gz ',
                     
                     '~/Documents/rnaseq/trimmed/', i, '/', i,'_reverse_paired.fastq.gz ', 
                     ' ~/Documents/rnaseq/trimmed/', i, '/', i, '_reverse_unpaired.fastq.gz ',
                     
                     'ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 ',
                     'TRAILING:3 MINLEN:36', sep='' ))
    
  }
  
  opt<- toupper(readline( 'Trim all datasets? (Y/N): '))
  
  if (opt == 'Y') {
    
    for ( i in files ) {
      
      trim(i)
      
    }
  }
  
  if (opt == 'N') {
    
    for ( i in files ) {
      option_1 <- toupper(readline( paste ('Trim dataset? (Y/N): ', i, ': ', sep='')))
      
      if (option_1 == 'Y') {
        setwd ( '~/Documents/rnaseq/tools/Trimmomatic-0.39' )
        trim(i)
        
      }
    }
  }
}


qc()

#RUN FASTQC AFTER TRIMMING
fastqc()




###########################################################################################################################################



#PREPARATION FOR ALIGNMENT 


prep_data <- function(){
  
  
  setwd ( '~/Documents/rnaseq/trimmed' )
  files <- list.dirs ( path = ".", full.names = FALSE, recursive = FALSE )
  
  for ( k in files ) {
    
    setwd ( paste('~/Documents/rnaseq/trimmed/', k, sep='' ))
    system( paste('rm -rf ~/Documents/rnaseq/trimmed/', k, '/*_unpaired*', sep='' ))
    
  }
}

prep_data()



###########################################################################################################################################



#RUN FASTQC ON BAM FILES ---> STAR 


fastqc <-function(file_path){
  
  setwd ( paste('~/Documents/', file_path, sep='' ))
  files <- list.dirs ( path = ".", full.names = FALSE, recursive = FALSE )
  
  for ( i in files ) {
    
    list <- list.files ( path = paste('~/Documents/', file_path, i, sep='' ), 
                         pattern = '*.bam$', all.files = FALSE,
                         full.names = TRUE, recursive = FALSE, ignore.case = FALSE, 
                         include.dirs = FALSE, no.. = FALSE )
    
    setwd('~/Documents/rnaseq/tools/FastQC')
    system( paste( './fastqc -t 12', list[[1]], sep = ' ' ))
    
  }
}


#HISAT2 
fastqc('hisat_bam_sorted/')


#STAR
fastqc('star/bam/')



###########################################################################################################################################



#BUILD GRCH38 GENOME ASSEMBLY/TRANSCRIPTS



build_hisat <- function(){
  
  setwd('~hisat2')
  system( './hisat2-build -p 12 ~/Documents/hisat2/assembly.fa genome')
  
}

build_star <- function(){
  
  setwd('~/STAR-2.7.8a/bin/Linux_x86_64')
  system( paste('./STAR --runThreadN 12 --runMode genomeGenerate --genomeDir ',
                '~/Documents/star --genomeFastaFiles ~/Documents/star/assembly.fa', sep='' ))
  
}

build_salmon <- function(){
  
  setwd('~/salmon/bin')
  system(paste('./salmon index -t ~/Documents/salmon/transcripts.fa ',
               '-i ~/Documents/salmon/transcripts_index -p 12', sep=''))
  
}


build_kallisto <- function(){
  
  setwd('~/kallisto/build/src')
  system(paste('./kallisto index -i ~/Documents/kallisto/genome/transcripts.idx ',
               '~/Documents/kallisto/transcripts.fa', sep=''))
  
}


build_bowtie2 <- function(){
  
  setwd('~/bowtie2')
  system( './bowtie2-build --threads 12 ~/Documents/star/assembly.fa ~/Documents/bowtie2/genome')
  
}


build_bwa <- function(){
  
  setwd('~/Documents/bwa')
  system('bwa index assembly.fa -b 1000000000000000')
  
}



#EXECUTION
build_hisat()
build_star()
build_salmon()
build_kallisto()
build_bowtie2()
build_bwa()



###########################################################################################################################################



#REFERENCE GENOME ALIGNMENT


align <- function(files, align_type){
  
  get_in <- function(opt){
    
    ll <- list.files( path = paste( '~/Documents/rnaseq/', opt, '/', i, sep='' ), 
                      pattern = '*fastq.gz$', all.files = FALSE, 
                      full.names = TRUE, recursive = FALSE, ignore.case = TRUE, 
                      include.dirs = FALSE, no.. = FALSE )
    return(ll)
    
  }
  
  setwd( '~/Documents/rnaseq/datasets')
  files <- list.dirs( path = ".", full.names = FALSE, recursive = FALSE )
  
  opt <- toupper(readline( 'Align all trimmed datasets (Y/N): ' ))
  
  if (opt == 'Y'){
    
    for( i in files ){
      
      files.in <- get_in('trimmed')
      align_type(i, files.in)
      
    }
  }
  
  if (opt == 'N'){
    
    for( i in files ){
      
      setwd( '~/Documents/rnaseq' )
      option <- toupper(readline( paste ('Dataset ', i, ' trimmed? (Y/N): ', i, ': ', sep='')))
      
      if (option == 'Y'){
        
        files.in <- get_in('trimmed')
        align_type(i, files.in)
        
      } else{
        
        files.in <- get_in('datasets')
        align_type(i, files.in)
        
      }
    }
  }
}


align_hisat <- function(i, files.in){
  
  setwd( '~/Documents/hisat2/')
  system( paste( 'hisat2 -p 12 -x genome -1 ', files.in[1], ' -2 ', files.in[2], 
                 ' -S ~/Documents/hisat2/bam/', i, '.bam | ', 
                 '~/samtools/bin/samtools sort -@ 12 -T temp -o ~/Documents/hisat/bam_sorted/', i, '.sorted.bam', sep='' ))
  
}


align_star <- function(i, files.in){
  
  
  setwd('~/STAR-2.7.8a/bin/Linux_x86_64')
  system( paste( './STAR --genomeDir /home/lebo/Documents/star --runThreadN 12 --readFilesIn ', 
                 files.in[1], ' ',files.in[2], ' --outFileNamePrefix ',
                 '/home/lebo/Documents/star/bam/', i, '/', i,
                 ' --outBAMsortingThreadN 12 --outSAMtype BAM SortedByCoordinate ',
                 '--limitBAMsortRAM 30000000000 --readFilesCommand zcat', sep=''))
  
}


align_salmon <- function(i, files.in){
  
  setwd('~/salmon/bin')
  system(paste('./salmon quant -p 12 -i ~/Documents/salmon/transcripts_index -l IU ',
               '-1 ', files.in[1], ' -2 ', files.in[2], ' --validateMappings ',
               '-o ~/Documents/salmon/transcripts_quant/', i, sep=''))
  
}


align_kallisto <- function(i, files.in){
  
  setwd(paste('~/Documents/kallisto/kallisto_quants/', i, sep=''))
  system(paste('~/kallisto/build/src/kallisto quant -t 12 -i ~/Documents/kallisto/genome/transcripts.idx ',
               '-b 30 -o kallisto_out --gtf genome.gtf ', files.in[1], ' ', files.in[2], sep=''))
  
}


align_bowtie <- function(i, files.in){
  
  setwd('~/bowtie2')
  system(paste('./bowtie2 -q -p 12 -x ~/Documents/bowtie2/genome -1 ', files.in[1], ' -2 ', files.in[2], ' | ', 
          '~/samtools/bin/samtools sort -@ 12 -T temp -o ~/Documents/bowtie2/bam/', i, '.sorted.bam', sep=''))
  
}



align_bwa <- function(i, files.in){
  
  system(paste('bwa mem -t 4 -R ~/Documents/bwa/genome ', files.in[1], ' ', files.in[2], ' | ',
               '~/samtools/bin/samtools sort -@ 12 -T temp -o ~/Documents/bwa/bam/', i, '.sorted.bam', sep=''))
}


#EXECUTION 
align(files, align_hisat)
align(files, align_star)
align(files, align_salmon)
align(files, align_kallisto)
align(files, align_bowtie)
align(files, align_bwa)


