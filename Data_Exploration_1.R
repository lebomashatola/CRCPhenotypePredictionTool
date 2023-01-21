
system('export LANG=en_US.UTF-8')


cbind.fill = function (...){

  nm = list(...)
  nm = lapply(nm, as.matrix)

  n = max(sapply(nm, nrow))
  do.call(cbind, lapply(nm, function (x)
    rbind(x, matrix(, n-nrow(x), ncol(x)))))

}


unzip_files = function(paths){

  setwd( paths )

  control = list.dirs ( path = ".", full.names = F, recursive = F )
  df_meta = data.frame()

  for ( file_1 in control ){

    setwd( paste(paths, file_1, sep='' ))

    files = list.files(path = ".", recursive = FALSE,
                        pattern = "\\.gz$",
                        full.names = FALSE)

    for ( file in files ){

      data = read.table(gzfile(file[[1]]))
      df_meta = cbind.fill( df_meta, data[,2])
    }

    rownames(df_meta) = data$V1
    all_genes = as.character(rownames(df_meta))
    all_genes = sub("\\.\\d+", "", all_genes)

    rownames(df_meta) = all_genes

  }

  df_meta = as.data.frame(df_meta)

  return(df_meta)

}

df_resis = unzip_files('/Users/lebomash/Documents/PhD_2022/Counts/Resistant/')
df_sensi = unzip_files('/Users/lebomash/Documents/PhD_2022/Counts/Sensitive/')



prep = function(vars){

  library(dplyr)

  all = vars %>% filter(row_number() <= n()-5)

  return(all)

}

df_resis = prep(df_resis)
df_sensi = prep(df_sensi)


write.csv(df_resis, '/Users/lebomash/Documents/PhD_2022/Counts/Resistant.csv')
write.csv(df_sensi, '/Users/lebomash/Documents/PhD_2022/Counts/Sensitive.csv')


gene_sets = function(){

  files = list.files(path = '/Users/lebomash/documents/phd_2022/counts/count_features',
                     pattern = '\\.txt$', full.names = TRUE)

  df_all = data.frame()

  for (i in files){

    df = read.delim(i, header=FALSE)
    df = as.data.frame(df)
    df = df[-(1:2), , drop = FALSE]
    df_all = rbind(df_all, df)

  }

  df_all = unique(df_all)
  df_all = df_all[,1]

  return(df_all)

}

genes = gene_sets()


gene_name = function(x){

  library(clusterProfiler)
  library(org.Hs.eg.db)

  eg = bitr(x, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db")
  select_genes = eg$ENSEMBL

  write.csv(select_genes, '/Users/lebomash/Documents/PhD_2022/Counts/Processed_Files/Biomarkers.csv')

}

biomarkers = gene_name(genes)


gene_sets = function(){

  files = list.files(path = '/Users/lebomash/documents/phd_2022/counts/count_features',
                     pattern = '\\.txt$', full.names = TRUE)

  df_all = data.frame()

  for (i in files){

    df = read.delim(i, header=FALSE)
    df = as.data.frame(df)
    df = df[-(1:2), , drop = FALSE]
    df_all = rbind(df_all, df)

  }

  df_all = unique(df_all)
  df_all = df_all[,1]

  return(df_all)

}

genes = gene_sets()


gene_name = function(x, dataset){

  library(clusterProfiler)
  library(org.Hs.eg.db)

  eg = bitr(x, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db")
  select_genes = eg$ENSEMBL

  dataset = dataset[rownames(dataset) %in% select_genes ,]

}

df_R = gene_name(genes, df_R)
df_S = gene_name(genes, df_S)
