



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

  z_scores = as.data.frame(sapply(all, function(all) (abs(all-mean(all))/sd(all))))
  rownames(z_scores) = rownames(all)
  no_outliers <- z_scores[!rowSums(z_scores>3), ]

  return(no_outliers)

}

df_resistance = prep(df_resis)
df_sensitive = prep(df_sensi)



plot_graph = function(dataset, opt){

  pdf(file = paste("/Users/lebomash/Documents/PhD_2022/Results/", opt, '.pdf',sep=''))

  affy::plotDensity(dataset, col='black',
                    lty=c(1:ncol(dataset)), xlab='Count',
                    main='Expression Distribution')

  dev.off()

}


plot_graph(df_resistance, 'Resistant_Origi')
plot_graph(df_sensitive, 'Sensitive_Origi')

plot_graph(df_res, 'Resistant_Log')
plot_graph(df_sen, 'Sensitive_Log')
