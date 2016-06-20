library(pcaMethods)
n_factor <- 40  # TODO
for (k in 0:(n_factor-1)){
  print(paste("Working on factor #", k))
  data = read.table(paste("./result/temp/f", toString(k), "_tissue_indiv.txt", sep=""))
  data[data=="NaN"] <- NA
  pc <- pca(data, nPcs=1, method="ppca")  # TODO
  loading <- loadings(pc)
  score <- scores(pc)
  # Individual
  write.table(loading, file=paste("./result/temp/f", toString(k), "_indiv.txt", sep=""), sep="\t", col.names = FALSE, row.names = FALSE)
  # Tissue
  write.table(score, file=paste("./result/temp/f", toString(k), "_tissue.txt", sep=""), sep="\t", col.names = FALSE, row.names = FALSE)
}
