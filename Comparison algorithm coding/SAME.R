install.packages("devtools")
devtools::install_github("yycunc/SAMEclustering")
##https://github.com/yycunc/SAMEclustering
library("SAMEclustering")
filepath=file.choose()
mat=readMat(filepath)
Data=mat$data ############Loading data
dim(Data)
##data_SAME$Zheng.expr[1:5, 1:5]
system.time({
cluster.result <- individual_clustering(inputTags = data_SAME$Zheng.expr, mt_filter = TRUE,
percent_dropout = 10, SC3 = TRUE, CIDR = TRUE, nPC.cidr = NULL, Seurat = TRUE, nGene_filter = FALSE,
nPC.seurat = NULL, resolution = 0.7, tSNE = TRUE, dimensions = 2, perplexity = 30, SIMLR = TRUE, diverse = TRUE, SEED = 123)
cluster.ensemble <- SAMEclustering(Y = t(cluster.result), rep = 3, SEED = 123)
})
cluster.ensemble
