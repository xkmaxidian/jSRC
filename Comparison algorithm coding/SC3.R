


library(Seurat)
library(SingleCellExperiment)
library(SC3)### http://bioconductor.org/packages/SC3/
library(ggplot2)
library(scater)

filepath=file.choose()
mat=readMat(filepath)
Data=mat$data ############Loading data
D<- Data
sce <- SingleCellExperiment(
    assays = list(
        counts = as.matrix(D),
        logcounts = log2(as.matrix(D) + 1)
    ), 
    colData = t(Da)
)

# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)
# remove features with duplicated names
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]

# define spike-ins
isSpike(sce, "ERCC") <- grepl("ERCC", rowData(sce)$feature_symbol)

sce <- sc3(sce, ks = 10, biology = TRUE)


sce <- sc3_kmeans(sce, ks = 3)
names(metadata(sce)$sc3$kmeans)
col_data <- colData(sce)
col_data[ , grep("sc3_", colnames(col_data))]

#sce <- sc3_calc_consens(sce)
#names(metadata(sce)$sc3$consensus)
#names(metadata(sce)$sc3$consensus$`3`)

#metadata(sce)$sc3$kmeans
#col_data <- colData(sce)
#col_data[ , grep("sc3_", colnames(col_data))]



