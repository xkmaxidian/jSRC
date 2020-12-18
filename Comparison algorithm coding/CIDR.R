
library(cidr)###https://github.com/VCCRI/CIDR
filepath=file.choose()
mat=readMat(filepath)
Data=mat$data ############Loading data
tags <- as.matrix(t(Data))
tags<-t(tags)
dim(Data)
 ##cols <- c(rep("RED",k), rep("BLUE",k), rep("GREEN",k))
 ## Standard principal component analysis.
## ltpm <- log2(t(t(tags)/colSums(tags))*1000000+1)
 ##pca <- prcomp(t(ltpm))
##plot(pca$x[,c(1,2)],col=cols,pch=1,xlab="PC1",ylab="PC2",main="prcomp")
 ## Use cidr to analyse the simulated dataset.
 ## The input for cidr should be a tag matrix.
system.time({
 sData <- scDataConstructor(tags)
 sData <- determineDropoutCandidates(sData)
 sData <- wThreshold(sData)
 sData <- scDissim(sData)
 sData <- scPCA(sData)
 sData <- nPC(sData)
 nCluster(sData)
 sData <- scCluster(sData)
})

 Use Adjusted Rand Index to measure the accuracy of the clustering output by cidr.
adjustedRandIndex(sData@clusters,cols)
