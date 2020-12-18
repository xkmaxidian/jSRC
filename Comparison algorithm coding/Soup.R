
library("devtools")
devtools::install_github("lingxuez/SOUP")
####https://github.com/lingxuez/SOUPR
library(SOUP)
filepath=file.choose()
mat=readMat(filepath)
Data=mat$data ############Loading data
counts=t(Data)
dim(counts)
log.expr = log2(scaleRowSums(counts)*(10^6) + 1)
dim(log.expr)

select.out = selectGenes(counts, type="count", n.cores=10)
select.genes = select.out$select.genes
select.genes = camp$select.genes
log.select.expr = log.expr[, colnames(log.expr) %in% select.genes]
dim(log.select.expr)
dim(log.select.expr)

system.time({
soup.out = SOUP(counts, Ks=3, type="log")
})

length(soup.out$memberships)
length(soup.out$centers)
 table(soup.out$major.labels)


