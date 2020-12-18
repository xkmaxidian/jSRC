library(umap)
filepath=file.choose()
mat=readMat(filepath)
Data=mat$data ############Loading data

system.time({
custom.settings = umap.defaults
custom.settings$n_components =50
custom.settings
iris.umap = umap(t(Data),config=custom.settings)
})
dim(iris.umap$layout)

###Import the features "iris.umap$layout" extracted by UMAP into kmeans to complete clustering
