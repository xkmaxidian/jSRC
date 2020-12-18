install.packages("devtools")
install_github("shibiaowan/SHARP")
library(SHARP)
filepath=file.choose()
mat=readMat(filepath)
scExp=mat$data ############Loading data
dim(scExp)
res = SHARP(scExp)
######More details: https://github.com/shibiaowan/SHARP
