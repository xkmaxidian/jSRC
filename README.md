# jSRC

Clustering scRNA-seq by joint sparse representation clustering (jSRC) Overview:

This is code to do jSRC: a flexible and accurate joint learning algorithm for clustering of single-cell RNA-sequencing data given in the "experiment" section of the paper:
Wenming Wu, Xiaoke Ma*. "jSRC: A Flexible and Accurate Joint Learning Algorithm for Clustering of Single-cell RNA-sequencing Data".

The coding here is a generalization of the algorithm given in the paper. jSRC is written in the MATLAB programming language. To use, please download the jSRC folder and follow the instructions provided in the “README.doc”.

Files:

main_jsrc.m - A script with a real scRNA-seq data to shows how to run the code.

Normalize_TPM.m- TPM normalization scRNA-seq count data X. Here, for scRNA-seq data, and normalize such that the total sum of counts, \sum_{i}X_{ij}, is 10^6 in each cell j, which is essentially the Transcript per Million (TPM) normalization.

%%It's worth noting that there are two main functions called jSRC, one in vector form and the other in matrix form. Call the matrix form when the data size is larger. When the matrix size is smaller, the two functions are used. The matrix function is shown here.
clustering_jsrcMatrix.m - The main function of jSRC to solve the variables by matrix solution form.

jsrc.m - The main function of jSRC to solve the variables by vector solution form.

Data_zheng.mat - A real scRNA-seq data used in the cell type clustering example. We retain the genes that are expressed in at least 10 cells for the dataset. The data Zheng dataset contains 500 human peripheral blood mononuclear cells (PBMCs) sequenced using GemCode platform, which consists of three cell types, CD56+ natural killer cells, CD19+ B cells and CD4+/CD25+ regulatory T cells.

Data

data$Zhengexpr.csv - Zheng dataset original expression data. The original data can be downloaded from 10X GENOMICS website. 

data$Zheng.celltype.csv - The cell type of Zheng dataset.

GSE51372_readCounts.txt - Ting dataset original expression data. 

Example:

Follow the steps below to run jSRC（also contained in the " main_jsrc.m" file）. Here use a real scRNA-seq data (Data_Zheng) set as an example.

%%% The input of jSRC is the data normalized by TPM to scRNA-seq count data.

%%%%====There are two main functions called jSRC, one in vector form and the other in matrix form. Call the matrix form when the data size is larger. When the matrix size is smaller, the two functions are used. The matrix function is shown here.

%%%%====Taking zheng dataset as an example

clear all;

clc;

load('Zheng_Data.mat') %%%Loading data

X=Data; %%%%%%%scRNA-seq count data, rows are genes, columns are cell sample

X=Normalize_TPM(X);%%%TPM normalization scRNA-seq count data X

k=5; %%%%%%%%%%%%%%k is the number of feature for dimension reduction,

[P,Z,A]=clustering_jsrcMatrix(X,k); %% Call the main function to solve the variables by matrix solution form

% [P,Z,A]=jsrc(X,k); %% Call the main function to solve the variables by vector solution form

%%% Clustering cell type: the block diagonal matrix R for the cells is constructed through strongest relationship based on matrix A, 

% %s % i.e.,  r{ij}^{*}=1, where j=argmax_{j} |a{.i}|, 0 otherwise. 

% % The similarity matrix built this way has ideally k2 connected components corresponding to the k2 clusters 

for i=1:size(X,2)
    
    r=zeros(1,size(X,2));
    
    ai=A(:,i);
    
    [s,t]=find(abs(ai)==max(abs(ai)));
    
    S(i)=s;
 
    if s>=i
        
        r(s+1)=1;
    
    else
       
       r(s)=1;
    
    end
    
    R(:,i)=r;

end

%%%%%%%%%%% Cluster grid map 

figure

spy(R)

Contact:

Please send any questions or found bugs to Xiaoke Ma xkma@xidian.edu.cn
