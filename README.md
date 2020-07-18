# jSRC
Clustering scRNA-seq by joint learning sparse representation-based clustering Overview:

This is code to do jSRC: a flexible and accurate joint learning algorithm for clustering of single-cell RNA-sequencing data given in the "experiment" section of the paper:
Wenming Wu, Xiaoke Ma*. "jSRC: A Flexible and Accurate Joint Learning Algorithm for Clustering of Single-cell RNA-sequencing Data".

The coding here is a generalization of the algorithm given in the paper. jSRC is written in the MATLAB programming language. To use, please download the jSRC folder and follow the instructions provided in the “README.doc”.

Files:

jsrc.m - The main function.

main_jsrc.m - A script with a real scRNA-seq data to shows how to run the code.

Data_zheng.mat - A real scRNA-seq data used in the cell type clustering example. We retain the genes that are expressed in at least 10 cells for the dataset. The data Zheng dataset contains 500 human peripheral blood mononuclear cells (PBMCs) sequenced using GemCode platform, which consists of three cell types, CD56+ natural killer cells, CD19+ B cells and CD4+/CD25+ regulatory T cells.


data$Zhengexpr.csv - Zheng dataset original expression data. The original data can be downloaded from 10X GENOMICS website. 

data$Zheng.celltype.csv - The cell type of Zheng dataset.

GSE51372_readCounts.txt - Ting dataset original expression data. 

Example:

Follow the steps below to run jSRC（also contained in the " main_jsrc.m" file）. Here use a real scRNA-seq data (Data_Zheng) set as an example.

clear all;

clc;

load('Data_zheng.mat') %%%Loading data

X=data; %%%%%%%scRNA-seq data, rows are genes, columns are cell sample

k=5; %%%%%%%%%%%%%%k is the number of feature for dimension reduction,

[P,Z,A]=jsrc(X,k); %% Call the main function to solve the variables 

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

Please send any questions or found bugs to Xiaoke Ma xkma@xidian.edu.cn
