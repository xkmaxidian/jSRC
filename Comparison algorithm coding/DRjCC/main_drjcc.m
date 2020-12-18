    clear all;
    clc;
    load('Data_zheng.mat')%Loading data
     X=data;%%%%%%%Rows are genes, columns are cell sample
    %==============Constructing a weight matrix==============
    %Preset value before constructing weight matrix
    options = [];
    option.Metric = 'Cosine';
    options.NeighborMode = 'KNN';%KNN
    options.k =5;%5 nearest neighbors
    options.WeightMode = 'Cosine';%Weights are 0 or 1, it can eplace with 'HeatKernel', 'Euclidean' 

    W = constructW(X',options);

    clear options;
    options = [];
    k1=100; k2=3;
    %%%%%%%%%%%%%%k1 is the number of feature for dimension reduction,
    %%%%%%%%%%%%%%k2 is the number of features for clustering in NMF(the number of cluster)
    [Z_final,A_final,B_final, F_final] = drjcc(X, k1, k2, W, options); %% Call the main function to solve the variables
%%%%%%%%%%% Clustering cell type label
    for e=1:size(F_final,2)
    v=F_final(:,e);
    ma=max(v);
    [s,t]=find(v==ma);
    l(e)=s;
    end
    %%%%%%%%%%%%%%==================Performance evaluation===============================
 