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


