clear all;
clc;
load('Zheng_Data.mat') %%%Loading data
X=Data; %%%%%%%scRNA-seq count data, rows are genes, columns are cell sample
% X=Normalize_TPM(X);%%%TPM normalization scRNA-seq count data X
k=5; %%%%%%%%%%%%%%k is the number of feature for dimension reduction,
k2=3;%%k is the number of clusters,
[mFea,nSmp]=size(X);
Winit = abs(rand(mFea,k));
Hinit = abs(rand(nSmp,k))'; 
maxiter=100;
[P,Z] = DR(X,Winit,Hinit,maxiter); %% Call the main function to solve the variables by matrix solution form
l=kmeans(Z',k2);

    