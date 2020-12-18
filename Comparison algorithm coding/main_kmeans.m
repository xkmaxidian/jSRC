clear all;
clc;
load('Data_zheng.mat')
X=data;
% X=Normalize_TPM(X);%%%TPM normalization scRNA-seq count data X
k2=3;
l=kmeans(X',k2);