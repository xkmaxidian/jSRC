%%%%%%%A script to shows how to run DR+SR.
clear all;
clc;
load('Data_zheng.mat')
X=data;
% X=Normalize_TPM(X);%%%TPM normalization scRNA-seq count data X
k=5; %%%%%%%%%%%%%%k is the number of feature for dimension reduction,
k2=3;%%k is the number of clusters,
[mFea,nSmp]=size(X);
Winit = abs(rand(mFea,k));
Hinit = abs(rand(nSmp,k))'; 
maxiter=100;
[P,Z] = DR(X,Winit,Hinit,maxiter); %% Call the main function to solve the variables by matrix solution form
lambda=0.1;
[A] = SR(Z, lambda);

 for i=1:size(A,2)
    r=zeros(1,500);
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
figure
spy(R)