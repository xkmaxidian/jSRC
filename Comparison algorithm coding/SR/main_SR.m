clear all;
clc;
load('Data_zheng.mat')
X=data;
% X=Normalize_TPM(X);%%%TPM normalization scRNA-seq count data X
lambda=0.1;
[A] = SR(X, lambda);

 for i=1:size(A,2)
    r=zeros(1,size(A,2));
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