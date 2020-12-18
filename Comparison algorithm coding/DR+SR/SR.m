function [A] = SR(X, lambda)
Dic = X;
flag = sqrt(sum(Dic.^2,1));
Dic = Dic./repmat(flag,size(Dic,1),1); % 字典原子归一化
for i = 1:size(X,2)
    xi = Dic(:, i);
    B=Dic;
    B(:,i)=[];
  [A(:, i), ~] = l1_ls(B, xi, lambda);%%%%%%Solving the sparse representation coefficient matrix A; 
%   Download "l1_ls.m" by https://github.com/cvxgrp/l1_ls 
end
end
