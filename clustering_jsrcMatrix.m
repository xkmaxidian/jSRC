function[P,Z,A]=clustering_jsrcMatrix(X,k)
%%%%%%% objective 
%          min{P,Z,A}||X-PZ||^2+
%          ||Z-ZA||^2+alpha*sum{i}||z.i-z.j||^2+beta*||A||1, s.t.
%          diag(A)=0, j=argmax_j|ai|1.
% You only need to provide the above two inputs.
% Notation:
% X ... (m x n) Normalized scRNA-seq data matrix 
%       m  ... number of  genes (feature size)   
%       n  ... number of cells    
% k ... number of feature for dimension reduction
%%%%Output:
%         P！！Project matrix
%         Z！！Feature matrix
%         A！！Coefficient matrix
% Coder Wenming Wu Email: wenmingwu55 at 163.com

%   version  -- July/2020, revision,--November/2020
%   maxIter The maximum number of iterations
%   tol1,tol2 Iteration error

%%%%%%%===========Initialization================
[d,n]=size(X);
[A1,B1,C1]=svds(X,k);
P=abs(A1*B1^0.5);
Z=abs(B1^0.5*C1');%%SVD initialization P and Z
A= rand(n,n);
B=A;
T= zeros(n,n);
I=eye(k,k);
II=eye(n,n);
%%%%%%===========Parameter settings===========
maxIter=100; %%Maximum number of iterations 
beta=0.6; %%Regularization parameter
alpha=0.2; %%Regularization parameter
tol1=1e-4;tol2=1e-4;
iter = 0; 
converged = 0;

%%%%%%%===========Update variables P,Z,A by iteration================
while ~converged  && iter < maxIter   
iter = iter + 1;
derta =50; 
%%%%%===========Update variable P ===========
P=P.*((X*Z')./(P*Z*Z'));

%%%%%===========Update variable Z ===========
for i=1:size(Z,2)
    ai=A(:,i);
    [j,v]=find(abs(ai)==max(abs(ai)));
    Zj(:,i)=Z(:,j);
end
Z=Z.*((P'*X+Z*A+alpha*Zj)./(P'*P*Z+Z+alpha*Z));

%%%%%===========Update variable A by ADMM===========
A=A.*((Z'*Z+derta*B-T)./(Z'*Z*A+derta*A));
for i=1:size(A,2)
    A(i,i)=0;
end
B=soft(A+T/derta,beta/derta);
%%%%%===========Update Lagrange multiplier T ===========
T=T+1.1*(A-B);

%%%%%%%%%%%%%%%===========Error===========
% Err(iter,:)=abs(y-awkl*Hwkl(:,1:size(X,2)))/20;
% Er(iter,:)=abs(mean(V-Wwkl*Hwkl));
% A(iter,:)=mean(Hwkl);

%%%%%%%%%%%%%%%===========Error===========
% temp = max ([norm(A-awk,2),norm(bk-bwk,2),norm(Wk-Wwk,'fro'),norm(Hk-Hwk,'fro')]);
% %     temp = muu*temp/norm(V,2);
% temp =temp/max([norm(ak,2),norm(bk,2),norm(Wk,'fro'),norm(Hk,'fro')]);
%     temp = max([(sqrt(L)*norm(ZK-Zkm1,'fro')),norm(WK-Wkm1,'fro'),norm(EK-Ekm1,'fro')]);
%     temp = muu*temp/norm(Y,'fro');
%     
%%%%%%%%%%%%%%%===========Error===========
temp1 = (norm((X - P*Z),'fro')+norm((Z -Z*A),'fro'))/sqrt(norm(X,'fro'));
% temp1 = (norm((X - P*Z),'fro')+norm((Z - xZ),'fro')) /max(norm(X,'fro'),norm(Z,'fro'))^2;
%     temp1 = norm( (V - WK*HK ), 'fro') / max([norm(WK*HK , 'fro'),norm(V,'fro')]);
if temp1 < tol1 
converged = 1;
end
%     disp(['temp1 ',num2str(temp1)]);
%     disp([' 亨旗肝方 ' num2str(iter) ' temp1 ' num2str(temp1) ]);
% t1(iter)=temp1;
end

end


%%%%%=========== Soft thresholding function called when solving B===========
function[y] = soft( x, T )
if sum( abs(T(:)) )==0
   y = x;
else
   y = max( abs(x) - T, 0);
   y = sign(x).*y;
end
end    