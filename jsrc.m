function[P,Z,A]=jsrc(X,k)

%%%%The model----------------
% min{P,Z,A}||X-PZ||^2+ ||Z-ZA||^2+alpha*sum{i}||z.i-z.j||^2+beta*||A||1
% s.t. j=argmaxj |aji|, diag(A)=0, i=1,2,...,n

% You only need to provide the above two inputs.
% Notation:
% X ... (m x n) scRNA-seq data matrix 
%       m  ... number of  genes (feature size)   
%       n  ... number of cells               
% k ... number of feature for dimension reduction

%%%%Output:
%         P！！Project matrix
%         Z！！Feature matrix
%         A！！Coefficient matrix
% Coder Wenming Wu Email: wenmingwu55 at 163.com

%   version  -- July/2020, 

%   maxIter The maximum number of iterations
%   tol1,tol2 Iteration error

%%%%%%%===========Initialization================

[m,n]=size(X);
P=rand(m,k); Z=rand(k,n);
A= rand(n-1,n);
B=A;
T= zeros(n-1,n);
I=eye(k,k);
II=eye(n-1,n-1);

%%%%%%===========Parameter settings===========
maxIter=80;  
tol1=1e-3;tol2=1e-3;
iter = 0; 
converged = 0;

%%%%%%%===========Update variables P,Z,A by iteration================
while ~converged  && iter < maxIter   

iter = iter + 1;
beta=0.6; %%%%%%%Regularization parameter
alpha=0.2; %%%%%%%Regularization parameter
derta =1; 
%%%%%===========Update variable P ===========
P=P.*((X*Z')./(P*Z*Z'));
%%%%%===========Update variable Z ===========
for i=1:size(Z,2)
    Sz=Z; %%%%%%%Sz, dictionary (columns correspond to atoms) in sparse representation
    zi=Sz(:,i);
    Sz(:,i)=[];
    xi=X(:,i);
    ai=A(:,i);
    [j,v]=find(abs(ai)==max(abs(ai)));
    zj=Sz(:,j);
    zi=zi.*((P'*xi+Sz*ai+alpha*zj)./((P'*P+I+alpha*I)*zi));
    Z(:,i)=zi;
end
%%%%%===========Update variable A by ADMM, where A=[a.1,a.2,...,a.n]===========
for i=1:size(A,2)
    Sz=Z;
    zi=Sz(:,i);
    Sz(:,i)=[];
    ai=A(:,i);
    bi=B(:,i);
    ti=T(:,i);
    ai=ai.*((Sz'*zi+derta*bi-ti)./((Sz'*Sz+derta*II)*ai));
    bi=soft(ai+ti/derta,beta/derta);
    A(:,i)=ai;
    B(:,i)=bi;
end
%%%%%===========Update Lagrange multiplier T ===========
T=T+1.618*derta*(A-B);

%%%%%%%%%%%%%%%===========Error===========
% Err(iter,:)=abs(y-awkl*Hwkl(:,1:size(X,2)))/20;
% Er(iter,:)=abs(mean(V-Wwkl*Hwkl));
% A(iter,:)=mean(Hwkl);

%%%%%%%%%%%%%%===========Error===========
% temp = max ([norm(A-awk,2),norm(bk-bwk,2),norm(Wk-Wwk,'fro'),norm(Hk-Hwk,'fro')]);
% %     temp = muu*temp/norm(V,2);
% temp =temp/max([norm(ak,2),norm(bk,2),norm(Wk,'fro'),norm(Hk,'fro')]);
%     temp = max([(sqrt(L)*norm(ZK-Zkm1,'fro')),norm(WK-Wkm1,'fro'),norm(EK-Ekm1,'fro')]);
%     temp = muu*temp/norm(Y,'fro');
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%===========Error===========
for i=1:size(A,2)
    Sz=Z;
    zi=Sz(:,i);
    Sz(:,i)=[];
    xai=A(:,i);
    xZ(:,i)=Sz*xai;
end
temp1 = (norm((X - P*Z),'fro')+norm((Z -xZ),'fro'))/norm(X,'fro');
if temp1 < tol1 
converged = 1;
end
%     disp(['temp1 ',num2str(temp1)]);
%     disp([' Number of iterations ' num2str(iter) ' temp1 ' num2str(temp1) ' temp ' num2str(temp)]);
t1(iter)=temp1;
end
        t1=(t1-min(t1))./(t1(1)-min(t1));
%         t=1:iter;
%         figure
%         plot(t,t1,'b-'),xlabel('Iteration times');ylabel('Error');
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