function [P,Z] = DR(X,Winit,Hinit,maxiter)
%%%%%%% objective 
%          min{P,Z}||X-PZ||^2
% You only need to provide the above three inputs.
% Notation:
% X ... (m x n) Normalized scRNA-seq data matrix 
%       Winit  ... Initial value of P    
%       Hinit  ... Initial value of Z   
% maxiter ... The maximum number of iterations
%%%%Output:
%         P！！Project matrix
%         Z！！Feature matrix
tol1=1e-5;
iter = 0; 
converged = 0;
P = Winit; Z = Hinit; 
while ~converged  && iter < maxiter 
    iter = iter + 1;
    %%%%%%%===========Update variables P,Z by iteration================
    P=P.*((X*Z')./(P*Z*Z'));
    Z=Z.*((P'*X)./(P'*P*Z));
    
    %%%%%%%%%%%%%%%===========Error===========
    temp1 = norm((X - P*Z),'fro');
    if temp1 < tol1 
        converged = 1;
    end
end
end

