function NX=Normalize_TPM(X)
%%%=========Here, for scRNA-seq data, and normalize such that the total sum of counts, X,
%%%====is 10^6 in each cell j, which is essentially the Transcript per Million (TPM) normalization.
        for i=1:size(X,1)
            x=X(i,:);
            n=find(x~=0);
            l=length(n);
            x=x./l;
           xX(i,:)=x;
        end
       ll=sum(xX);
        NX=xX./ll;
        NX=NX*1e6;
end