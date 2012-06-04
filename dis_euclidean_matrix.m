function d=dis_euclidean_matrix(A,B)
       [hA,wA]=size(A);
       [hB,wB]=size(B);
       
       %check for same dimensions
             if wA ~= wB,  error(' second dimension of A and B must be the same'); end
             for k=1:wA
                  C{k}= repmat(A(:,k),1,hB);
                  D{k}= repmat(B(:,k),1,hA);
             end
             S=zeros(hA,hB);
             for k=1:wA
                  S=S+(C{k}-D{k}').^2;
             end
             d=sqrt(S);
end