function ARI = Calculate_Cluster_AdjustedRandIndex(c1,c2)
% Calculate the "Hubert & Arabie adjusted Rand index" of two partitions
% where c1,c2 are vectors listing the class membership

% Based on the version written by David Corney (D.Corney@cs.ucl.ac.uk)

C=Contingency(c1,c2);	%form contingency matrix

n=sum(sum(C));
nis=sum(sum(C,2).^2);		%sum of squares of sums of rows
njs=sum(sum(C,1).^2);		%sum of squares of sums of columns

t1=nchoosek(n,2);		%total number of pairs of entities
t2=sum(sum(C.^2));	%sum over rows & columnns of nij^2
t3=.5*(nis+njs);

%Expected index (for adjustment)
nc=(n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1));

A=t1+t2-t3;		%no. agreements


if t1==nc
   ARI=0;			%avoid division by zero; if k=1, define Rand = 0
else
   ARI=(A-nc)/(t1-nc);		%adjusted Rand - Hubert & Arabie 1985
end

end

function Cont=Contingency(Mem1,Mem2)
%CONTINGENCY Form contigency matrix for two vectors
% C=Contingency(Mem1,Mem2) returns contingency matrix for two
% column vectors Mem1, Mem2. These define which cluster each entity 
% has been assigned to.


if nargin < 2 | min(size(Mem1)) > 1 | min(size(Mem2)) > 1
   error('Contingency: Requires two vector arguments')
   return
end

Cont=zeros(max(Mem1),max(Mem2));

for i = 1:length(Mem1);
   Cont(Mem1(i),Mem2(i))=Cont(Mem1(i),Mem2(i))+1;
end

end


