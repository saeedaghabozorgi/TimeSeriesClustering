function cluster = do_Hierarchical(Distance, Method, k)
%DCAGG Performs agglomerative clustering
% Cluster = DCAGG(Distance, Method, k) where Distance is
% square dissimilarity matrix, with Inf on leading diagonal,
% Method is one of 'single', 'complete' or 'centroid',
% and k is number of clusters. Cluster is a cell array showing
% which entities belong to which cluster.

%(C) David Corney (2000)   		D.Corney@cs.ucl.ac.uk

%If you want a nice hierarchical graph, try Jang's hierarchical
%clustering: http://neural.cs.nthu.edu.tw/jang/mlbook/demo/


%demo
% D=[0 0;1 1; 2 1; 0 4; 2 5; 3 5; 3 4; 4 5; 4 4];
% Distance=dis_euclidean_matrix(D,D);
% Distance(1:10:end)=Inf;
% Method='single';
% k=3;
%---

% if nargin < 3
%    error('DCAGG requires three arguments');
%    help(mfilename)
%    return
% end

[R,C] = size(Distance);

cluster=num2cell(1:R);		%Put every point is in its own cluster

for i = 2:(R-k+1)
   
   %Find closest clusters
   [MinRow, IdxDow] = min(Distance);
   [temp, MinJ] = min(MinRow);
   MinI = IdxDow(MinJ);
   
	if MinI > MinJ
      t=MinI;
      MinI=MinJ;
      MinJ=t;
   end
   
   
   %Merge cluster j into cluster i, then delete j
   cluster{MinI} = [cluster{MinI} cluster{MinJ}];
   cluster(MinJ) = [];
   
   %Calculate new Distance matrix
   switch Method
   case 'single'
      Distance(:, MinI) = min(Distance(:, MinI), Distance(:, MinJ)); 
      Distance(MinI, :) = min(Distance(MinI, :), Distance(MinJ, :)); 
   case 'complete'
      Distance(:, MinI) = max(Distance(:, MinI), Distance(:, MinJ)); 
      Distance(MinI, :) = max(Distance(MinI, :), Distance(MinJ, :)); 
   case 'centroid'
      Distance(:, MinI) = (Distance(:, MinI) + Distance(:, MinJ))/2; 
      Distance(MinI, :) = (Distance(MinI, :) + Distance(MinJ, :))/2;
   otherwise
      error('Unsupported Method in DCAGG!');
   end
   
   Distance(MinJ, :) = [];
   Distance(:, MinJ) = [];
   Distance(MinI, MinI) = inf;
   
end


