function y=Clustering_kMeans_2dim()
clear all
close all
m=rand(100,2);
k=5;
y=kMeansCluster(m,k)
end

function y=kMeansCluster(m,k,isRand)

if nargin<3,        isRand=0;   end
if nargin<2,        k=1;        end
if nargin<1,        k=9;
    load('data\CAST_randData_2+500_inc.mat');
    m=D;
end


[maxRow, maxCol]=size(m)
if maxRow<=k,
    y=[m, 1:maxRow]
else
    
    % initial value of centroid
    if isRand,
        p = randperm(size(m,1));      % random initialization
        for i=1:k
            c(i,:)=m(p(i),:)
        end
    else
        for i=1:k
            c(i,:)=m(i,:)        % sequential initialization
        end
    end
    
    temp=zeros(maxRow,1);   % initialize as zero vector
    itr=1;
    while 1,
        itr=itr+1;
        d=DistMatrix(m,c);  % calculate objcets-centroid distances
        [z,g]=min(d,[],2);  % find group matrix g
        if g==temp,
            break;          % stop the iteration
        else
            temp=g;         % copy group matrix to temporary variable
        end
        for i=1:k
            f=find(g==i);
            if f            % only compute centroid if f is not empty
                c(i,:)=mean(m(find(g==i),:),1);
            end
        end
    end
    
    y=[m,g];
    
    hold off;
    cc=hsv(200);
    center=c;
    C=g;
    D=m;
    scatter(D(:,1),D(:,2),[],'b','c');
    if ~isempty(C)
        hold on;
        for j=1:max(C);
            plot(D((C==j),1),D((C==j),2),'g+--');
            hold on;
            scatter(D((C==j),1),D((C==j),2),[],cc(j*10,:),'filled');
            hold on;
            scatter(center(j,1),center(j,2),[],'rs','filled');
            hold on;
            text(center(j,1)+.01,center(j,2), strcat('  ',num2str(j)) ,'FontSize',8,'color','red');
        end
        for i=1:length(D);
            text(D(i,1)+.01,D(i,2), strcat('  ',num2str(i)) ,'FontSize',8);
        end
    end
end
end

function d=DistMatrix(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%
% DISTMATRIX return distance matrix between points in A=[x1 y1 ... w1] and in B=[x2 y2 ... w2]
% Copyright (c) 2005 by Kardi Teknomo,  http://people.revoledu.com/kardi/
%
% Numbers of rows (represent points) in A and B are not necessarily the same.
% It can be use for distance-in-a-slice (Spacing) or distance-between-slice (Headway),
%
% A and B must contain the same number of columns (represent variables of n dimensions),
% first column is the X coordinates, second column is the Y coordinates, and so on.
% The distance matrix is distance between points in A as rows
% and points in B as columns.
% example: Spacing= dist(A,A)
% Headway = dist(A,B), with hA ~= hB or hA=hB
%          A=[1 2 3; 4 5 6; 2 4 6; 1 2 3]; B=[4 5 1; 6 2 0]
%          dist(A,B)= [ 4.69   5.83;
%                       5.00   7.00;
%                       5.48   7.48;
%                       4.69   5.83]
%
%          dist(B,A)= [ 4.69   5.00     5.48    4.69;
%                       5.83   7.00     7.48    5.83]
%%%%%%%%%%%%%%%%%%%%%%%%%%%
[hA,wA]=size(A);
[hB,wB]=size(B);
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