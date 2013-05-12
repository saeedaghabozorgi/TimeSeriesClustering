function [Classes,Centres,FinalDistance]=do_kMedoids(Data,k,InitCentres,MaxIters)

if nargin < 3
    InitCentres = ChooseInitialCentres(Data,k);		%randomly choose starting point (where needed)
end
Centres=InitCentres;
OldCentres=Centres;
if nargin<4
    MaxIters=500;
end

[R,C]=size(Data);

DataSq=repmat(sum(Data.^2,2),1,k);	%sum squared data - save re-calculating repeatedly later
%Do we need DataSq? It's constant, and we're minimsing things...


for i = 1:MaxIters
    
    Dist = DataSq + repmat(sum((Centres.^2)',1),R,1) - 2.*(Data*(Centres'));   %i.e. d^2 = (x-c)^2 = x^2 + c^2 -2xc
    
    [D,Centre]=min(Dist,[],2);		%label of nearest centre for each point
    
    for j=1:k
        idx=find(Centre==j);
        if length(idx)>0
            Centres(j,:)=medoid(Data(idx,:));
        end
    end
    clf;
    scatter(Data(:,1),Data(:,2),[],'b','c');
    hold on;
    scatter(Centres(:,1),Centres(:,2),[],'b','filled');
    
    Change=sum(sum(abs(OldCentres-Centres)));
    if Change < 1e-10	%Have we converged yet?
        break
    end
    OldCentres=Centres;
    
end
[FinalDistance,Classes]=min(Dist,[],2);		%label points one last time
end

function [medoid]=medoid(Data)
dis=dis_euclidean_matrix(Data,Data);
dis=dis^2;
Error=sum(dis);
[s,idx]=min(Error);
medoid=Data(idx,:);
end

