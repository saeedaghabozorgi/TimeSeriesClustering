function clustering_Hybrid_3Level_2dim(Data,k)

load('..\data\CAST_randData_10_1.mat');
%D=rand(100,2);
%D=[1 1; 2 1; 1 2; 2 2; 3 2; 3 3;2 4; 4 4;4 1];

clf;
scatter(D(:,1),D(:,2),[],'b','c');

% D=[];
% hFig = figure(1);
% axis([0 2 0 2])
% hold off
% for i=1:20
%     [x,y] = ginput(1);
%     D=[D;x,y];
%     hold on;
%     scatter(D(:,1),D(:,2),[],'b','c');
% end
% save('..\data\CAST_randData_10_3.mat', 'D');


cc=hsv(100);
c=[];
k=4;
[c,Centres,FinalDistance]=do_kMeans(D,k);
%[c,Centres,FinalDistance]=do_kMedoids(D,k);
for Cluster_num=1:k
    hold on;
    scatter(D(c==Cluster_num,1),D(c==Cluster_num,2),[],cc(Cluster_num*5,:),'filled');
    %   pause;
%     vv=(c==Cluster_num);
%     for i=1:max(vv);
%         plot(D((vv==i),1),D((vv==i),2),'g+--');
%     end
end

center=[];



%----------Level 2--CAST(raw)--------------------------------------------
clusterCount=size(unique((c)),1);

T=0;
for i=1:clusterCount
    newData=[];
    newData=find(c(:,1)==i);
    dis=dis_euclidean_matrix(D(newData,:),D(newData,:));
    dis=squareform(dis);
    sim=1./(1+dis);
    meansim(i,:)=[mean(sim) size(newData,1)];
end
T=sum((meansim(:,1).*meansim(:,2)))./size(D,1);
center=[];
for i=1:clusterCount

    newData=[];
    newData=find(c(:,1)==i);
    scatter(D(newData,1),D(newData,2),[],'b','filled');
    temp_c=[];
    temp_c= do_ECAST_2dim(D(newData,:),-1);
    clustCoun=max(temp_c);
    c(newData,2)=temp_c;
    data=D(newData,:);
    for Cluster_num=1:clustCoun
        hold on;
        mems=[];
        mems=find(temp_c==Cluster_num);
        if ~isempty(mems)
            scatter(data(mems,1),data(mems,2),[],cc(Cluster_num*10,:),'filled');
        end
        %   pause;
        %             vv=temp_c==Cluster_num;
        %             for i=1:max(vv);
        %                 plot(data((vv==i),1),data((vv==i),2),'g+--');
        %             end
    end
    
end
%-------
% to map raw objects to new clusters
c(:,3)=c(:,1)*1000+c(:,2);
[x,y]=sort(c);
clsNum=1;
c(y(1,3),4)=clsNum;
for i=2:size(c,1)
    if ( x(i,3) ~= x(i-1,3))
        clsNum=clsNum+1;
        c(y(i,3),4)=clsNum;
    else
        c(y(i,3),4)=clsNum;
    end
end

l2_clusterCount=max(c(:,4));

% to make the prototypes
prototypes_clus=[];
center=[];
newpoints=[];
Prot_2lcluster=[];
counter=1;
% clf;
% scatter(D(:,1),D(:,2),[],'b');
minp =4;
for i=1:l2_clusterCount;
    mems=find(c(:,4)==i);
    % find rep
    if(size(mems,1)>minp)
        dis=dis_euclidean_matrix(D(mems,:),D(mems,:));
        representative=rep(dis,round(size(mems,1)/minp));
        prototypes=D(mems(representative),:);
        prototypes=shrink_tocenter(D(mems,:),prototypes,0.2);
    elseif (size(mems,1)>0) 
        prototypes=mean(D(mems,:),1);
    end
        prototypes_clus{counter}=prototypes;
        newpoints=[newpoints; prototypes_clus{counter}];
        Prot_2lcluster=[Prot_2lcluster; repmat([i],size(prototypes_clus{counter},1),1)];
        hold on;
        scatter(prototypes_clus{counter}(:,1),prototypes_clus{counter}(:,2),[],cc(50,:),'filled');
        counter=counter+1;
end
%--------------------------------------
%remove outliers
for i=1:length(prototypes_clus)
   if size( prototypes_clus{i},1)==1
      newpoints (Prot_2lcluster==i,:)=[];
      Prot_2lcluster(Prot_2lcluster==i)=[];
      prototypes_clus{i}=[];
   end
end


%--Level 3-------------------------------------------------

Distance=dis_euclidean_matrix(newpoints,newpoints);
Distance(1:size(newpoints,1)+1:end)=Inf;


for i=1:size(newpoints,1)-1
    for j=i+1:size(newpoints,1)
        if(Prot_2lcluster(i)==Prot_2lcluster(j))
            Distance(i,j)=0;
            Distance(j,i)=0;
        end
    end
end

cluster = do_Hierarchical(Distance, 'single', k)

for i=1:k
    for j=1:size(cluster{i},2)
        prot_num =cluster{i}(j);
        clus2=Prot_2lcluster(prot_num);
        pp=find(c(:,4)==clus2);
        
    
        
        c(pp,5)=i;
    end
end

% hold on;
% scatter(centre(:,1),centre(:,2),[],cc(50,:),'filled');
% if l2_clusterCount>k
%   %  [c2,Centres,FinalDistance]=do_kMeans(centre,k)
%   c2=  do_Hierarchical(Distance, 'single', k);
%     for j=1:k
%         aa=find(c2==j);
%         for i=1:length(aa)
%             pp=find(c(:,4)==aa(i));
%             c(pp,5)=j;
%         end
%     end
% else
%     c(:,5)=c(:,4);
% end
clf;
scatter(D(:,1),D(:,2),[],'b','c');

final_c=c(:,5);
for Cluster_num=1:max(c(:,5));
    hold on;
    scatter(D(final_c==Cluster_num,1),D(final_c==Cluster_num,2),[],cc(Cluster_num*10,:),'filled');
    %   pause;
    %     vv=final_c==Cluster_num;
    %     for i=1:max(vv);
    %         plot(D((vv==i),1),D((vv==i),2),'g+--');
    %     end
end

end

function representative=rep(dis,num)
a=zeros(size(dis,1),1);
representative=[];
U=[1:size(dis,1)]';
[m,inx]=MaxMat(dis,U);
u=U(inx);
representative(end+1,1)=u;
U(inx,:)=[];
for i=1:size(U,1)
    a(U(i),1)=a(U(i),1)+ dis(U(i),u);
end
for i=1:size(representative,1)
    a(representative(i),1)=a(representative(i),1)+ dis(u,u);
end
while size(representative,1)<num
    [~,inx]= max(a(U),[],1);
    u=U(inx);
    representative(end+1,1)=u;
    U(inx,:)=[];
    for i=1:size(U,1)
        a(U(i),1)=a(U(i),1)+ dis(U(i),u);
    end
    for i=1:size(representative,1)-1
        a(representative(i),1)=a(representative(i),1)+ dis(representative(i),u);
    end
end
end

function prototype=shrink_tocenter(D,prototype,center_shrink)
centre(1,:)=mean(D,1);
centre=repmat(centre,size(prototype,1),1);
prototype=((1-center_shrink)*prototype+center_shrink*centre);
end

function [m,i]=MaxMat(d,U)
d=tril(d,-1);
d=d(:,U);
m=max(d);
[m,i]=max(m,[],2);
end