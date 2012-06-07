function [details]=clustering_Hybrid_3Level(nor_traj, k,p)
%% ----------Level 1----------------------------------------------
disp('level 1');
c=[];
[c,itr]= do_kModes_time(nor_traj,k,'Euclid',0,'SAX','alphabet_size',8,'compression_ratio',6);
%------------------
% Evaluation
%  [SSEP,SSEC,RI,purity,BCubed,ConEntropy,f_measure,jacard,FM,quality]= do_Evaluate(p,c,nor_traj,[],[]);
%  details_l1(dataset_no,:)=[k,rows,data_len,SSEP,SSEC,RI,purity,BCubed,ConEntropy,f_measure,jacard,FM,quality];
 Plot_time_series_luminate(0,0,c,p,[],nor_traj,[],k,2,0.5,1);

%% ----------Level 2--CAST(raw)--------------------------------------------
disp('level 2');
clusterCount=max(c);
for i=1:clusterCount
    newData=find(c(:,1)==i);
    if length(newData)<5
        c(newData,2)=0;
    else
           temp_c= do_CAST_time (nor_traj(newData),'DTW',-1,'SAX','alphabet_size',8,'compression_ratio',2,'dtw_bound',0.8);
      %  temp_c= do_CAST_time (nor_traj(newData),'DTW',-1,'RAW','dtw_bound',.9);
        c(newData,2)=temp_c;
    end
end
%-------
% to map raw objects to new clusters
c(:,3)=c(:,1)*1000+c(:,2);
[x,y]=sort(c);
clsNum=1;
c(y(1,3),4)=clsNum;
for i=2:length(c)
    if ( x(i,3) ~= x(i-1,3))
        clsNum=clsNum+1;
        c(y(i,3),4)=clsNum;
    else
        c(y(i,3),4)=clsNum;
    end
end
l2_clusterCount=max(c(:,4));
%--evaluation---------------------------------------------------------
%qual_2lev=Calculate_Cluster_correct_ratio(c(:,4),p);
%     N_reduction_2lev=1-(l2_clusterCount/rows);
%     purity2=Calculate_Cluster_Purity(c(:,4),p,1);
%     level2_details(dataset_no,:)=[details_l1(dataset_no,7),qual_2lev,l2_clusterCount,purity2,N_reduction_2lev];
%    Plot_time_series(0,0,c(:,4),p,[],nor_traj,t_traj,l2_clusterCount,2,2);
 Plot_time_series_luminate(0,0,c(:,4),p,[],nor_traj,[],l2_clusterCount,2,0.5,2);
%  [SSEP,SSEC,RS,purity,BCubed,ConEntropy,fm]= do_Evaluate(p,cc,nor_traj,class_center,center);

%% --Level 3-------------------------------------------------
% to make the prototypes
disp('level 3');
center=[];
for i=1:l2_clusterCount;
    center{i}=centre_mediod(c(:,4),i,nor_traj,'RAW','DTW','dtw_bound',1);
    weight(i,1)=length(find(c(:,4)==i));
    %      center{i}=centre_mediod(c(:,4),i,nor_traj,'SAX','DTW','alphabet_size',8,'compression_ratio',2,'dtw_bound',1);
end
  Plot_time_series_luminate(0,0,c(:,4),p,center,nor_traj,[],l2_clusterCount,2,0.2,2);
%-----------
if l2_clusterCount>k
    %  k=5
     %  [c3,Z]=do_Hierarchical_time(center,k,'DTW','average',-1,'dtw_bound',0.1);
     %  [c3,Z]=do_Hierarchical_time(center,k,'DTW','average',-1,'rep','SAX','alphabet_size',8,'compression_ratio',4,'dtw_bound',0.1);
    %   [c3,Z]=do_Hierarchical_time(center,k,'DTW','complete',-1,'dtw_bound',1);
    %   [c3,Z]=do_Hierarchical_time(center,k,'DTW','single',-1,'dtw_bound',0.9);
    %   [c3,itr]= do_kMeans_time (center,k,'DTW',0,'RAW','dtw_bound',1);
       [c3,~]= do_kMediod_time (center,k,'DTW',0,'SAX','alphabet_size',8,'compression_ratio',4,'dtw_bound',0.3,'weight',weight);
   %  [c3,~]= do_kMediod_time (center,k,'DTW',0,'RAW','dtw_bound',1);
    for j=1:k
        l3_mems=find(c3==j);
        for i=1:length(l3_mems)
            sub_members=find(c(:,4)==l3_mems(i));
            c(sub_members,5)=j;
        end
    end
else
    c(:,5)=c(:,4);
end

%--evaluation-------------------------------------------------------------
%h= dendrogram(Z);
% Plot_time_series(0,0,c(:,5),p,[],nor_traj,t_traj,k,3,2);
%--------------------
%-- to show the centers
for j=1:length(center)
    sub_members=find(c(:,4)==j);
    xx=p(sub_members);
    pp(j,1)=mode(xx); %  to find more frequent in this cluster
end
 Plot_time_series_luminate(0,0,c3,pp,[],center,[],k,2,0.5,3);
%--------------------
  Plot_time_series_luminate(0,0,c(:,5),p,[],nor_traj,[],k,2,0.5,4);
clus_center={};
for i=1:k
    clus_center{i}=centre_mean(c(:,5),i,nor_traj);
end
[SSEP,SSEC,RI,purity,BCubed,ConEntropy,f_measure,jacard,FM,quality]= do_Evaluate(p,c(:,5),nor_traj,[],[]);

details=[SSEP,SSEC,RI,purity,BCubed,ConEntropy,f_measure,jacard,FM,quality];
end



function mmean=centre_mean(c,clusterNum,nor_traj)

members=find(c(:,1)==clusterNum);
n=length(nor_traj{1});
cluster_mem=zeros(size(members,1),n);
for j=1:size(members,1)
    if (members(j) > length(nor_traj))
        o=1;
    end
    cluster_mem(j,:)=nor_traj{members(j)};
end
mmean=mean(cluster_mem,1);
end


function medoid=centre_mediod(c,clusterNum,nor_traj,rep,dis_method,varargin)
t=find(c(:,1)==clusterNum);
if isempty (t)
    medoid=[];
elseif length(t)==1
    medoid=nor_traj{t(1)};
else
    
    switch rep
        case 'SAX'
            nor_traj_rep=represent_TS(nor_traj(t),rep,varargin{:});
        case 'RAW'
            nor_traj_rep=nor_traj(t);
    end
    
    %find distance of objects in cluster
    dis=Mtx_Distance(nor_traj_rep,nor_traj_rep,'same',dis_method,varargin{:});
    
    %find the SSE
    dis=dis+dis';
    dis=dis.^2;
    Error=sum(dis);
    [s,m]=min(Error);
    medoid=nor_traj{t(m)};
end
end

function XXX()
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
    counter=counter+1;
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

