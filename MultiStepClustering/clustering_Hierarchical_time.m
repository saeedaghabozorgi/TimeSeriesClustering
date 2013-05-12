function clustering_Hierarchical_time()

fprintf('Importing data .. \n\r');
[nor_traj,t_traj]=ImportBankData(1,200,'sequence',1,256);

%----------------------------------------

%dis=Mtx_SAX_APX_Distance(SAX_nor_traj,SAX_nor_traj,'same');
%dis=Mtx_SAX_Min_Distance(SAX_nor_traj,SAX_nor_traj,'same');
dis=Mtx_Euclid_Distance(nor_traj,nor_traj,'same');
dis= squareform(dis+dis');

Zh = linkage(dis,'average');
c = cluster(Zh,'maxclust',10);
clusterCount=max(c)

for j=1:clusterCount
    newOffMeanCenter{j}=centre_mean(c,j,nor_traj);
    cl_count_h(i,1)=size(find(c==i),1);
end

disp(cl_count_h);
Plot_time_series(4,3,c,newOffMeanCenter,nor_traj,t_traj,clusterCount);


end

function mmean=centre_mean(c,clusterNum,nor_traj)
t=find(c(:,1)==clusterNum);
n=length(nor_traj{1});
cluster_mem=zeros(1,n);
for j=1:size(t,1)
    cluster_mem(j,:)=nor_traj{t(j)};
    
end
mmean=mean(cluster_mem,1);
end





