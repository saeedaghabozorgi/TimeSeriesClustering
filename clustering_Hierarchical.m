function clustering_Hierarchical()
%demo
% D=[0 0;1 1; 2 1; 0 4; 2 5; 3 5; 3 4; 4 5; 4 4];
% dis=dis_euclidean_matrix(D,D);
% dis=squareform(dis);
% 
% cluster = do_Hierarchical(dis, 'single', 3)
%---


fprintf('Importing data .. \n\r');
[nor_traj,t_traj]=ImportBankData(1,200,'sequence',1,256);
%----------------------------------------
fprintf('Dimensionaly reduction of  data .. \n\r');

data_len      = 256;
global nseg;
nseg          = 8;

global alphabet_size;
alphabet_size = 32;

for i=1:length(nor_traj)
    SAX_nor_traj{i}= timeseries2symbol(nor_traj{i}, data_len, nseg, alphabet_size);
end

%----------------------------------------
cc=hsv(100);


fprintf('Hierarchal clustering ... \n\r');
fprintf('Hierarchal clustering ... Distance Calculation ... \n\r');
dis=Mtx_SAX_APX_Distance(SAX_nor_traj,SAX_nor_traj,'same');
%dis=Mtx_SAX_Min_Distance(SAX_nor_traj,SAX_nor_traj,'same');
%dis=Mtx_Euclid_Distance(nor_traj,nor_traj,'same');
dis= squareform(dis+dis');


Zh = linkage(dis,'average');
c = cluster(Zh,'maxclust',10);
clusterCount=10;

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





