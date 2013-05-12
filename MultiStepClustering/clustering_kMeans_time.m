function clustering_kMeans_time

k=15;
rows=500;
fprintf('Importing data .. \n\r');
%[nor_traj,t_traj]=Import_Data_SwedishLeaf(1,rows,'sequence',2,129);
%[nor_traj,t_traj]=Import_Data_SwedishLeaf(1,rows,'sequence',2,129);
file_name='data\SwedishLeaf_TRAIN.txt';
[nor_traj,t_traj]=Import_Data_UCR(1,rows,2,129,file_name);
%demo_sax(nor_traj{3});
%----------------------------------------
fprintf('Dimensionaly reduction of  data .. \n\r');

data_len      = 128;

global nseg;
nseg          = 64;

global alphabet_size;
alphabet_size = 6;


%compression_ratio: original_data_len / symbolic_len
compression_ratio= data_len/nseg;

 for i=1:length(nor_traj)
     SAX_nor_traj{i}= timeseries2symbol(nor_traj{i}, data_len, nseg, alphabet_size);
 end

%----------------------------------------
maxRow=length(nor_traj);
for i=1:k
    center{i}=nor_traj{i};
    center_SAX{i}=SAX_nor_traj{i};
end

temp=zeros(maxRow,1);   % initialize as zero vector
itr=1;
while 1,
    itr=itr+1;
    % dis=Mtx_SAX_APX_Distance(SAX_nor_traj,center_SAX,'');
     dis=Mtx_SAX_Min_Distance(SAX_nor_traj,center_SAX,'',compression_ratio);
    %   dis=Mtx_Euclid_Distance(nor_traj,center,'cell_not_same');
    
    [z,c]=min(dis,[],2);  % find group matrix g
    if (c==temp | itr==100),
        break;          % stop the iteration
    else
        temp=c;         % copy group matrix to temporary variable
    end
    for i=1:k
        center{i}=centre_mean(c,i,nor_traj);
        center_SAX{i}= timeseries2symbol(center{i}, data_len, nseg, alphabet_size);
    end
end

clusterCount=max(c);

for j=1:clusterCount
    cl_count_h(j,1)=size(find(c==j),1);
end


% Import the file

train_data = importdata(file_name);
p=train_data(:,1);
for i=1:k
    train_cluster_center{i}=centre_mean(p,i,nor_traj);
end


% TRAIN_class_labels = newData1(:,1);     % Pull out the class labels.
% TRAIN_cluster_number=length(unique(TRAIN_class_labels));
SSEP(1,1)=calculate_Cluster_SSE(p,nor_traj,train_cluster_center);
SSEC(1,1)=calculate_Cluster_SSE(c,nor_traj,center);
RS(1,1)= calculate_Cluster_Rand_statistic(c,p);
purity(1,1)= Calculate_Cluster_Purity(c,p,2);


% Plot_time_series(2,2,c,center,nor_traj,t_traj,clusterCount,1);
% Plot_time_series(2,2,p,train_cluster_center,nor_traj,t_traj,clusterCount,2);
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


