function clustering_kMeans_time_Raw
ab=[];

    k=15;
    fprintf('Importing data .. \n\r');
    nor_traj=[];
    rows=500;
    %[nor_traj,t_traj,TRAIN_class_labels]=Import_Data_SwedishLeaf(1,rows,'sequence',2,129);
    [nor_traj,t_traj]=Import_Data_Bank(1,300,'sequence',1,360);
    
    %----------------------------------------
    data_len  = 360;
    %----------------------------------------
    
    maxRow=length(nor_traj);
    for i=1:k
        center{i}=nor_traj{i};
    end
    
    temp=zeros(maxRow,1);   % initialize as zero vector
    itr=1;
    while 1,
        itr=itr+1;
        % dis=Mtx_SAX_APX_Distance(SAX_nor_traj,center_SAX,'');
        %dis=Mtx_SAX_Min_Distance(nor_traj,center,'');
        dis=Mtx_Euclid_Distance(nor_traj,center,'cell_not_same');
        
        [z,c]=min(dis,[],2);  % find group matrix g
        if (c==temp | itr==100),
            break;          % stop the iteration
        else
            temp=c;         % copy group matrix to temporary variable
        end
        for i=1:k
            center{i}=centre_mean(c,i,nor_traj);
        end
    end
    clusterCount=max(c);
Plot_time_series(0,0,c,center,nor_traj,t_traj,clusterCount,1,1);
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


