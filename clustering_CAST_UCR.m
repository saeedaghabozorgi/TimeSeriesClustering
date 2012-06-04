function clustering_CAST_UCR

fnames = dir('..\data\dataset UCR\All train\*');
tradata = cell(1,length(fnames));

for k=3:length(fnames)
    
    fname = fnames(k).name;
    files_name{k-2}=fname;
end
result=[];
for dataset_no=16:length(files_name)
    file_name=['..\data\dataset UCR\All train\' files_name{dataset_no}];
    disp(file_name);
    train_data = importdata(file_name);
    TRAIN_class_labels = train_data(:,1);     % Pull out the class labels.
    k=length(unique(TRAIN_class_labels));
    rows=size(train_data,1);
    data_len= size(train_data,2)-1;
    [nor_traj,t_traj]=Import_Data_UCR(1,rows,2,data_len+1,file_name);
    c=[];
%     c= do_CAST_time (nor_traj);
    c= do_CAST_time (nor_traj,'Euclid',-1,-1,-1);
    clusterCount=max(c)
  
    center=[];
    for i=1:clusterCount
        center{i}=centre_mean(c,i,nor_traj);
    end
    

%     Plot_time_series(9,6,c,center,nor_traj,t_traj,k,2);
    
    p=train_data(:,1);
    if(min(p)==0)
        p=p+1;
    end
    for i=1:k
        train_cluster_center{i}=centre_mean(p,i,nor_traj);
    end
    
    
    SSEP=calculate_Cluster_SSE(p,nor_traj,train_cluster_center);
    SSEC=calculate_Cluster_SSE(c(:,1),nor_traj,center);
    RS= calculate_Cluster_Rand_statistic(c(:,1),p);
    purity= Calculate_Cluster_Purity(c(:,1),p,1);
    
    
    details(dataset_no,1)=k;
    details(dataset_no,2)=rows;
    details(dataset_no,3)=data_len;
    details(dataset_no,4)=SSEP;
    details(dataset_no,5)=SSEC;
    details(dataset_no,6)=RS;
    details(dataset_no,7)=purity;

end
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

