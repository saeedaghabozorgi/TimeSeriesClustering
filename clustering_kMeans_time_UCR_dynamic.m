function clustering_kMeans_time_UCR_dynamic
% to test the trend of centers in each point of time
fnames = dir('..\data\dataset UCR\All train\*');
tradata = cell(1,length(fnames));

for k=3:length(fnames)
    fname = fnames(k).name;
    files_name{k-2}=fname;
end


for dataset_no=3:length(files_name)
    file_name=['..\data\dataset UCR\All train\' files_name{dataset_no}];
    disp(file_name);
    train_data = importdata(file_name);
    TRAIN_class_labels = train_data(:,1);     % Pull out the class labels.
    k=length(unique(TRAIN_class_labels));
    rows=size(train_data,1);
    data_len= size(train_data,2)-1;
    k=5;
    data_c=train_data(:,2:end);
    % data_c=data_normalize(data_c,'range');
    
    hold off;
    mycenter=zeros(k,1);
    for i=1:data_len+1
        data=data_c(:,1:i);
        [Classes,Centres,FinalDistance]=do_kMeans(data,k);
        mycenter=Centres;
        figure(1);
        plot(mycenter');
        pause;
    end
    
    [B,IX] = sort(mycenter);
    figure(2);
    plot(B');
    [nor_traj,t_traj]=Import_Data_UCR(1,rows,2,data_len+1,file_name);
    
    c=[];
    [c,type]= do_kMeans_time (nor_traj,k,'Euclid',0)  ;
    
    
    
    center=[];
    for i=1:k
        center{i}=centre_mean(c,i,nor_traj);
    end
    %--------------------------
    dis=Mtx_Euclid_Distance(nor_traj,nor_traj,'same');
    dis=dis+dis';
    dis=zscore(dis);
    [RV,C,I,RI]=VAT(dis);
    
    imshow(RV);
    
    %  Plot_time_series(0,0,c,center,nor_traj,t_traj,k,1);
    
    %--------------------------
    
    p=train_data(:,1);
    if(min(p)==0)
        p=p+1;
    end
    for i=1:k
        train_cluster_center{i}=centre_mean(p,i,nor_traj);
    end
    
    %   Plot_time_series(0,0,p,train_cluster_center,nor_traj,t_traj,k,2);
    
    SSEP=Calculate_Cluster_SSE(p,nor_traj,train_cluster_center);
    SSEC=Calculate_Cluster_SSE(c,nor_traj,center);
    RS= Calculate_Cluster_Rand_statistic(c,p);
    purity= Calculate_Cluster_Purity(c,p,2);
    
    
    
    
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
t=find(c(:,1)==clusterNum);
n=length(nor_traj{1});
cluster_mem=zeros(1,n);
for j=1:size(t,1)
    cluster_mem(j,:)=nor_traj{t(j)};
end
mmean=mean(cluster_mem,1);
end


