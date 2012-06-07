function clustering

fnames = dir('..\data\dataset UCR\All train\*');

for k=3:length(fnames)
    fname = fnames(k).name;
    files_name{k-2}=fname;
end
for dataset_no=1:length(files_name)
    file_name=['..\data\dataset UCR\All train\' files_name{dataset_no}];

    train_data = importdata(file_name);
    TRAIN_class_labels = train_data(:,1);     % Pull out the class labels.
    
    p=train_data(:,1);
    if(min(p)==-1)
        p(p==-1)=2;
    end
    if(min(p)==0)
        p(p==0)=2;
    end
    
    k=length(unique(TRAIN_class_labels));
    rows=size(train_data,1);
    data_len= size(train_data,2)-1;
    [nor_traj,~]=Import_Data_UCR(1,rows,2,data_len+1,file_name);
    cluster_count(dataset_no)=k;
    pp{dataset_no}=p;
    ds{dataset_no}=nor_traj;
  %   details(dataset_no,:)=clustering_Hybrid_3Level(nor_traj,k,p);
end
%  matlabpool open 8
% parfor dataset_no=1:length(files_name)
%      disp(file_name);
%      details(dataset_no,:)=clustering_Hybrid_3Level(ds{dataset_no},cluster_count(dataset_no),pp{dataset_no});
% end
% matlabpool close


fileID = fopen('result.txt','w');
for dataset_no=1:length(files_name)
%      disp(files_name.);
    files_name(dataset_no)
    details(dataset_no,:)=clustering_Hybrid_3Level(ds{dataset_no},cluster_count(dataset_no),pp{dataset_no});
    fprintf(fileID,'%6.2f %12.8f\n',details(dataset_no,:));
    fprintf(fileID,'\r\n');
end
fclose(fileID);


a=1;