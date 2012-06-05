function clustering

fnames = dir('..\data\dataset UCR\All train\*');
tradata = cell(1,length(fnames));

for k=3:length(fnames)
    fname = fnames(k).name;
    files_name{k-2}=fname;
end
result=[];
for dataset_no=3:3 %length(files_name)
    file_name=['..\data\dataset UCR\All train\' files_name{dataset_no}];
    disp(file_name);
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
    [nor_traj,t_traj]=Import_Data_UCR(1,rows,2,data_len+1,file_name);
    
     details(dataset_no,:)=clustering_Hybrid_3Level(nor_traj,k,p);
 
end

