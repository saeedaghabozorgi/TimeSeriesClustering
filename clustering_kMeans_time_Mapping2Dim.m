function clustering_kMeans_time_Mapping2Dim
% to test the trend of centers in each point of time
fnames = dir('..\data\dataset UCR\All train\*');
tradata = cell(1,length(fnames));

for k=3:length(fnames)
    fname = fnames(k).name;
    files_name{k-2}=fname;
end


for dataset_no=1:length(files_name)
    file_name=['..\data\dataset UCR\All train\' files_name{dataset_no}];
    disp(file_name);
    train_data = importdata(file_name);
    TRAIN_class_labels = train_data(:,1);     % Pull out the class labels.
    k=length(unique(TRAIN_class_labels));
    rows=size(train_data,1);
    data_len= size(train_data,2)-1;
    k=10
    data_c=train_data(:,2:end);
    data_c=data_normalize(data_c,'range');
    
    
    
    [Classes,Centres,FinalDistance]=do_kMeans(data_c,k);
    
    mappedX= mds(data_c,2);
    hold off;
    cc=hsv(k);
    for i=1:k
        scatter(mappedX(Classes==i,1),mappedX(Classes==i,2),12,cc(i,:))
        hold on
    end
    
    
    scatter(mappedX(:,1),mappedX(:,2),5,'r');
    
end

end

