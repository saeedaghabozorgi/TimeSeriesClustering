function clustering

fnames = dir('..\data\dataset UCR\All train\*');

for k=3:length(fnames)
    fname = fnames(k).name;
    files_name{k-2}=fname;
end
disp('Reading data ..');
for dataset_no=3:3 %length(files_name)
    file_name=['..\data\dataset UCR\All train\' files_name{dataset_no}];
    disp([files_name(dataset_no)] );
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

% -------------------------
%   matlabpool open 8
%  parfor dataset_no=14:length(files_name)
%      files_name(dataset_no);
%       details(dataset_no,:)=clustering_Hybrid_3Level(ds{dataset_no},cluster_count(dataset_no),pp{dataset_no});
%  end
%  matlabpool close





fileID = fopen('result.txt','a');
for dataset_no=3:length(files_name)
    disp(['-------------------',files_name(dataset_no)] );
    details(dataset_no,:)=clustering_Hybrid_3Level(ds{dataset_no},cluster_count(dataset_no),pp{dataset_no});
    
   % [c,~]= do_kMediod_time (ds{dataset_no},cluster_count(dataset_no),'Euclid',0,'RAW');
%     [c,~]=do_kMediod_time(ds{dataset_no},cluster_count(dataset_no),'SAXminDis',0,'SAX','alphabet_size',8,'compression_ratio',4);
%     [SSEP,SSEC,RI,purity,BCubed,ConEntropy,f_measure,jacard,FM,quality]= do_Evaluate(pp{dataset_no},c,ds{dataset_no},[],[]);
%     details(dataset_no,:)=[SSEP,SSEC,RI,purity,BCubed,ConEntropy,f_measure,jacard,FM,quality];

    fprintf(fileID,'dataset_no: %d \n',dataset_no);
    dlmwrite('result.txt',details(dataset_no,:) ,'-append','delimiter', '\t','newline','pc');
end
fprintf(fileID,'---------------------------------------------------------------------------------------------------------\n');
fprintf(fileID,'\n');
fclose(fileID);
a=1;