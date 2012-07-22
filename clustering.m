function clustering

fnames = dir('..\data\dataset UCR\All train\*');

for k=3:length(fnames)
    fname = fnames(k).name;
    files_name{k-2}=fname;
end
disp('Reading data ..');
for dataset_no=3:3%length(files_name)
    file_name=['..\data\dataset UCR\All train\' files_name{dataset_no}];
    disp([files_name(dataset_no)] );
    train_data = importdata(file_name);
    TRAIN_class_labels = train_data(:,1);     % Pull out the class labels.
    
    p=train_data(:,1);
    if(min(p)==-1)
        p(p==-1)=max(p)+1;
    end
    if(min(p)==0)
        p(p==0)=max(p)+1;
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

%%
fileID = fopen('result.txt','a');
for dataset_no=3:3%length(files_name)
    disp(['-------------------',files_name(dataset_no)] );
                parameter={'l1_dis_method','SAXAPX','l1_dtw_bound',1,'l1_rep','SAX','l1_alphabet_size',8,'l1_compression_ratio',6};
     parameter=[parameter,  'l2_dis_method','DTW','l2_dtw_bound',.07,'l2_rep','RAW','l2_alphabet_size',8,'l2_compression_ratio',6];
     parameter=[parameter,  'l3_dis_method','DTW','l3_dtw_bound',0.7,'l3_rep','SAX','l3_alphabet_size',8,'l3_compression_ratio',6,'l3_alg','k-medoid'];
     details(dataset_no,:)=clustering_Hybrid_3Level(ds{dataset_no},cluster_count(dataset_no),pp{dataset_no},parameter{:});
%    Plot_time_series(2,1,pp{dataset_no},pp{dataset_no},[],ds{dataset_no},[],cluster_count(dataset_no),3,0);
    
%     dis=Mtx_Distance(ds{dataset_no},ds{dataset_no},'same','dis_method','Euclid','rep','RAW');
%     x=squareform(dis);
%     figure(1);
%     h = normplot(x);
%     figure(2);
%     hist(x)
    %     [c,~]= do_kMediod_time (ds{dataset_no},cluster_count(dataset_no),0,'dis_method','Euclid','rep','RAW');
    %     Plot_time_series_luminate(0,0,c,pp{dataset_no},[],ds{dataset_no},[],cluster_count(dataset_no),2,0.5,1);
    
    %     [c,~]= do_kMediod_time (ds{dataset_no},cluster_count(dataset_no),0,'dis_method','SAXminDis','rep','SAX','alphabet_size',8,'compression_ratio',4);
     %     [c,Z]=do_Hierarchical_time(ds{dataset_no},cluster_count(dataset_no),'average',-1,'dis_method','SAXminDis','rep','SAX','alphabet_size',8,'compression_ratio',4,'dtw_bound',0.8);

       %   [SSEP,SSEC,RI,purity,BCubed,ConEntropy,f_measure,jacard,FM,quality]= do_Evaluate(pp{dataset_no},c,ds{dataset_no},[],[]);
%          details(dataset_no,:)=[SSEP,SSEC,RI,purity,BCubed,ConEntropy,f_measure,jacard,FM,quality];
%          h= dendrogram(Z);
%       Plot_time_series_luminate(0,0,c, pp{dataset_no},[],ds{dataset_no},[],cluster_count(dataset_no),2,0.5,4);
%     fprintf(fileID,'dataset_no: %d \n',dataset_no);
%     dlmwrite('result.txt',details(dataset_no,:) ,'-append','delimiter', '\t','newline','pc');
end
fprintf(fileID,'---------------------------------------------------------------------------------------------------------\n');
fprintf(fileID,'\n');
fclose(fileID);
a=1;