function clustering_ds_UCR

fnames = dir('..\data\dataset UCR\All train\*');

for k=3:length(fnames)
    fname = fnames(k).name;
    files_name{k-2}=fname;
end
disp('Reading data ..');
for dataset_no=1:length(files_name)
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
    [nor_traj_raw,~]=Import_Data_UCR(1,rows,2,data_len+1,file_name);
    cluster_count(dataset_no)=k;
    pp{dataset_no}=p;
    ds{dataset_no}=nor_traj_raw;
    %   details(dataset_no,:)=clustering_Hybrid_3Level(nor_traj,k,p);
end

%%
fileID = fopen('result.txt','a');
for dataset_no=3:length(files_name)
    disp(['-------------------',files_name(dataset_no)] );
 %   details(dataset_no,:)=evaluate_distance( ds{dataset_no});
 %   details(dataset_no,:)=clustering_Hybrid_3Level(ds{dataset_no},cluster_count(dataset_no),pp{dataset_no});
calculate_DTW_ground_truth(ds{dataset_no},dataset_no);
    
end
fprintf(fileID,'---------------------------------------------------------------------------------------------------------\n');
fprintf(fileID,'\n');
fclose(fileID);

end

function calculate_DTW_ground_truth(nor_traj_raw,dataset_no)
dist=[];
dist=Mtx_Distance(nor_traj_raw,nor_traj_raw,'same','Org', 'dis_method','DTW','dtw_bound',1,'rep','RAW');
filename=['result_dis_',num2str(dataset_no),'.mat'];
save(filename, 'dist')
% dlmwrite(filename,dist1 ,'-append','delimiter','\t','newline','pc');
end

function do_others()
%plot_histogram();
    %     [c,~]= do_kMediod_time (ds{dataset_no},cluster_count(dataset_no),0,'dis_method','Euclid','rep','RAW');
    %     Plot_time_series_luminate(0,0,c,pp{dataset_no},[],ds{dataset_no},[],cluster_count(dataset_no),2,0.5,1);
    
    %     [c,~]= do_kMediod_time (ds{dataset_no},cluster_count(dataset_no),0,'dis_method','SAXminDis','rep','SAX','alphabet_size',8,'compression_ratio',4);
    %     [c,Z]=do_Hierarchical_time(ds{dataset_no},cluster_count(dataset_no),'average',-1,'dis_method','SAXminDis','rep','SAX','alphabet_size',8,'compression_ratio',4,'dtw_bound',0.8);
    
    %   [SSEP,SSEC,RI,purity,BCubed,ConEntropy,f_measure,jacard,FM,quality]= do_Evaluate(pp{dataset_no},c,ds{dataset_no},[],[]);
    %          details(dataset_no,:)=[SSEP,SSEC,RI,purity,BCubed,ConEntropy,f_measure,jacard,FM,quality];
    %          h= dendrogram(Z);
    %       Plot_time_series_luminate(0,0,c, pp{dataset_no},[],ds{dataset_no},[],cluster_count(dataset_no),2,0.5,4);
    %     fprintf(fileID,'dataset_no: %d \n',dataset_no);
end


function plot_histogram()
    %    Plot_time_series(2,1,pp{dataset_no},pp{dataset_no},[],ds{dataset_no},[],cluster_count(dataset_no),3,0);
    
    %     dis=Mtx_Distance(ds{dataset_no},ds{dataset_no},'same','dis_method','Euclid','rep','RAW');
    %     x=squareform(dis);
    %     figure(1);
    %     h = normplot(x);
    %     figure(2);
    %     hist(x)
end

function result=evaluate_distance(nor_traj_raw)
%-------------------------------------------
dist1=Mtx_Distance(nor_traj_raw,nor_traj_raw,'same','Org', 'dis_method','Euclid','dtw_bound',1,'rep','RAW','alphabet_size',4,'compression_ratio',8);
% representation
nor_traj=represent_TS(nor_traj_raw,'SAX','alphabet_size',4,'compression_ratio',8);


dist2=Mtx_Distance(nor_traj,nor_traj,'same','Org', 'dis_method','SAXDIST','dtw_bound',1,'rep','SAX','alphabet_size',4,'compression_ratio',8);
Avg_Diff1=sum(abs(squareform(dist1-dist2)))/length(squareform(dist1));
NTighetness1=((dist1-abs(dist1-dist2))./dist1);
NTighetness1(1:length(NTighetness1)+1:length(NTighetness1)*length(NTighetness1))=0;
Avg_NTighetness1= mean(mean(NTighetness1));

dist3=Mtx_Distance(nor_traj,nor_traj,'same','Org', 'dis_method','SAXminDis','dtw_bound',1,'rep','SAX','alphabet_size',4,'compression_ratio',8);
Avg_Diff2=sum(abs(squareform(dist1-dist3)))/length(squareform(dist3));
NTighetness2=((dist1-abs(dist1-dist3))./dist1);
NTighetness2(1:length(NTighetness2)+1:length(NTighetness2)*length(NTighetness2))=0;
Avg_NTighetness2= mean(mean(NTighetness2));

dist4=Mtx_Distance(nor_traj,nor_traj,'same','Org', 'dis_method','Euclid','dtw_bound',1,'rep','SAX','alphabet_size',4,'compression_ratio',8);
Avg_Diff3=sum(abs(squareform(dist1-dist4)))/length(squareform(dist4));
NTighetness3=((dist1-abs(dist1-dist4))./dist1);
NTighetness3(1:length(NTighetness3)+1:length(NTighetness3)*length(NTighetness3))=0;
Avg_NTighetness3= mean(mean(NTighetness3));

result=[Avg_Diff1,Avg_Diff2,Avg_Diff3,Avg_NTighetness1,Avg_NTighetness2,Avg_NTighetness3]

dlmwrite('result.txt',result ,'-append','delimiter','\t','newline','pc');
%-------------------------------------------------------------
end

