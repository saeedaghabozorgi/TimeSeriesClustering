function clustering_ds_UCR
foldpath='..\dataset UCR';
nameFolds = dir(foldpath);
for k=3:length(nameFolds)
    fold_name = nameFolds(k).name;
    folderes_name{k-2}=fold_name;
    files_name{k-2}=[fold_name,'_TEST'];
end

disp('Reading data ..');
for dataset_no=1:length(files_name)
    file_name=[foldpath,'\',folderes_name{dataset_no},'\',files_name{dataset_no}];
    disp([files_name(dataset_no)] );
    train_data = importdata(file_name);
    
    %----labels from data------
    TRAIN_class_labels = train_data(:,1);
    if(min(TRAIN_class_labels)==-1)
        TRAIN_class_labels(TRAIN_class_labels==-1)=max(TRAIN_class_labels)+1;
    end
    if(min(TRAIN_class_labels)==0)
        TRAIN_class_labels(TRAIN_class_labels==0)=max(TRAIN_class_labels)+1;
    end
    
    %----labels from DTW_k-medoids------
    %     label_filename=['..\dataset UCR\prepared_k-medoids_label\ds_UCR_',num2str(dataset_no),'_DTW_label.mat'];
    %     load(label_filename,'label');
    %     TRAIN_class_labels=label;
    
    
    
    p=TRAIN_class_labels;
    
    k=length(unique(TRAIN_class_labels));
    rows=size(train_data,1);
    data_len= size(train_data,2)-1;
    [nor_traj_raw,~]=Import_Data_UCR(1,rows,2,data_len+1,file_name);
    cluster_count(dataset_no)=k;
    pp{dataset_no}=p;
    ds{dataset_no}=nor_traj_raw;
    
    
    dist_mtx_DTW_file=[foldpath,'\',folderes_name{dataset_no},'\',files_name{dataset_no},'_dismat_DTW.mat'];
    load(dist_mtx_DTW_file,'dismat');
    dismat=dismat+dismat';
    dist_mtx_DTW{dataset_no}=dismat;
end

%%

for dataset_no=1:length(files_name)
    disp(['-------------------',files_name(dataset_no)] );
    %   details(dataset_no,:)=evaluate_distance( ds{dataset_no});
    %     for i=1: 20
    %         def(i,:)= clustering_Hybrid_3Level(ds{dataset_no},cluster_count(dataset_no),pp{dataset_no},dist_mtx_DTW{dataset_no});
    %     end
    %   details{dataset_no}= def;
  
    
      evaluate_distance(ds{dataset_no})


    
%details{dataset_no}=clustering_anytime2(ds{dataset_no},cluster_count(dataset_no),pp{dataset_no},dist_mtx_DTW{dataset_no});
    %   details(dataset_no,:)=clustering_Hybrid_2Level_graphbased(ds{dataset_no},cluster_count(dataset_no),pp{dataset_no},dist_mtx_DTW{dataset_no});
   %    details(dataset_no,:)=clustering_Hybrid_3Level(ds{dataset_no},cluster_count(dataset_no),pp{dataset_no},dist_mtx_DTW{dataset_no});
    %     [det,hurestic_param]=clustering_Hybrid_3Level_heuristic(ds{dataset_no},cluster_count(dataset_no),pp{dataset_no},dist_mtx_DTW{dataset_no});
    %     details(dataset_no,:)=det;
    %     param{dataset_no}=hurestic_param;
    
    
    %   calculate_DTW_matrix_paralele(ds{dataset_no},[foldpath,'\',folderes_name{dataset_no},'\',files_name{dataset_no},'_dismat_DTW.mat']);
    %   calculate_DTW_ground_truth(ds{dataset_no},dataset_no);
    %details(dataset_no,:)=clustering_Hybrid_3Level_anytime3(ds{dataset_no},cluster_count(dataset_no),pp{dataset_no});
    %  plot_histogram(dist_mtx_DTW{dataset_no},ds{dataset_no});
 %   save('details.mat','details');

end
end

function test_2ltc()
cluster_count=[20	2	10	2	4	9	11	2	4	4	4	7	12	8	2	4	5	7	6]';

details(dataset_no,:)=clustering_Hybrid_2Level(ds{dataset_no},cluster_count(dataset_no),pp{dataset_no},dist_mtx_DTW{dataset_no});
end

function do_hier_shape()

   x=[ 2  19 23 28 16 20 15 22 17]';
ds{dataset_no}=ds{dataset_no}(x);
pp{dataset_no}=pp{dataset_no}(x);
dist_mtx{dataset_no}=Mtx_Distance(ds{dataset_no},ds{dataset_no},'same','Org','dis_method','Euclid');
%dist_mtx{dataset_no}=dist_mtx_DTW{dataset_no}(x,x);
[c,Z]=do_Hierarchical_time(ds{dataset_no},cluster_count(dataset_no),'single',-1,dist_mtx{dataset_no});
[SSEP,SSEC,RI,ARI,purity,ConEntropy,f_measure,jacard,FM,NMI,CSM,quality]= do_Evaluate( pp{dataset_no},c,ds{dataset_no},[],[]);
figure(1);
details=[SSEP,SSEC,RI,ARI,purity,ConEntropy,f_measure,jacard,FM,NMI,CSM,quality];
[H,T,perm] =  dendrogram(Z,'orientation','left');

plot_time_sereis(perm,pp{dataset_no},ds{dataset_no})
Plot_time_series_luminate(0,0,c,pp{dataset_no},[],ds{dataset_no},[],cluster_count(dataset_no),2,0.5,1);


end

function calculate_DTW_matrix(nor_traj_raw,dataset_no)
dist=[];
dist=Mtx_Distance(nor_traj_raw,nor_traj_raw,'same','Org', 'dis_method','DTW','dtw_bound',1,'rep','RAW');
filename=['result_dis_',num2str(dataset_no),'.mat'];
save(filename, 'dist')
% dlmwrite(filename,dist1 ,'-append','delimiter','\t','newline','pc');
end

function find_label_k_medoid_DTW(dist_mtx_DTW,k,dataset_no)
dist_mtx_DTW=squareform(dist_mtx_DTW);
[label,~]= do_kMedoids_keogh(k,dist_mtx_DTW);
filename=['..\dataset UCR\prepared_k-medoids_label\ds_UCR_',num2str(dataset_no),'_DTW_label.mat'];
save(filename, 'label')
end

function calculate_DTW_matrix_paralele(nor_traj,path)
matlabpool open 2
poolSize = matlabpool('size');
if poolSize == 0
    error('parallel:demo:poolClosed', ...
        'This demo needs an open MATLAB pool to run.');
end
fprintf('This demo is running on %d MATLABPOOL workers.\n', ...
    matlabpool('size'));

x=nor_traj;
data_n = length(x);
dismat=zeros(data_n,data_n);
parfor  i = 1:data_n
    disp(num2str(i));
    dismat(i,:)=rowcalc(i,data_n,x);
end
save(path, 'dismat')
matlabpool close
end

function dismatrix=rowcalc(row,data_n,x)
dismatrix=zeros(1,data_n);
a=x{row};
for j = row:data_n,
    if row ~= j
        b=x{j};
        dis=dis_dtw3(a,b,size(a,2));
        dismatrix(1, j) =dis;
    end
end
end

function do_others()
%plot_histogram();
%     [c,~]= do_kMediod_time (ds{dataset_no},cluster_count(dataset_no),0,'dis_method','Euclid','rep','RAW');
%     Plot_time_series_luminate(0,0,c,pp{dataset_no},[],ds{dataset_no},[],cluster_count(dataset_no),2,0.5,1);

%     [c,~]= do_kMediod_time (ds{dataset_no},cluster_count(dataset_no),0,'dis_method','SAXminDis','rep','SAX','alphabet_size',8,'compression_ratio',4);
[c,Z]=do_Hierarchical_time(ds{dataset_no},cluster_count(dataset_no),'average',-1,'dis_method','SAXminDis','rep','SAX','alphabet_size',8,'compression_ratio',4,'dtw_bound',0.8);

%   [SSEP,SSEC,RI,purity,BCubed,ConEntropy,f_measure,jacard,FM,quality]= do_Evaluate(pp{dataset_no},c,ds{dataset_no},[],[]);
%          details(dataset_no,:)=[SSEP,SSEC,RI,purity,BCubed,ConEntropy,f_measure,jacard,FM,quality];
%          h= dendrogram(Z);
%       Plot_time_series_luminate(0,0,c, pp{dataset_no},[],ds{dataset_no},[],cluster_count(dataset_no),2,0.5,4);
%     fprintf(fileID,'dataset_no: %d \n',dataset_no);
end


function plot_histogram(dis,data)
%    Plot_time_series(2,1,pp{dataset_no},pp{dataset_no},[],ds{dataset_no},[],cluster_count(dataset_no),3,0);
if isempty(dis)
    dis=Mtx_Distance(data,data,'same','dis_method','Euclid','rep','RAW');
end
Nor = dis - min( dis(:) );
if max( Nor(:) ) ~= 0
    dis = Nor / max( Nor(:) );
else
    dis=Nor;
end
dis=squareform(dis);
figure(1);
h = normplot(dis);
figure(2);
hist(dis)

%------------------------------------------
delete('test2.txt');
 nor_traj=represent_TS(ds{dataset_no},'SAX','alphabet_size',6,'compression_ratio',4);
        A = importdata('test2.txt');
        nor_traj_raw=ds{dataset_no};
        figure(1);
         hold off;
%         plot (nor_traj_raw{1});
%         hold on;
%         plot (nor_traj_raw{2},'r');
%         plot (nor_traj_raw{3},'g');
%         print(gcf,'-dpng',['TS',int2str(dataset_no),'.png']);
       
         hist(A,50)
         print(gcf,'-dpng',['dis_all',int2str(dataset_no),'.png']);
    %    h = normplot(A) 


end

function result=evaluate_distance(nor_traj_raw)
%-------------------------------------------
dist1=Mtx_Distance(nor_traj_raw,nor_traj_raw,'same','Org', 'dis_method','Euclid','dtw_bound',1,'rep','RAW','alphabet_size',4,'compression_ratio',8);
dist1=squareform(dist1);
% representation
nor_traj=represent_TS(nor_traj_raw,'SAX','alphabet_size',4,'compression_ratio',8);


dist2=Mtx_Distance(nor_traj,nor_traj,'same','Org', 'dis_method','SAXDIST','dtw_bound',1,'rep','SAX','alphabet_size',4,'compression_ratio',8);
dist2=squareform(dist2);
Avg_Diff1=sum(abs(dist1-dist2))/length(dist1);
NTighetness1=((dist1-abs(dist1-dist2))./dist1);
NTighetness1(find(isnan(NTighetness1)))=1;
Avg_NTighetness1= (mean(NTighetness1));

dist3=Mtx_Distance(nor_traj,nor_traj,'same','Org', 'dis_method','SAXminDis','dtw_bound',1,'rep','SAX','alphabet_size',4,'compression_ratio',8);
dist3=squareform(dist3);
Avg_Diff2=sum(abs((dist1-dist3)))/length((dist3));
NTighetness2=((dist1-abs(dist1-dist3))./dist1);
NTighetness2(find(isnan(NTighetness2)))=1;
Avg_NTighetness2= (mean(NTighetness2));

dist4=Mtx_Distance(nor_traj,nor_traj,'same','Org', 'dis_method','Euclid','dtw_bound',1,'rep','SAX','alphabet_size',4,'compression_ratio',8);
dist4=squareform(dist4);
Avg_Diff3=sum(abs((dist1-dist4)))/length((dist4));
NTighetness3=((dist1-abs(dist1-dist4))./dist1);
NTighetness3(find(isnan(NTighetness3)))=1;
Avg_NTighetness3= mean((NTighetness3));


dist4=Mtx_Distance(nor_traj,nor_traj,'same','Org', 'dis_method','SAX_MEANDIST','dtw_bound',1,'rep','SAX','alphabet_size',4,'compression_ratio',8);
dist4=squareform(dist4);
Avg_Diff4=sum(abs(dist1-dist4))/length(dist1);
NTighetness4=((dist1-abs(dist1-dist4))./dist1);
NTighetness4(find(isnan(NTighetness4)))=1;
Avg_NTighetness4= (mean(NTighetness4));

result=[Avg_Diff1,Avg_Diff2,Avg_Diff3,Avg_Diff4,Avg_NTighetness1,Avg_NTighetness2,Avg_NTighetness3,Avg_NTighetness4]

dlmwrite('result.txt',result ,'-append','delimiter','\t','newline','pc');
%-------------------------------------------------------------
end

