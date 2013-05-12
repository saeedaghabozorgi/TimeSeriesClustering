function clustering_time_UCR
%this function is used for SAX and RAW
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
    
    compression_ratio=2;
    alphabet_size = 8;
    
    
    data_len=(floor(data_len/compression_ratio))*compression_ratio;
    [nor_traj,t_traj]=Import_Data_UCR(1,rows,2,data_len+1,file_name);
    
    nseg = data_len/compression_ratio;
    SAX_nor_traj=[];
    for i=1:length(nor_traj)
        SAX_nor_traj{i}= rep_SAX(nor_traj{i}, data_len, nseg, alphabet_size);
    end
    xxx=1000;
    for rr=1:xxx
        rr
        [c,itr]= do_kMediod_time(SAX_nor_traj,k,'SAXminDis',1,alphabet_size,compression_ratio,'SAX');
        %  [c]=do_Hierarchical_time(SAX_nor_traj,k,'SAXminDis','average',-1,alphabet_size,compression_ratio);
        %  [c,itr]= do_kMeans_time (nor_traj,k,'Euclid',0,-1,-1,'RAW')  ;
        ccc(:,rr)=c;
    end
    for i=1:xxx
        for j=1:xxx
            sim(i,j)=Calculate_Cluster_Rand_statistic(ccc(:,i),ccc(:,j));
        end
    end
    
    ct=[1:1:xxx]
    clust_num=0;
    while ~isempty(ct)
        item=ct(1,1);
        clust_num=clust_num+1
        clus{clust_num}=item;
        ct(1)=[];
        i=1;
        while i<=length(ct)
            if sim(item,ct(1,i))==1
                 clus{clust_num}=[clus{clust_num};ct(1,i)]
                 ct(i)=[];
                
            else
            i=i+1
            end
        end    
    end
    dis=sim*(-1)+1;
dis=squareform(dis)
    Zh = linkage(dis,'average');
    H = dendrogram(Zh)
    ct = cluster(Zh,'cutoff',0.01);
    Plot_time_series_luminate(0,0,c,p,[],nor_traj,t_traj,k,1,2);
    
    % to sort all the class labels
    if(length(unique(c))<k)
        [f,g]=sort(unique(c));
        for i=1:max(g)
            c(c==f(i))=g(i);
        end
    end
    %---------------------------------
    for i=1:k
        clus_center{i}=centre_mean(c,i,nor_traj);
    end
    %Plot_time_series(0,0,c,cluster_center,nor_traj,t_traj,k,2,1)
    
    for i=1:k
        class_center{i}=centre_mean(p,i,nor_traj);
    end
    % Plot_time_series(0,0,p,class_center,nor_traj,t_traj,k,2,1)
    %--Original K-mediod------------------------------
    %      p=[];
    %      [p,itr]= do_kMediod_time (nor_traj,k,'Euclid',0,-1,-1,'RAW')  ;
    
    [SSEP,SSEC,RI,purity,BCubed,ConEntropy,f_measure,jacard,FM,quality]= do_Evaluate(p,c,nor_traj,-1,-1);
    details(dataset_no,:)=[k,rows,data_len,SSEP,SSEC,RI,purity,BCubed,ConEntropy,f_measure,jacard,FM,quality];
end

end


%%
function medoid=medoid(c,clusterNum,nor_traj)
t=find(c(:,1)==clusterNum);
if isempty(t)
    mmean=[];
else
    dis=Mtx_SAX_Min_Distance(nor_traj(t),nor_traj(t),'same');
    sum_dis=sum(dis);
    [s,m]=min(sum_dis);
    mmean=nor_traj{t(m)};
end
medoid=mmean;
end

function mmean=centre_mean(c,clusterNum,nor_traj)
t=find(c(:,1)==clusterNum);
if isempty(t)
    p = randperm(length(nor_traj));      % random initialization
    mmean= nor_traj{p(1)};
else
    n=length(nor_traj{1});
    cluster_mem=zeros(1,n);
    for j=1:size(t,1)
        cluster_mem(j,:)=nor_traj{t(j)};
    end
    mmean=mean(cluster_mem,1);
end
end
% function mmean=centre_mean(c,clusterNum,nor_traj)
% t=find(c(:,1)==clusterNum);
% n=length(nor_traj{1});
% cluster_mem=zeros(1,n);
% for j=1:size(t,1)
%     cluster_mem(j,:)=nor_traj{t(j)};
% end
% mmean=mean(cluster_mem,1);
% end