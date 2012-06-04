function Evaluating_approx_Hierarichal_time

fnames = dir('..\data\dataset UCR\All train\*');
tradata = cell(1,length(fnames));

for k=3:length(fnames)
    
    fname = fnames(k).name;
    files_name{k-2}=fname;
end


for dataset_no=2:length(files_name)
    file_name=['..\data\dataset UCR\All train\' files_name{dataset_no}];
    disp(file_name);
    train_data = importdata(file_name);
    TRAIN_class_labels = train_data(:,1);     % Pull out the class labels.
    k=length(unique(TRAIN_class_labels));
    rows=size(train_data,1);
    data_len= size(train_data,2)-1;
    compression_ratio=2;
    alphabet_size = 6;
    
    data_len=(floor(data_len/compression_ratio))*compression_ratio;
    [nor_traj,t_traj]=Import_Data_UCR(1,rows,2,data_len+1,file_name);
    
    nseg = data_len/compression_ratio;
    SAX_nor_traj=[];
    for i=1:length(nor_traj)
        SAX_nor_traj{i}= timeseries2symbol(nor_traj{i}, data_len, nseg, alphabet_size);
    end
    
    %----------Level 1----------------------------------------------
    
         dis=Mtx_Euclid_Distance(nor_traj,nor_traj,'same');
        dis= squareform(dis+dis');
        Zh = linkage(dis,'average');
    
        Ap_dis=Mtx_SAX_Max_Distance(SAX_nor_traj,'SAXminDis','same',alphabet_size,compression_ratio);
        Ap_dis= squareform(Ap_dis+Ap_dis');
        Ap_Zh = linkage(Ap_dis,'average');
 
    

        
    for k1=2:k
        k1
        c=[];
        cutoff=0.4+0.025 *k1
       % cp = cluster(Zh,'maxclust',k1);
        cp =cluster(Zh, 'distance',cutoff);
        cp_count=max(cp)
        cpcop = cophenet(Zh,dis);
        
        
        c = cluster(Ap_Zh,'maxclust',cp_count);
%          c =cluster(Ap_Zh,'cutoff',cutoff);
%         c_count=max(c);

        cop = cophenet(Ap_Zh,Ap_dis);
         
        
        RS= calculate_Cluster_Rand_statistic(c,cp);
        purity= Calculate_Cluster_Purity(c,cp,2);
        [BCubed,f_precision,f_recall]=Calculate_Cluster_BCubed_precision(c,cp);
        entropy=Calculate_cluster_entropy(c,cp);
        fm=Calculate_Cluster_F_measure(c,cp);
        comp_cop = cophenet(Zh,Ap_dis)
        
        
        details(k1,1)=k;
        details(k1,2)=rows;
        details(k1,3)=data_len;
        details(k1,6)=RS;
        details(k1,7)=purity;
        details(k1,8)=BCubed;
        details(k1,9)=f_precision;
        details(k1,10)=f_recall;
        details(k1,11)=entropy;
        details(k1,12)=fm;
        details(k1,13)=cop;
        details(k1,14)=cpcop;
        details(k1,15)=comp_cop;
        details(k1,16)=cp_count;   
        details(k1,17)=cutoff;
        % Import the file
        %     p=train_data(:,1);
        %     if(min(p)==0)
        %         p=p+1;
        %     end
        %
        %     RS= calculate_Cluster_Rand_statistic(c,p);
        %     purity= Calculate_Cluster_Purity(c,p,2);
        %     sim= Calculate_Cluster_Similarity(c,p);
        
        
        
    end
    
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

function slider_callback1(src,eventdata,arg1)
val = get(src,'Value');
set(arg1,'Position',[0 -val*5 1 5])
end